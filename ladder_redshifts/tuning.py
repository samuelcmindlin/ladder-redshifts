"""Tuning optimiser: where should a spectral scan's tunings go to get the most confident, correct redshifts?

Receiver models (from C. De Breuck, "ALMA Band 2 as a redshift machine", ESO, Feb 2026):
  * Band 2 now   (67-116 GHz): one tuning = 4 contiguous spectral windows = 7.5 GHz block; 7 tunings cover the band.
  * Band 2 + WSU (2x bandwidth): one tuning = two 16 GHz sidebands (IF 2-18 GHz) around the LO, 4 GHz gap;
    2 tunings cover the band.
  * Band 3 now   (84-116 GHz, IF 4-8 GHz): one tuning = two 3.75 GHz sidebands at LO -/+ 4.125-7.875 GHz, as in
    the Band 3 scans of data/ (the classic 5-tuning scan covers 84.2-114.9 GHz).
Bakx 2026 (arXiv:2604.01089) studies the same question for Band 3 and for Band 2 with 8 and 14 GHz of bandwidth,
using the tool of Bakx & Dannerbauer 2022; the two studies assume different receivers and have not been compared.
Population: SPT-like redshifts (Reuter+20: median 3.9) and the atlas's SLED/[CI]/width distributions.
Time models: (i) each tuning gets the same time (depth fixed, total time grows with N);
             (ii) total time fixed, so per-line S/N falls as 1/sqrt(N).
Method: greedy selection of tuning positions on a 1 GHz LO grid, scoring each candidate set with the SLED
likelihood engine (robust prior, SPT-like redshift prior) on one simulated population; the chosen layouts
are then scored on an independent population (new seed), and those held-out fractions are the ones reported.
Writes results/redshift_tuning.json and figures/redshift_tuning.png.
Usage: python -m ladder_redshifts.tuning [N_sources]  (default 200; ~55 min on 2 cores)
"""
import sys as _sys
from .paths import AUDIO, FIG, RESULTS, DATA  # noqa: E402
import json
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .ladder import (Grid, Engine, sled_flux, simulate, matched_filter, templates_grid, prior_spt,  # noqa: E402
                    TOL_LN, CONF_WIN, CO, detection_index)

INK, INK2, GRID_C, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
COL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300']
BAND2 = (67.0, 116.0)
BAND3 = (84.0, 116.0)


def tuning_windows(kind, lo_ghz):
    """Frequency windows of one tuning. lo_ghz = lower edge of the block (Band 2 now) or LO frequency."""
    if kind == 'b2_now':
        return [(lo_ghz, lo_ghz + 7.5)]
    if kind == 'b2_wsu':
        return [(lo_ghz - 18.0, lo_ghz - 2.0), (lo_ghz + 2.0, lo_ghz + 18.0)]
    if kind == 'b3_now':                    # IF 4-8 GHz: two 1.875 GHz windows per sideband, centred in the IF
        return [(lo_ghz - 7.875, lo_ghz - 4.125), (lo_ghz + 4.125, lo_ghz + 7.875)]
    raise ValueError(kind)


def clip(windows, band):
    out = []
    for a, b in windows:
        a, b = max(a, band[0]), min(b, band[1])
        if b - a > 0.2:
            out.append((a, b))
    return out


STEP = 1.0                                   # LO / block-start grid (GHz)
HELD_OUT = {'b2_now_per_tuning': ('b2_now', BAND2, False), 'b2_now_fixed_total': ('b2_now', BAND2, True),
            'b2_wsu_per_tuning': ('b2_wsu', BAND2, False), 'b2_wsu_fixed_total': ('b2_wsu', BAND2, True),
            'b3_now_per_tuning': ('b3_now', BAND3, False)}


def candidates(kind):
    if kind == 'b2_now':
        return [round(x, 2) for x in np.arange(67.0, 108.51, STEP)]
    if kind == 'b2_wsu':
        return [round(x, 2) for x in np.arange(69.0, 114.01, STEP)]
    if kind == 'b3_now':
        return [round(x, 2) for x in np.arange(92.0, 108.01, STEP)]      # both sidebands inside the band


class Sim:
    """Simulated sources, stored as separate noiseless-signal and noise spectra, so that the per-line S/N can be
    rescaled (time model ii) without re-simulating: d(a) = a*d_sig + d_noise. Each candidate layout filters
    only its own observed pixels."""

    def __init__(self, n, rho, seed=11):
        rng = np.random.default_rng(seed)
        self.grid = Grid(62.0, 122.0)
        z = np.clip(rng.normal(3.9, 1.1, n), 1.5, 7.0)
        Jc = np.clip(np.exp(rng.normal(np.log(3.7), 0.35, n)), 1.8, 8)
        rci = np.exp(rng.normal(np.log(0.29), 0.18 * np.log(10), n))
        fwhm = np.clip(np.exp(rng.normal(np.log(450), 0.4, n)), 150, 1500)
        self.z = z
        self.dsig, self.dnoise = [], []

        class Zero:
            def normal(self, mu, s, k):
                return np.zeros(k)
        for i in range(n):
            f = sled_flux(Jc[i], rci[i], scatter_dex=0.1, rng=rng)
            amp = rho * f / f[3]
            self.dsig.append(simulate(self.grid, z[i], amp, fwhm[i], Zero()))
            self.dnoise.append(rng.normal(0, 1, self.grid.npix))
        self.dsig = np.array(self.dsig); self.dnoise = np.array(self.dnoise)
        self.T, self.lw = templates_grid(robust=True)

    def score(self, windows, scale=1.0):
        """Fraction confident (>=0.95) & correct, fraction confidently wrong, fraction with >=1 line >5 sigma."""
        if not windows:
            return 0.0, 0.0, 0.0
        cov = self.grid.mask(windows)
        e = Engine(self.grid, cov, self.T, logw=self.lw, zmin=0.5, zmax=8.0)
        cc = cw = det = 0
        for i in range(len(self.z)):
            Y = e.filter(scale * self.dsig[i] + self.dnoise[i])
            zs, c, _ = e.posterior(e.loglike(Y), prior=prior_spt)
            okz = abs(np.log1p(zs) - np.log1p(self.z[i])) < TOL_LN
            cc += okz and c >= 0.95; cw += (not okz) and c >= 0.95
            det += bool((Y[detection_index()][cov] > 5).any())
        n = len(self.z)
        return cc / n, cw / n, det / n


_SIM = None


def _eval(args):
    w, scale = args
    return _SIM.score(w, scale)


def greedy(sim, kind, nmax, band, fixed_total_time=False, log=None):
    global _SIM
    _SIM = sim
    import multiprocessing as mp
    chosen, path = [], []
    for k in range(1, nmax + 1):
        cands = [lo for lo in candidates(kind) if lo not in chosen]
        scale = 1.0 / np.sqrt(k) if fixed_total_time else 1.0
        jobs = [(clip(sum((tuning_windows(kind, x) for x in chosen + [lo]), []), band), scale) for lo in cands]
        with mp.get_context('fork').Pool(2) as pool:
            res = pool.map(_eval, jobs)
        best = None
        for lo, s in zip(cands, res):
            if best is None or s[0] > best[1][0] + 1e-9:
                best = (lo, s)
        chosen.append(best[0])
        path.append(dict(n_tunings=k, positions=list(chosen), confident_correct=best[1][0],
                         confidently_wrong=best[1][1], any_line_5sigma=best[1][2]))
        if log:
            log(f'  {kind} {"fixed-total" if fixed_total_time else "per-tuning"} k={k}: add {best[0]} GHz -> '
                f'confident&correct {best[1][0]:.3f}, confidently wrong {best[1][1]:.3f}, >=1 line {best[1][2]:.3f}')
    return path


def main():
    n = int(_sys.argv[1]) if len(_sys.argv) > 1 else 200
    t0 = time.time()
    log = lambda m: print(f'[{time.time() - t0:5.0f}s] {m}', flush=True)
    out = {'n_sources': n, 'population': 'SPT-like z ~ N(3.9, 1.1) clipped to 1.5-7', 'scenarios': {}}
    for label, rho in (('bright (CO(4-3) S/N 12 per tuning)', 12.0), ('faint (CO(4-3) S/N 5 per tuning)', 5.0)):
        log(f'simulating {n} sources, {label}')
        sim = Sim(n, rho)
        sc = {}
        # reference strategies
        refs = {'Band 3 5-tuning scan (84.2-114.9)': [(84.2, 114.9)],
                'Band 2 full, 7 tunings (67-116)': [BAND2],
                'Band 2 SV-like, 5 tunings (76.9-111)': [(76.9, 111.0)]}
        sc['references'] = {}
        for nm, w in refs.items():
            ntun = 5 if '5' in nm else 7
            s1 = sim.score(w, 1.0); s2 = sim.score(w, 1.0 / np.sqrt(ntun))
            sc['references'][nm] = dict(per_tuning=s1, fixed_total=s2, n_tunings=ntun)
            log(f'  ref {nm}: per-tuning {s1[0]:.3f} (wrong {s1[1]:.3f}); fixed total time {s2[0]:.3f}')
        sc['b2_now_per_tuning'] = greedy(sim, 'b2_now', 5, BAND2, False, log)
        sc['b2_now_fixed_total'] = greedy(sim, 'b2_now', 4, BAND2, True, log)
        sc['b2_wsu_per_tuning'] = greedy(sim, 'b2_wsu', 2, BAND2, False, log)
        sc['b2_wsu_fixed_total'] = greedy(sim, 'b2_wsu', 2, BAND2, True, log)
        sc['b3_now_per_tuning'] = greedy(sim, 'b3_now', 4, BAND3, False, log)
        # the layouts were chosen on `sim`; score them again on an independent population (new seed), so the
        # reported fractions are not inflated by the selection
        log(f'  held-out check: simulating {n} new sources')
        sim_ho = Sim(n, rho, seed=12)
        for key, (kind, band, fixed) in HELD_OUT.items():
            for q in sc[key]:
                w = clip(sum((tuning_windows(kind, x) for x in q['positions']), []), band)
                s = sim_ho.score(w, 1.0 / np.sqrt(q['n_tunings']) if fixed else 1.0)
                q['held_out'] = dict(confident_correct=s[0], confidently_wrong=s[1], any_line_5sigma=s[2])
            log(f"  held out {key}: " + ', '.join(f"{q['held_out']['confident_correct']:.3f}" for q in sc[key]))
        for nm, w in refs.items():
            r = sc['references'][nm]
            r['held_out'] = dict(per_tuning=sim_ho.score(w, 1.0), fixed_total=sim_ho.score(w, 1.0 / np.sqrt(r['n_tunings'])))
        out['scenarios'][label] = sc
    (RESULTS / 'redshift_tuning.json').write_text(json.dumps(out, indent=1, default=float))
    plot(out, FIG / 'redshift_tuning.png')
    log('done')


def plot(out, path):
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.8), facecolor=SURF)
    for j, (label, sc) in enumerate(out['scenarios'].items()):
        ax = axs[j]
        for key, nm, c, mk in (('b2_now_per_tuning', 'Band 2 now, each tuning same depth', COL[0], 'o'),
                               ('b2_now_fixed_total', 'Band 2 now, fixed total time', COL[1], 's'),
                               ('b2_wsu_per_tuning', 'Band 2 + WSU, same depth', COL[2], '^'),
                               ('b2_wsu_fixed_total', 'Band 2 + WSU, fixed total time', COL[3], 'D'),
                               ('b3_now_per_tuning', 'Band 3 now, same depth', COL[4], 'v')):
            p = sc[key]
            ax.plot([q['n_tunings'] for q in p], [q.get('held_out', q)['confident_correct'] for q in p], color=c, marker=mk, lw=2,
                    ms=6, label=nm)
        refs = sc['references']
        r7 = refs['Band 2 full, 7 tunings (67-116)']; r5 = refs['Band 3 5-tuning scan (84.2-114.9)']
        r7 = dict(r7.get('held_out', r7), n_tunings=7); r5 = dict(r5.get('held_out', r5), n_tunings=5)
        ax.plot([7], [r7['per_tuning'][0]], marker='o', ms=9, color=COL[0], mfc='none', mew=2, ls='none',
                label='Band 2 full band (7), same depth')
        ax.plot([7], [r7['fixed_total'][0]], marker='s', ms=9, color=COL[1], mfc='none', mew=2, ls='none',
                label='Band 2 full band (7), fixed total time')
        ax.plot([5], [r5['per_tuning'][0]], marker='v', ms=9, color=COL[4], mfc='none', mew=2, ls='none',
                label='Band 3 standard 5-tuning scan')
        ax.set_xticks(range(1, 8)); ax.set_xlim(0.6, 7.4)
        ax.set_ylim(0, 1); ax.set_xlabel('number of tunings (greedy-optimal positions)', color=INK2, fontsize=9)
        ax.set_ylabel('confident & correct redshifts (independent population)', color=INK2, fontsize=9)
        ax.set_title(f'{"AB"[j]}. {label}', fontsize=10, color=INK, loc='left')
        if j == 0:
            ax.legend(fontsize=7, frameon=False, loc='lower right', ncol=1)
        else:
            ax.text(0.98, 0.02, 'legend as in A', transform=ax.transAxes, ha='right', fontsize=7.5, color=INK2)
    ax = axs[2]
    sc = list(out['scenarios'].values())[0]
    for k, q in enumerate(sc['b2_now_per_tuning'][:4]):
        for lo in q['positions']:
            ax.add_patch(plt.Rectangle((lo, k + 0.65), 7.5, 0.7, color=COL[0], alpha=0.75, lw=0))
    for k, q in enumerate(sc['b2_wsu_per_tuning']):
        for lo in q['positions']:
            for a, b in clip(tuning_windows('b2_wsu', lo), BAND2):
                ax.add_patch(plt.Rectangle((a, k + 5.65), b - a, 0.7, color=COL[2], alpha=0.75, lw=0))
    for J in range(2, 8):
        for z in (2, 3, 4, 5, 6):
            pass
    ax.set_xlim(66, 117); ax.set_ylim(0.4, 7.6)
    ax.set_yticks([1, 2, 3, 4, 6, 7]); ax.set_yticklabels(['1 tuning', '2', '3', '4', 'WSU 1', 'WSU 2'], fontsize=8)
    ax.set_xlabel('frequency (GHz)', color=INK2, fontsize=9)
    ax.set_title('C. Greedy-optimal Band 2 tunings (bright case)', fontsize=10, color=INK, loc='left')
    for a in axs:
        a.set_facecolor(SURF)
        for sp in ('top', 'right'):
            a.spines[sp].set_visible(False)
        a.tick_params(colors=INK2, labelsize=8); a.grid(color=GRID_C, lw=0.6); a.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(path, dpi=130, facecolor=SURF)


if __name__ == '__main__':
    if '--plot-only' in _sys.argv:
        plot(json.loads((RESULTS / 'redshift_tuning.json').read_text()), FIG / 'redshift_tuning.png')
    else:
        main()
