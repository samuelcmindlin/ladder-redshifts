"""Which single follow-up tuning settles an ambiguous real redshift? The tuning optimiser applied to a real source.

Takes a real spectrum (data/spectra/<SRC>.csv), runs the engine to get the candidate redshifts and their
posterior mass, then for every ALMA tuning in Bands 3-7 simulates the follow-up under each candidate as the
truth, added to the real data, and asks whether the engine then returns a confident (>= 0.95) and correct z.

Model (deliberately simple, stated so it can be argued with):
  - a tuning = two 3.75 GHz sidebands (four 1.875 GHz windows) centred in the receiver's IF range: LO -/+
    4.125-7.875 GHz in Bands 3, 4, 5 and 7 (IF 4-8 GHz) and LO -/+ 5.375-9.125 GHz in Band 6 (IF 4.5-10 GHz),
    for every LO that keeps both sidebands inside the band, less the frequencies the real data already cover;
  - per-channel noise in the new tuning = depth x the median noise of the windows measured at the source
    (depth 1 = the same noise per channel as the existing scan; depth 3 = 3x noisier, roughly a ninth of the
    integration time, ignoring the change of system temperature between bands);
  - predicted lines scaled in flux from the observed line (its matched-filter S/N x the noise where it sits) by
    the CO/[CI] ladder (ladder.sled_flux) for low excitation (Jc = 2.0, the least excited template, i.e. the
    hardest case for high-J lines) and typical excitation (Jc = 3.7), with [CI]/CO(1-0) = 0.15 (the middle of
    the robust prior); line width = --fwhm km/s;
  - a tuning can separate two candidates only if it covers a line that one of them predicts and the other does
    not (no line of the other within 2 x DV_MAX = 300 km/s). A line both predict - a shared harmonic, such as
    CO(9-8) at z = 2.93 and CO(6-5) at z = 1.62 - shows that one of them is right, not which: its position
    differs only by the ladder's anharmonicity, which the engine deliberately ignores, and its strength
    relative to the observed line is no safe guide. When a tuning leaves the true candidate unseparated from a
    rival, its success probability for that candidate is 0 and no simulation is run.
The figure shows, per LO, the posterior-weighted probability of a confident and correct answer (averaged over
noise draws) and the predicted lines of each hypothesis.

Writes results/redshift_followup_<SRC>.json and figures/redshift_followup_<SRC>.png.
Usage: python -m ladder_redshifts.followup [--src SPT0314-44] [--fwhm 400] [--nreal 10] [--step 1.0] [--plot-only]
(Noise draws are skipped where the noise-free follow-up leaves the true z below 0.5 posterior; P is then 0.)
"""
from .paths import FIG, RESULTS, DATA  # noqa: E402
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from .ladder import line_list, sled_flux, C_KMS, PIX_KMS, TOL_LN, DV_MAX  # noqa: E402
from .batch_real import read_spectrum, to_grid_auto, read_truth, PRIORS  # noqa: E402
from .realspec import run_engine  # noqa: E402

INK, INK2, GRID_C, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
COL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300']
BANDS = {'Band 3': (84.0, 116.0), 'Band 4': (125.0, 163.0), 'Band 5': (158.0, 211.0),
         'Band 6': (211.0, 275.0), 'Band 7': (275.0, 373.0)}          # sky frequency (GHz)
IF_RANGE = {'Band 3': (4.0, 8.0), 'Band 4': (4.0, 8.0), 'Band 5': (4.0, 8.0),
            'Band 6': (4.5, 10.0), 'Band 7': (4.0, 8.0)}             # intermediate frequency (GHz)
SIDEBAND_GHZ = 3.75                                                   # two 1.875 GHz windows per sideband
CONDITIONS = {'depth 1, low excitation': (1.0, 2.0), 'depth 3, low excitation': (3.0, 2.0),
              'depth 1, typical excitation': (1.0, 3.7)}
R_CI = 0.15


def sideband_offsets(bname):
    """(near, far) distance of a sideband's edges from the LO: SIDEBAND_GHZ wide, centred in the band's IF."""
    mid = 0.5 * sum(IF_RANGE[bname])
    return mid - SIDEBAND_GHZ / 2, mid + SIDEBAND_GHZ / 2


def lo_grid(bname, step=1.0):
    """LO frequencies on a `step` GHz grid for which both sidebands lie inside the band."""
    far = sideband_offsets(bname)[1]
    lo0 = np.ceil((BANDS[bname][0] + far) / step - 1e-9) * step
    return np.arange(lo0, BANDS[bname][1] - far + 1e-9, step)


def tuning(lo_ghz, bname, covered):
    """Windows of one tuning of band `bname` at LO `lo_ghz`, less the frequencies already `covered`."""
    band, (near, far) = BANDS[bname], sideband_offsets(bname)
    out = []
    for a, b in [(lo_ghz - far, lo_ghz - near), (lo_ghz + near, lo_ghz + far)]:
        a, b = max(a, band[0]), min(b, band[1])
        pieces = [(a, b)]
        for c0, c1 in covered:                                   # remove what the real data already cover
            nxt = []
            for p0, p1 in pieces:
                if c1 <= p0 or c0 >= p1:
                    nxt.append((p0, p1)); continue
                if c0 > p0:
                    nxt.append((p0, c0))
                if c1 < p1:
                    nxt.append((c1, p1))
            pieces = nxt
        out += [(p0, p1) for p0, p1 in pieces if p1 - p0 > 0.2]
    return out


def separates(wins, z_a, z_b, nus):
    """True if the windows cover a line that only one of the two redshifts predicts: a line of one with no line
    of the other within 2 x DV_MAX (each line may sit DV_MAX from its predicted position, so closer lines cannot
    be told apart by position)."""
    fa, fb = nus / (1 + z_a), nus / (1 + z_b)
    for f1, f2 in ((fa, fb), (fb, fa)):
        for f in f1:
            if any(a <= f <= b for a, b in wins) and np.min(np.abs(np.log(f / f2))) * C_KMS > 2 * DV_MAX:
                return True
    return False


def channels(a, b):
    n = int(np.log(b / a) / (PIX_KMS / C_KMS)) + 1
    return a * np.exp(np.arange(n) * PIX_KMS / C_KMS)


def synth_windows(wins, z_true, amp_obs, j_obs, Jc, sig_new, fwhm, rng):
    """New windows (same dict format as read_spectrum) with predicted lines + noise, in mJy.
    amp_obs: matched-filter S/N the observed line (index j_obs) would have at the new windows' noise sig_new,
    i.e. its flux / sig_new; the other lines follow from the template."""
    names, nus = line_list()
    flux = sled_flux(Jc, R_CI)
    rel = flux / flux[j_obs]
    s = fwhm / 2.3548 / PIX_KMS
    knorm = 1.0 / np.sqrt(np.sqrt(np.pi) * s)
    out = []
    for k, (a, b) in enumerate(wins):
        f = channels(a, b)
        y = rng.normal(0, 1, len(f)) if rng is not None else np.zeros(len(f))
        for nu, r in zip(nus, rel):
            fo = nu / (1 + z_true)
            if a - 0.05 * a < fo < b + 0.05 * b:
                v = np.log(f / fo) * C_KMS / PIX_KMS                 # offset in pixels
                y += amp_obs * r * knorm * np.exp(-0.5 * (v / s) ** 2)
        out.append(dict(spw=f'new{k}', freq=f, S=y * sig_new, sig=np.full(len(f), sig_new),
                        cont=np.zeros(len(f)), kind='measured'))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--src', default='SPT0314-44')
    ap.add_argument('--fwhm', type=float, default=400.0)
    ap.add_argument('--nreal', type=int, default=10)
    ap.add_argument('--step', type=float, default=1.0)
    ap.add_argument('--seed', type=int, default=8)
    ap.add_argument('--plot-only', action='store_true', help='re-plot from results/redshift_followup_<SRC>.json')
    a = ap.parse_args()
    spec = read_spectrum(DATA / 'spectra' / f'{a.src}.csv')
    if a.plot_only:
        o = json.loads((RESULTS / f'redshift_followup_{a.src}.json').read_text())
        hyp = [dict(z=z, line=ln, mass=m) for z, ln, m in o['current']['candidates']]
        covered = sorted((float(w['freq'].min()), float(w['freq'].max())) for w in spec)
        plot(a, o['scan'], hyp, o['predicted_lines'], covered, None)
        print(f'{a.src}: re-plotted the saved follow-up scan')
        for r in o['best'][:3]:
            print(f"  {r['band']} LO {r['lo_ghz']:6.1f} GHz  " + '  '.join(f"{k}: {r[k]['mc']:.2f}" for k in CONDITIONS))
        return
    t = read_truth().get(a.src, {})
    pname = t.get('prior') or 'spt'
    prior = PRIORS[pname]
    r0 = run_engine(*to_grid_auto(spec), prior=prior)
    cands = [(z, m) for z, m in r0['candidates'] if m >= 0.01]
    f_det, snr_det = max(r0['detections'], key=lambda d: d[1])
    names, nus = line_list()
    hyp = []
    for z, m in cands:
        j = int(np.argmin(np.abs(np.log(nus / (1 + z) / f_det))))
        hyp.append(dict(z=z, mass=m, line=names[j], j=j))
    print(f'{a.src}: current engine z = {r0["sled"][0]:.4f} (conf {r0["sled"][1]:.3f}); brightest line '
          f'{f_det:.3f} GHz at matched-filter S/N {snr_det:.1f}')
    for h in hyp:
        print(f'  candidate z = {h["z"]:.4f} ({h["line"]}), posterior mass {h["mass"]:.3f}')
    meas = [w for w in spec if w.get('kind', 'measured') == 'measured']
    sig_obs = float(np.median(np.concatenate([w['sig'] for w in meas])))
    # flux of the observed line in the engine's units: its matched-filter S/N x the noise where it sits (which
    # is not the median noise: this line lies where two tunings overlap). Predicted lines are scaled from it.
    grid0, _, sig0, _ = to_grid_auto(spec)
    sig_line = float(sig0[int(round(grid0.pix(f_det)))])
    amp = snr_det * sig_line
    covered = sorted((float(w['freq'].min()), float(w['freq'].max())) for w in spec)
    rng = np.random.default_rng(a.seed)
    rows = []
    for bname, band in BANDS.items():
        for lo in lo_grid(bname, a.step):
            wins = tuning(lo, bname, covered)
            if not wins:
                continue
            res = dict(band=bname, lo_ghz=round(float(lo), 2),
                       windows=[(round(float(x), 2), round(float(y), 2)) for x, y in wins])
            for ename, (depth, Jc) in CONDITIONS.items():
                sig_new = sig_obs * depth
                p_asimov, p_mc = 0.0, 0.0
                per = {}
                for h in hyp:
                    def outcome(rr):
                        okz = abs(np.log1p(rr['sled'][0]) - np.log1p(h['z'])) < TOL_LN
                        return float(okz and rr['sled'][1] >= 0.95)
                    if not all(separates(wins, h['z'], g['z'], nus) for g in hyp if g is not h):
                        per[f'{h["z"]:.4f}'] = dict(asimov=0.0, mc=0.0, unseparated=True)
                        continue                      # only lines shared with a rival: cannot settle it
                    new = synth_windows(wins, h['z'], amp / sig_new, h['j'], Jc, sig_new, a.fwhm, None)
                    ra = run_engine(*to_grid_auto(spec + new), prior=prior)
                    pa = outcome(ra)
                    pt = ra['post'][np.abs(ra['lnz'] - np.log1p(h['z'])) < 0.006].sum()
                    if pt < 0.5:                      # noise-free data leave the truth below 0.5: skip the draws
                        pm = 0.0
                    else:
                        hits = [outcome(run_engine(*to_grid_auto(spec + synth_windows(
                            wins, h['z'], amp / sig_new, h['j'], Jc, sig_new, a.fwhm, rng)), prior=prior))
                            for _ in range(a.nreal)]
                        pm = float(np.mean(hits))
                    per[f'{h["z"]:.4f}'] = dict(asimov=pa, mc=pm)
                    p_asimov += h['mass'] * pa; p_mc += h['mass'] * pm
                tot = sum(h['mass'] for h in hyp)
                res[ename] = dict(asimov=round(p_asimov / tot, 3), mc=round(p_mc / tot, 3), per_hypothesis=per)
            rows.append(res)
    ks = list(CONDITIONS)
    best = sorted(rows, key=lambda r: (-r[ks[1]]['mc'], -r[ks[0]]['mc'], -r[ks[2]]['mc'], r['lo_ghz']))[:8]
    print(f'\nbest single tunings (posterior-weighted P(confident & correct), {a.nreal} noise draws):')
    for r in best:
        print(f"  {r['band']} LO {r['lo_ghz']:6.1f} GHz  windows {r['windows']}  " +
              '  '.join(f"{k}: {r[k]['mc']:.2f}" for k in ks))
    for k in ks:
        n_dec = sum(r[k]['mc'] >= 0.9 for r in rows)
        print(f'  {k}: {n_dec} of {len(rows)} LO settings (step {a.step} GHz) reach P >= 0.9')
    # discriminating lines: predicted by one hypothesis, not within 1000 km/s of a line of another
    lines = []
    for h in hyp:
        for n, nu in zip(names, nus):
            fo = nu / (1 + h['z'])
            if 84 < fo < 373:
                shared = any(abs(np.log(fo * (1 + g['z']) / nu2)) < 1000 / C_KMS
                             for g in hyp if g is not h for nu2 in nus)
                lines.append(dict(z=round(h['z'], 4), line=n, freq_ghz=round(float(fo), 3), shared=bool(shared)))
    out = dict(source=a.src, prior=pname, conditions=CONDITIONS, fwhm_kms=a.fwhm, nreal=a.nreal, r_ci=R_CI,
               current=dict(z=round(r0['sled'][0], 4), confidence=round(r0['sled'][1], 3),
                            candidates=[(round(h['z'], 4), h['line'], round(h['mass'], 3)) for h in hyp]),
               observed_line=dict(freq_ghz=round(f_det, 3), snr=round(snr_det, 1)),
               predicted_lines=lines, best=best, scan=rows)
    (RESULTS / f'redshift_followup_{a.src}.json').write_text(json.dumps(out, indent=1))
    plot(a, rows, hyp, lines, covered, r0)


def plot(a, rows, hyp, lines, covered, r0):
    fig, ax = plt.subplots(figsize=(14, 5.2), facecolor=SURF)
    hpal = [COL[0], COL[1], COL[4], COL[3], COL[5]]            # COL[2] is the depth-3 curve
    hcol = {round(h['z'], 4): hpal[i % len(hpal)] for i, h in enumerate(hyp)}
    styles = {'depth 1, low excitation': dict(color=INK, lw=1.6),
              'depth 3, low excitation': dict(color=COL[2], lw=1.4, ls='--'),
              'depth 1, typical excitation': dict(color=INK2, lw=0.8, ls=':')}
    for bname, band in BANDS.items():
        ax.axvspan(band[0], band[1], color=GRID_C, alpha=0.35, lw=0)
        ax.text((band[0] + band[1]) / 2, 1.15, bname, ha='center', fontsize=8.5, color=INK2)
        rr = [r for r in rows if r['band'] == bname]
        for k, st in styles.items():
            if rr and k in rr[0]:
                ax.plot([r['lo_ghz'] for r in rr], [r[k]['mc'] for r in rr], **st)
    for c0, c1 in covered:
        ax.axvspan(c0, c1, ymin=0, ymax=0.035, color=INK2, lw=0)
    for ln in lines:
        c = '#a3a29d' if ln['shared'] else hcol[ln['z']]
        ax.plot([ln['freq_ghz']] * 2, [1.03, 1.10], color=c, lw=1.0 if ln['shared'] else 2.4)
    ax.set_xlim(84, 373); ax.set_ylim(-0.03, 1.2)
    ax.set_xlabel('LO frequency of the follow-up tuning (GHz). Ticks at the top mark predicted lines at their sky '
                  'frequency \u03bd; a tuning catches a line when LO = \u03bd \u00b1 4.1\u20137.9 GHz '
                  '(5.4\u20139.1 in Band 6), hence the paired plateaus.', fontsize=8.5, color=INK2)
    ax.set_ylabel('P(confident & correct),\nweighted by current posterior', fontsize=9, color=INK2)
    hz = ', '.join(f"z = {h['z']:.3f} ({h['line']}, {h['mass']:.2f})" for h in hyp)
    ax.set_title(f"{a.src}: which single follow-up tuning settles the redshift? Now: {hz}", fontsize=10,
                 color=INK, loc='left', pad=24)
    handles = [Line2D([], [], label='same noise per channel as the Band 3 scan, low excitation', **styles[
        'depth 1, low excitation']),
               Line2D([], [], label='3\u00d7 noisier (\u2248 1/9 of the time), low excitation', **styles[
                   'depth 3, low excitation']),
               Line2D([], [], label='same noise, typical excitation', **styles['depth 1, typical excitation'])]
    handles += [Patch(color=hcol[round(h['z'], 4)], label=f"line predicted only if z = {h['z']:.3f}") for h in hyp]
    handles += [Patch(color='#a3a29d', label='line predicted by more than one'),
                Patch(color=INK2, label='already observed (Band 3 scan)')]
    ax.legend(handles=handles, fontsize=7.8, frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.13),
              ncol=3)
    ax.set_facecolor(SURF)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.tick_params(colors=INK2, labelsize=8)
    fig.tight_layout()
    fig.savefig(FIG / f'redshift_followup_{a.src}.png', dpi=120, facecolor=SURF)


if __name__ == '__main__':
    main()
