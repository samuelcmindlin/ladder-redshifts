"""Octave-error atlas: how often each redshift estimator is right, wrong, or confidently wrong,
for ALMA-style frequency setups, as a function of source brightness and redshift.

Synthetic population (illustrative parameters; see ladder.py and docs/method.md):
  z ~ U(1, 7); rho = S/N the source's CO(4-3) line would have, log-uniform 2-64;
  CO SLED r_J1 = exp(-(J-1)/Jc), Jc lognormal (median 3.7, 0.35 in ln), 0.1 dex per-line scatter;
  [CI](1-0) L' ratio lognormal (median 0.29, 0.18 dex); FWHM lognormal (median 450 km/s);
  photometric redshift with ln(1+z) scatter 0.20.
The default seed (20261112) was first used after the changes of 12 November 2026 (velocity-offset tolerance),
so these sources played no part in developing the method (development used seeds 7 and 20261102).
Writes results/redshift_atlas.json and figures/redshift_atlas.png.
Usage: python -m ladder_redshifts.atlas [N_sources]   (default 2500; quick mode uses 600)
"""
import sys as _sys
from .paths import AUDIO, FIG, RESULTS, DATA  # noqa: E402
import json
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .ladder import (Grid, Engine, line_list, sled_flux, simulate, matched_filter, templates_grid,  # noqa: E402
                    detect, naive_spacing, sudoku_rule, classify, TOL_LN, detection_index)

SETUPS = {
    'B3 single tuning': [(89.5, 93.25), (101.5, 105.25)],
    'B3 scan': [(84.2, 114.9)],
    'B4 scan (2 mm)': [(125.0, 163.0)],
    'B2 full (WSU)': [(67.0, 116.0)],
    'B3 + B4 scans': [(84.2, 114.9), (125.0, 163.0)],
    'B2 + B4': [(67.0, 116.0), (125.0, 163.0)],
}
EST = ['naive spacing', 'sudoku rule', 'flat sum', 'sled', 'sled + photo-z']
RHO_BINS = [2, 4, 8, 16, 32, 64]
CONF = 0.95

COL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300']
MK = ['o', 's', '^', 'D', 'v', 'P']
INK, INK2, GRID_C, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'


def population(n, rng):
    z = rng.uniform(1, 7, n)
    rho = np.exp(rng.uniform(np.log(2), np.log(64), n))
    Jc = np.clip(np.exp(rng.normal(np.log(3.7), 0.35, n)), 1.8, 8)
    rci = np.exp(rng.normal(np.log(0.29), 0.18 * np.log(10), n))
    fwhm = np.clip(np.exp(rng.normal(np.log(450), 0.4, n)), 150, 1500)
    zph = (1 + z) * np.exp(rng.normal(0, 0.20, n)) - 1
    return z, rho, Jc, rci, fwhm, zph


def pop_weight(z, kind):
    if kind == 'SPT-like':        # Reuter+20: 1.9 < z < 6.9, median 3.9
        w = np.exp(-0.5 * ((z - 3.9) / 1.1) ** 2)
    else:                         # AS2UDS-like: median 2.61, 68% range 1.8-3.4 (Dudzevicute+20)
        w = np.exp(-0.5 * (np.log((1 + z) / 3.61) / 0.226) ** 2) / (1 + z)
    return w / w.sum()


def run(n, seed=20261112):
    rng = np.random.default_rng(seed)
    grid = Grid(66.0, 164.0)
    covers = {k: grid.mask(v) for k, v in SETUPS.items()}
    T, lw = templates_grid()
    eng = {k: (Engine(grid, c, T, logw=lw), Engine(grid, c, np.ones(T.shape[1]))) for k, c in covers.items()}
    z, rho, Jc, rci, fwhm, zph = population(n, rng)
    out = {k: {e: {'z': np.full(n, np.nan), 'conf': np.zeros(n)} for e in EST} for k in SETUPS}
    nlines = {k: np.zeros(n, int) for k in SETUPS}
    names, nus = line_list()
    t0 = time.time()
    for i in range(n):
        f = sled_flux(Jc[i], rci[i], scatter_dex=0.1, rng=rng)
        snr = rho[i] * f / f[3]
        d = simulate(grid, z[i], snr, fwhm[i], rng)
        for k, cov in covers.items():
            e_sled, e_flat = eng[k]
            Y = e_sled.filter(d)                            # only this setup's observed pixels enter
            nlines[k][i] = sum(1 for nu in nus if cov[int(round(min(max(grid.pix(nu / (1 + z[i])), 0), grid.npix - 1)))])
            L = e_sled.loglike(Y)
            for name, zp in (('sled', None), ('sled + photo-z', zph[i])):
                zs, c, _ = e_sled.posterior(L, z_phot=zp)
                out[k][name]['z'][i], out[k][name]['conf'][i] = zs, c
            zs, c, _ = e_flat.posterior(e_flat.loglike(Y))
            out[k]['flat sum']['z'][i], out[k]['flat sum']['conf'][i] = zs, c
            yd = Y[detection_index()]
            dets = detect(yd, cov)
            zn = naive_spacing(dets, grid)
            if zn is not None:
                out[k]['naive spacing']['z'][i], out[k]['naive spacing']['conf'][i] = zn, 1.0
            zr, _ = sudoku_rule(dets, grid, cov, yd)
            if zr is not None:
                out[k]['sudoku rule']['z'][i], out[k]['sudoku rule']['conf'][i] = zr, 1.0
        if (i + 1) % 250 == 0:
            print(f'  {i + 1}/{n} sources  ({time.time() - t0:.0f} s)', flush=True)
    return dict(z=z, rho=rho, Jc=Jc, fwhm=fwhm, out=out, nlines=nlines, seed=seed)


def summarize(R):
    z, rho, out = R['z'], R['rho'], R['out']
    S = {'setups': SETUPS, 'estimators': EST, 'rho_bins': RHO_BINS, 'confidence_threshold': CONF,
         'n_sources': int(len(z)), 'seed': R.get('seed'), 'by_rho': {}, 'by_z': {}, 'calibration': {},
         'reliability': {}, 'taxonomy': {}, 'population': {}}
    for k in SETUPS:
        S['by_rho'][k], S['by_z'][k], S['taxonomy'][k], S['population'][k] = {}, {}, {}, {}
        for e in EST:
            zs, c = out[k][e]['z'], out[k][e]['conf']
            ok = np.abs(np.log1p(np.nan_to_num(zs, nan=-0.5)) - np.log1p(z)) < TOL_LN
            conf = c >= CONF
            rows = []
            for lo, hi in zip(RHO_BINS[:-1], RHO_BINS[1:]):
                m = (rho >= lo) & (rho < hi)
                rows.append(dict(rho=[lo, hi], n=int(m.sum()),
                                 confident_correct=float((ok & conf & m).sum() / m.sum()),
                                 confident_wrong=float((~ok & conf & m).sum() / m.sum()),
                                 map_correct=float((ok & m).sum() / m.sum())))
            S['by_rho'][k][e] = rows
            mz = rho >= 8
            zrows = []
            for lo in np.arange(1, 7, 0.5):
                m = mz & (z >= lo) & (z < lo + 0.5)
                zrows.append(dict(z=[float(lo), float(lo + 0.5)], n=int(m.sum()),
                                  confident_correct=float((ok & conf & m).sum() / max(m.sum(), 1)),
                                  confident_wrong=float((~ok & conf & m).sum() / max(m.sum(), 1))))
            S['by_z'][k][e] = zrows
            labs = {}
            for zi, zt, cc in zip(zs, z, c):
                if np.isnan(zi):
                    continue
                lab = classify(float(zi), float(zt))
                if lab != 'correct':
                    key = lab + (' (confident)' if cc >= CONF else '')
                    labs[key] = labs.get(key, 0) + 1
            S['taxonomy'][k][e] = dict(sorted(labs.items(), key=lambda t: -t[1]))
            pops = {}
            for kind in ('SPT-like', 'AS2UDS-like'):
                w = pop_weight(z, kind)
                pr = {}
                for lo, hi in ((4, 8), (8, 16), (16, 32)):
                    m = (rho >= lo) & (rho < hi)
                    ww = w * m; ww = ww / ww.sum()
                    pr[f'{lo}-{hi}'] = dict(confident_correct=float((ww * (ok & conf)).sum()),
                                            confident_wrong=float((ww * (~ok & conf)).sum()))
                pops[kind] = pr
            S['population'][k][e] = pops
            if e in ('flat sum', 'sled', 'sled + photo-z'):
                edges = [0, 0.5, 0.8, 0.9, 0.95, 0.99, 1.0001]
                cal = []
                for lo, hi in zip(edges[:-1], edges[1:]):
                    m = (c >= lo) & (c < hi)
                    if m.sum() >= 10:
                        cal.append(dict(conf=[lo, min(hi, 1.0)], n=int(m.sum()), mean_conf=float(c[m].mean()),
                                        accuracy=float(ok[m].mean())))
                S['calibration'].setdefault(k, {})[e] = cal
                # how often a confident answer is right (all S/N together); compare with the stated confidence
                S['reliability'].setdefault(k, {})[e] = {
                    f'>={thr}': dict(n=int((c >= thr).sum()), fraction_correct=float(ok[c >= thr].mean()) if (c >= thr).any() else None,
                                     mean_stated=float(c[c >= thr].mean()) if (c >= thr).any() else None,
                                     n_wrong=int((~ok & (c >= thr)).sum()))
                    for thr in (0.95, 0.99)}
        S['by_z'][k]['mean_lines_in_band'] = [float(R['nlines'][k][(z >= lo) & (z < lo + 0.5)].mean())
                                              for lo in np.arange(1, 7, 0.5)]
    return S


def style(ax):
    ax.set_facecolor(SURF)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(color=GRID_C, lw=0.6); ax.set_axisbelow(True)


def plot(S, R, path):
    rc = np.sqrt(np.array(RHO_BINS[:-1]) * np.array(RHO_BINS[1:]))
    fig, axs = plt.subplots(2, 3, figsize=(15, 8.6), facecolor=SURF)
    ax = axs[0, 0]
    for j, k in enumerate(SETUPS):
        y = [r['confident_correct'] for r in S['by_rho'][k]['sled']]
        ax.plot(rc, y, color=COL[j], marker=MK[j], ms=6, lw=2, label=k)
    ax.set_xscale('log'); ax.set_xticks(rc); ax.set_xticklabels([f'{a}-{b}' for a, b in zip(RHO_BINS[:-1], RHO_BINS[1:])])
    ax.set_ylim(0, 1); ax.set_xlabel('source brightness: S/N its CO(4-3) line would have', color=INK2, fontsize=9)
    ax.set_ylabel('fraction with a confident, correct z', color=INK2, fontsize=9)
    ax.set_title('A. Setups (SLED likelihood, uniform z 1-7)', color=INK, fontsize=10, loc='left')
    ax.legend(fontsize=7.5, frameon=False, loc='lower right', bbox_to_anchor=(1.0, 0.1))

    for col, (key, title) in enumerate((('confident_correct', 'B. Estimators, Band 2 full (WSU): confident & correct'),
                                        ('confident_wrong', 'C. Estimators, Band 2 full (WSU): confidently WRONG'))):
        ax = axs[0, col + 1]
        for j, e in enumerate(EST):
            y = np.array([r[key] for r in S['by_rho']['B2 full (WSU)'][e]])
            if col == 1 and e == 'naive spacing':
                ax.plot(rc, np.minimum(y, 0.075), color=COL[j], marker=MK[j], ms=6, lw=2, ls='--', label=e + ' (off scale)')
                ax.text(rc[-1], 0.072, f'naive spacing: {y[1]:.0%}-{y.max():.0%}\n(off scale)', fontsize=7.5,
                        color=INK2, ha='right', va='top')
            else:
                ax.plot(rc, y, color=COL[j], marker=MK[j], ms=6, lw=2, label=e)
        ax.set_xscale('log'); ax.set_xticks(rc); ax.set_xticklabels([f'{a}-{b}' for a, b in zip(RHO_BINS[:-1], RHO_BINS[1:])])
        ax.set_xlabel('S/N of CO(4-3)-equivalent', color=INK2, fontsize=9)
        ax.set_title(title, color=INK, fontsize=10, loc='left')
        ax.set_ylim(0, 1 if col == 0 else 0.08)
        ax.legend(fontsize=7.5, frameon=False, loc='upper left' if col == 0 else 'center left')

    ax = axs[1, 0]
    z, out = R['z'], R['out']
    for j, e in enumerate(('flat sum', 'sled')):
        zs, c = out['B3 scan'][e]['z'], out['B3 scan'][e]['conf']
        d = np.log1p(zs) - np.log1p(z)
        bad = np.abs(d) >= TOL_LN
        ax.hist(np.exp(d[bad]), bins=np.exp(np.linspace(np.log(0.2), np.log(5), 161)), color=COL[j + 2],
                alpha=0.75, label=f'{e}: {bad.mean():.0%} wrong best guess\n(incl. sources with no line in band)', histtype='stepfilled' if j == 0 else 'step', lw=1.5)
    for r, lab in ((2, '2'), (1.5, '3/2'), (4 / 3, '4/3'), (0.5, '1/2'), (2 / 3, '2/3'), (0.75, '3/4'), (3, '3'), (1 / 3, '1/3')):
        ax.axvline(r, color=INK2, lw=0.6, ls=':')
        ax.text(r, ax.get_ylim()[1] * 0.97, lab, fontsize=7, color=INK2, ha='center', va='top')
    ax.set_xscale('log'); ax.set_xlim(0.2, 5)
    ax.set_xticks([0.25, 0.5, 1, 2, 4]); ax.set_xticklabels(['1/4', '1/2', '1', '2', '4'])
    ax.set_xlabel('(1 + z_est) / (1 + z_true) for wrong answers (Band 3 scan)', color=INK2, fontsize=9)
    ax.set_title('D. Wrong answers land on harmonic ratios', color=INK, fontsize=10, loc='left')
    ax.legend(fontsize=7.5, frameon=False, loc='center left')

    ax = axs[1, 1]
    ax.plot([0.4, 1], [0.4, 1], color=INK2, lw=1, ls='--')
    for j, (k, e) in enumerate((('B3 scan', 'flat sum'), ('B3 scan', 'sled'), ('B2 full (WSU)', 'sled'),
                                ('B3 scan', 'sled + photo-z'))):
        cal = S['calibration'][k][e]
        ax.plot([q['mean_conf'] for q in cal], [q['accuracy'] for q in cal], color=COL[[2, 3, 1, 4][j]],
                marker=MK[j], ms=6, lw=2, label=f'{e} - {k}')
    ax.set_xlim(0.4, 1.01); ax.set_ylim(0.0, 1.01)
    ax.set_xlabel('stated confidence (posterior mass on the best solution)', color=INK2, fontsize=9)
    ax.set_ylabel('actual fraction correct', color=INK2, fontsize=9)
    ax.set_title('E. Calibration: does 95% mean 95%?', color=INK, fontsize=10, loc='left')
    ax.legend(fontsize=7.5, frameon=False)

    ax = axs[1, 2]
    zc = np.arange(1, 7, 0.5) + 0.25
    for j, k in enumerate(SETUPS):
        y = [r['confident_correct'] for r in S['by_z'][k]['sled']]
        ax.plot(zc, y, color=COL[j], marker=MK[j], ms=5, lw=2, label=k)
    ax.set_ylim(0, 1); ax.set_xlabel('true redshift', color=INK2, fontsize=9)
    ax.set_ylabel('confident & correct (S/N >= 8)', color=INK2, fontsize=9)
    ax.set_title('F. Where each setup works (colours as in A)', color=INK, fontsize=10, loc='left')
    from matplotlib.ticker import NullFormatter
    for a in axs.flat:
        style(a)
        a.xaxis.set_minor_formatter(NullFormatter())
    fig.suptitle('Redshift from a CO ladder as missing-fundamental pitch estimation: synthetic octave-error atlas',
                 color=INK, fontsize=12, x=0.01, ha='left')
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=130, facecolor=SURF)


def save_raw(R, path):
    arr = dict(z=R['z'], rho=R['rho'], Jc=R['Jc'], fwhm=R['fwhm'], seed=np.array(R.get('seed', -1)))
    for k in SETUPS:
        arr[f'nlines|{k}'] = R['nlines'][k]
        for e in EST:
            arr[f'z|{k}|{e}'] = R['out'][k][e]['z']; arr[f'conf|{k}|{e}'] = R['out'][k][e]['conf']
    np.savez_compressed(path, **arr)


def load_raw(path):
    a = np.load(path)
    R = dict(z=a['z'], rho=a['rho'], Jc=a['Jc'], fwhm=a['fwhm'], out={}, nlines={},
             seed=int(a['seed']) if 'seed' in a.files else None)
    for k in SETUPS:
        R['nlines'][k] = a[f'nlines|{k}']
        R['out'][k] = {e: {'z': a[f'z|{k}|{e}'], 'conf': a[f'conf|{k}|{e}']} for e in EST}
    return R


if __name__ == '__main__':
    if '--plot-only' in _sys.argv:
        R = load_raw(RESULTS / 'redshift_atlas_raw.npz')
        plot(summarize(R), R, FIG / 'redshift_atlas.png'); print('re-plotted figures/redshift_atlas.png from saved results'); _sys.exit()
    n = int(_sys.argv[1]) if len(_sys.argv) > 1 else 2500
    R = run(n)
    save_raw(R, RESULTS / 'redshift_atlas_raw.npz')
    S = summarize(R)
    (RESULTS / 'redshift_atlas.json').write_text(json.dumps(S, indent=1))
    plot(S, R, FIG / 'redshift_atlas.png')
    print(f'\n{n} synthetic sources; confidence threshold {CONF}')
    for k in SETUPS:
        print(f'\n== {k}')
        for e in EST:
            rows = S['by_rho'][k][e]
            cc = ' '.join(f"{r['confident_correct']:.2f}" for r in rows)
            cw = ' '.join(f"{r['confident_wrong']:.3f}" for r in rows)
            print(f'  {e:15s} confident-correct by S/N bin: {cc}   confidently-wrong: {cw}')
    print('\nError taxonomy, Band 3 scan:')
    for e in EST:
        print(f'  {e:15s}', dict(list(S["taxonomy"]["B3 scan"][e].items())[:6]))
    print('\nCalibration (mean confidence -> accuracy):')
    for k in ('B3 scan', 'B2 full (WSU)'):
        for e in ('flat sum', 'sled', 'sled + photo-z'):
            print(f'  {k:14s} {e:15s}', [(round(q['mean_conf'], 2), round(q['accuracy'], 2), q['n']) for q in S['calibration'][k][e]])
    print('\nSPT-like population, S/N 8-16 (confident-correct / confidently-wrong):')
    for k in SETUPS:
        print(f'  {k:16s}', {e: (round(S['population'][k][e]['SPT-like']['8-16']['confident_correct'], 2),
                                 round(S['population'][k][e]['SPT-like']['8-16']['confident_wrong'], 3)) for e in EST})
