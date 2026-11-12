"""Known octave traps, reconstructed: does the ladder likelihood flag them as ambiguous BEFORE the
decisive follow-up data, and resolve them after?

Scenarios use the published frequency coverage and line identifications; line strengths are
ILLUSTRATIVE where the papers (as fetched) do not give them. See docs/results.md.

  HerBS-108 (z-GAL, Cox et al. 2023, Sec. 3.1): lines at 97.75 and 146.6 GHz fit CO(2-1)/(3-2)
      or CO(4-3)/(6-5); z_phot = 2.4 sits between. NOEMA 3 mm main 75.874-106.868 GHz, 2 mm
      127.377-158.380 GHz; a line at 73.3 GHz in an extra 70.689-78.433 GHz setup (CO(3-2))
      gave z = 3.7166.
  SPT0551-50 (Strandet et al. 2016, Sec. 3.2): one Band 3 line, published as CO(3-2) at
      z = 2.1232 (with a false optical CIV "marker"); corrected to CO(4-3), z = 3.164, by APEX
      [CII] and CO(8-7).
  SPT0452-50 (Riechers 2025, ApJ 980, 59; arXiv:2401.05487): CO lines read as CO(1-0)+CO(3-2), z = 2.011;
      really CO(2-1)+CO(6-5), z = 5.016. CO(5-4) was missed through edge-channel flagging;
      archival [CI] windows at z = 2.011 were empty; [CII] (29 sigma) settled it.
Writes results/redshift_traps.json and figures/redshift_traps.png.
"""
import sys as _sys
from .paths import AUDIO, FIG, RESULTS, DATA  # noqa: E402
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .ladder import (Grid, Engine, line_list, simulate, matched_filter, templates_grid, TOL_LN,  # noqa: E402
                    CONF_WIN, prior_spt, prior_herschel)

INK, INK2, GRID_C, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'

# Windows are (lo GHz, hi GHz, relative noise). Line strengths are in the same flux units, so a
# line's S/N = strength / noise of its window. Values marked (assumed) are illustrative.
B3 = (84.2, 114.9, 1.0)
SCENARIOS = [
    dict(name='HerBS-108', prior=prior_herschel, prior_name='Herschel-like', z_true=3.7166, z_alt=(1 + 3.7166) / 2 - 1, fwhm=500, z_phot=2.4,
         note='z-GAL quotes z = 1.25 for the CO(2-1)/(3-2) reading; the line frequencies imply 1.358 = (1+3.7166)/2 - 1.',
         base={'CO(4-3)': 8.0},
         variants=[('low excitation, weak [CI]', {'CO(6-5)': 6.4, '[CI](1-0)': 0.8, 'CO(3-2)': 7.4}),
                   ('median excitation, weak [CI]', {'CO(6-5)': 10.7, '[CI](1-0)': 0.8, 'CO(3-2)': 5.9}),
                   ('high excitation, weak [CI]', {'CO(6-5)': 13.5, '[CI](1-0)': 0.8, 'CO(3-2)': 5.2}),
                   ('counterfactual: typical [CI]', {'CO(6-5)': 10.7, '[CI](1-0)': 6.1, 'CO(3-2)': 5.9})],
         steps=[('3 mm + 2 mm scans', [(75.874, 106.868, 1.0), (127.377, 158.380, 1.0)]),
                ('+ extra 3 mm setup', [(75.874, 106.868, 1.0), (127.377, 158.380, 1.0), (70.689, 78.433, 1.0)])]),
    dict(name='SPT0551-50', prior=prior_spt, prior_name='SPT-like', z_true=3.164, z_alt=(1 + 3.164) * 3 / 4 - 1, fwhm=450, z_phot=None,
         note='Published z = 2.1232 (CO(3-2) plus a false optical CIV marker); corrected z = 3.164.',
         base={'CO(4-3)': 6.0, '[CII]': 120.0},
         variants=[('low excitation', {'CO(8-7)': 0.54 * 6}), ('median excitation', {'CO(8-7)': 1.36 * 6}),
                   ('high excitation', {'CO(8-7)': 2.27 * 6})],
         steps=[('Band 3 scan (one line)', [B3]),
                ('+ 1 mm window', [B3, (217.0, 225.0, 2.0)]),
                ('+ [CII] window', [B3, (217.0, 225.0, 2.0), (452.0, 460.0, 3.0)])]),
    dict(name='SPT0452-50', prior=prior_spt, prior_name='SPT-like', z_true=5.016, z_alt=(1 + 5.016) / 2 - 1, fwhm=450, z_phot=None,
         note='Published z = 2.011 (CO(1-0)+CO(3-2)); corrected z = 5.016. Fluxes (Jy km/s) from Riechers 2025 '
              'Table 4 where given: CO(2-1) 0.96 +/- 0.12, CO(5-4) 2.80 +/- 0.67, [CII] 54.8 +/- 4.0; CO(6-5) '
              'not quoted in what we fetched (varied); archival [CI] window depth unknown (assumed 2.0).',
         base={'CO(2-1)': 0.96, 'CO(5-4)': 2.80, '[CII]': 54.8},
         variants=[('CO(6-5) = 2.0 Jy km/s', {'CO(6-5)': 2.0}), ('CO(6-5) = 3.1', {'CO(6-5)': 3.1}),
                   ('CO(6-5) = 5.0', {'CO(6-5)': 5.0})],
         steps=[('ATCA + Band 3, CO(5-4) flagged', [(37.3, 39.3, 0.12), (84.2, 95.4, 0.67), (96.2, 115.0, 0.67)]),
                ('+ empty [CI] windows for z=2.011', [(37.3, 39.3, 0.12), (84.2, 95.4, 0.67), (96.2, 115.0, 0.67),
                                                       (161.5, 165.5, 2.0), (266.8, 270.8, 2.0)]),
                ('+ re-reduced Band 3 (CO(5-4))', [(37.3, 39.3, 0.12), (84.2, 115.0, 0.67), (161.5, 165.5, 2.0),
                                                    (266.8, 270.8, 2.0)]),
                ('+ [CII] window', [(37.3, 39.3, 0.12), (84.2, 115.0, 0.67), (161.5, 165.5, 2.0),
                                    (266.8, 270.8, 2.0), (312.0, 320.0, 4.0)])]),
]


def mass_near(lnz, p, z):
    return float(p[np.abs(lnz - np.log1p(z)) < CONF_WIN].sum())


def run(nreal=150, seed=11, ROBUST=True):
    rng = np.random.default_rng(seed)
    grid = Grid(30.0, 480.0)
    names, nus = line_list(include_cii=True)
    T, lw = templates_grid(include_cii=True, robust=ROBUST)
    results, examples = [], []
    for sc in SCENARIOS:
        engines = []
        for lab, w in sc['steps']:
            sig = grid.noise(w)
            # SLED likelihood (noise-weighted) and the equal-weight baseline (equal weights on line S/N)
            engines.append((lab, Engine(grid, grid.mask(w), T, include_cii=True, zmax=12.0, sigma=sig, logw=lw),
                            Engine(grid, grid.mask(w), np.ones(len(names)), include_cii=True, zmax=12.0)))
        # noise per pixel for simulation: union of all windows of the last step (a superset)
        sig_sim = grid.noise(sc['steps'][-1][1])
        res_sc = dict(name=sc['name'], z_true=sc['z_true'], z_alt=sc['z_alt'], note=sc['note'], fwhm_kms=sc['fwhm'],
                      redshift_prior=sc['prior_name'],
                      nreal=nreal, variants=[])
        for vi, (vlab, over) in enumerate(sc['variants']):
            fl = dict(sc['base']); fl.update(over)
            amp = np.array([fl.get(n, 0.0) for n in names])
            stats = {lab: {k: [] for k in ('ok', 'conf', 'pt', 'pa', 'fok', 'fconf', 'ppz', 'pt_flatprior')} for lab, _, _ in engines}
            for r in range(nreal):
                d = simulate(grid, sc['z_true'], amp, sc['fwhm'], rng, include_cii=True, sigma=sig_sim)
                ex_r = []
                for lab, e, ef in engines:
                    Y = e.filter(d)                       # only the pixels observed at this step enter
                    L = e.loglike(Y)
                    zs, c, p = e.posterior(L, prior=sc['prior'])
                    st = stats[lab]
                    st['ok'].append(abs(np.log1p(zs) - np.log1p(sc['z_true'])) < TOL_LN); st['conf'].append(c)
                    st['pt'].append(mass_near(e.lnz, p, sc['z_true'])); st['pa'].append(mass_near(e.lnz, p, sc['z_alt']))
                    zf, cf, _ = ef.posterior(ef.loglike(Y), prior=sc['prior'])
                    st['pt_flatprior'].append(mass_near(e.lnz, e.posterior(L)[2], sc['z_true']))
                    st['fok'].append(abs(np.log1p(zf) - np.log1p(sc['z_true'])) < TOL_LN); st['fconf'].append(cf)
                    if sc['z_phot'] is not None:
                        st['ppz'].append(mass_near(e.lnz, e.posterior(L, z_phot=sc['z_phot'], prior=sc['prior'])[2], sc['z_true']))
                    ex_r.append((e.lnz, p))
                if vi == 1:                       # baseline (median) variant: average posterior over draws
                    if r == 0:
                        acc = [(lz, pp / nreal) for lz, pp in ex_r]
                    else:
                        acc = [(lz, a + pp / nreal) for (lz, a), (_, pp) in zip(acc, ex_r)]
            if vi == 1:
                examples.append(acc)
            rows = []
            for lab, _, _ in engines:
                st = {k: np.array(v, float) for k, v in stats[lab].items()}
                ok, conf = st['ok'].astype(bool), st['conf'] >= 0.95
                row = dict(step=lab, confident_correct=float((ok & conf).mean()),
                           confident_wrong=float((~ok & conf).mean()),
                           flagged_ambiguous=float((~conf).mean()),
                           median_p_true=float(np.median(st['pt'])), median_p_alt=float(np.median(st['pa'])),
                           median_p_true_uniform_prior=float(np.median(st['pt_flatprior'])),
                           flat_confident_wrong=float((~st['fok'].astype(bool) & (st['fconf'] >= 0.95)).mean()))
                if len(st['ppz']):
                    row['median_p_true_with_photz'] = float(np.median(st['ppz']))
                rows.append(row)
            res_sc['variants'].append(dict(variant=vlab, line_strengths=fl, steps=rows))
        results.append(res_sc)
    return results, examples


def clusters(lnz, p, floor=1e-3):
    """Group the posterior into candidate solutions: (z of peak, total mass)."""
    keep = p > floor * p.max()
    out, i = [], 0
    idx = np.where(keep)[0]
    if not len(idx):
        return out
    groups = np.split(idx, np.where(np.diff(idx) > int(2 * TOL_LN / (lnz[1] - lnz[0])))[0] + 1)
    for g in groups:
        lo, hi = max(g[0] - 20, 0), min(g[-1] + 20, len(p))
        out.append((float(np.exp(lnz[lo + np.argmax(p[lo:hi])]) - 1), float(p[lo:hi].sum())))
    return out


def plot(results, examples, path):
    nrow = len(results); ncol = max(len(r['variants'][0]['steps']) for r in results)
    fig, axs = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 2.9 * nrow), facecolor=SURF, squeeze=False)
    for i, (res, ex) in enumerate(zip(results, examples)):
        steps = res['variants'][0]['steps']
        for j in range(ncol):
            ax = axs[i, j]
            ax.set_facecolor(SURF)
            if j >= len(steps):
                ax.axis('off'); continue
            lnz, p = ex[j]
            for zc, m in clusters(lnz, p):
                if m < 0.005:
                    continue
                is_true = abs(np.log1p(zc) - np.log1p(res['z_true'])) < CONF_WIN
                is_alt = abs(np.log1p(zc) - np.log1p(res['z_alt'])) < CONF_WIN
                c = AQUA if is_true else (ORANGE if is_alt else BLUE)
                ax.bar(1 + zc, m, width=0.06 * (1 + zc), color=c, edgecolor=SURF, lw=0.5)
                ax.text(1 + zc, m + 0.03, f'{m:.2f}', ha='center', fontsize=7, color=INK2)
            ax.plot([1 + res['z_true']], [-0.06], marker='^', color=AQUA, ms=8, clip_on=False)
            ax.plot([1 + res['z_alt']], [-0.06], marker='^', color=ORANGE, ms=8, clip_on=False)
            ax.set_xscale('log'); ax.set_xlim(1.5, 13); ax.set_ylim(0, 1.15)
            ax.set_xticks([2, 3, 4, 6, 8, 12]); ax.set_xticklabels(['1', '2', '3', '5', '7', '11'])
            from matplotlib.ticker import NullFormatter
            ax.xaxis.set_minor_formatter(NullFormatter())
            st = steps[j]
            ax.set_title(f"{res['name']}: {st['step']}", fontsize=8.5, color=INK, loc='left')
            for sp in ('top', 'right'):
                ax.spines[sp].set_visible(False)
            ax.tick_params(colors=INK2, labelsize=7.5)
            ax.grid(axis='y', color=GRID_C, lw=0.6); ax.set_axisbelow(True)
            if j == 0:
                ax.set_ylabel('posterior mass', fontsize=8, color=INK2)
            ax.set_xlabel('redshift z (log 1+z axis)', fontsize=8, color=INK2)
    from matplotlib.patches import Patch
    axs[0, ncol - 1].legend(handles=[Patch(color=AQUA, label='mass on the true z (triangle marks it)'),
                                     Patch(color=ORANGE, label='mass on the octave-family alternative (triangle)'),
                                     Patch(color=BLUE, label='mass on other candidates')],
                            fontsize=8, frameon=False, loc='center')
    fig.suptitle('Known octave traps, reconstructed: posterior mass on each candidate redshift before and after the decisive '
                 'data (mean over noise draws; second listed variant of each case)', fontsize=11, color=INK, x=0.01, ha='left')
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=130, facecolor=SURF)


def save_examples(examples, path):
    arr = {}
    for i, ex in enumerate(examples):
        for j, (lz, pp) in enumerate(ex):
            arr[f'{i}|{j}|lnz'] = lz; arr[f'{i}|{j}|p'] = pp
    np.savez_compressed(path, **arr)


def load_examples(path):
    a = np.load(path)
    keys = sorted({tuple(map(int, k.split('|')[:2])) for k in a.files})
    out = {}
    for i, j in keys:
        out.setdefault(i, []).append((a[f'{i}|{j}|lnz'], a[f'{i}|{j}|p']))
    return [out[i] for i in sorted(out)]


if __name__ == '__main__':
    if '--plot-only' in _sys.argv:
        R = json.loads((RESULTS / 'redshift_traps.json').read_text())
        plot(R['robust_prior'], load_examples(RESULTS / 'redshift_traps_examples.npz'), FIG / 'redshift_traps.png')
        print('re-plotted figures/redshift_traps.png from saved results'); _sys.exit()
    n = int(_sys.argv[1]) if len(_sys.argv) > 1 else 150
    results, examples = run(n, ROBUST=True)
    save_examples(examples, RESULTS / 'redshift_traps_examples.npz')
    narrow, _ = run(n, ROBUST=False)
    (RESULTS / 'redshift_traps.json').write_text(json.dumps(
        {'robust_prior': results, 'narrow_prior_v0': narrow,
         'note': 'robust_prior adds a weak-[CI] template (L ratio 0.05, prior 0.2); narrow_prior_v0 was the first version.'},
        indent=1))
    plot(results, examples, FIG / 'redshift_traps.png')
    for title, R in (('ROBUST [CI] PRIOR (adopted)', results), ('NARROW [CI] PRIOR (first version, for the record)', narrow)):
        print(f"\n################ {title}")
        for res in R:
            print(f"\n== {res['name']}: true z = {res['z_true']}, octave-family alternative z = {res['z_alt']:.4f}")
            print(f"   {res['note']}  [redshift prior: {res['redshift_prior']}]")
            for v in res['variants']:
                print(f"  -- {v['variant']}")
                for st in v['steps']:
                    extra = f", with photo-z {st['median_p_true_with_photz']:.2f}" if 'median_p_true_with_photz' in st else ''
                    print(f"     {st['step']:36s} ambiguous {st['flagged_ambiguous']:.2f} | confident&correct "
                          f"{st['confident_correct']:.2f} | confidently WRONG {st['confident_wrong']:.2f} "
                          f"(flat sum {st['flat_confident_wrong']:.2f}) | P(true) {st['median_p_true']:.2f}{extra}, "
                          f"P(alt) {st['median_p_alt']:.2f} [P(true), flat prior to z=12: {st['median_p_true_uniform_prior']:.2f}]")
