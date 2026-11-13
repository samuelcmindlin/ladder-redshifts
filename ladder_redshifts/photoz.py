"""Each source's published photometric redshift as its prior: the information the survey team had.

Eleven of the eighteen SPT sources of 2016.1.00672.S show one line in the 3 mm scan. With the survey's redshift
distribution as the only prior, the engine ranks the published redshift first for ten of them but is confident
of none. The survey team had more: a photometric redshift for every source, from the dust continuum (Reuter et
al. 2020, ApJ 902, 78, Table F.1; copied for our sources into data/photoz.csv by scripts/photoz_calibration.py).
Here each source's photometric redshift is the prior, as a Gaussian in z, in three versions:

  photoz        z_phot +/- error as printed in the 2020 paper (the main case);
  photoz_2024   the collaboration's 2024 revision of the same table, refitted after the redshift of SPT0452-50
                was corrected;
  photoz_recal  the 2020 values shifted and widened by the amounts that make them match the spectroscopic
                redshifts of the 60 SPT sources that are not in this sample (RECAL below; the 2020 values run
                high by about half their error bar).

Four ways of using a photometric redshift are run on every real spectrum that has one:
  - the engine, with the Gaussian as its prior;
  - the guided rule (ladder.guided_rule) with the same Gaussian: it answers when the lines single out one
    solution, or when one remaining solution holds 95% of the prior weight;
  - a window rule: the guided rule with a top-hat prior of half-width 0.13 (1 + z_phot), the photometric
    accuracy the High-z Sudoku code assumes by default. It answers when one solution lies inside the window;
  - the template cross-correlation multiplied by the Gaussian (ladder.xcorr_redshift), after Yun et al. 2015.
Both rules are run twice: counting a missing [CI] line against a redshift (as in the rest of the repository),
and counting CO lines only.

The noise-added test (degrade.py) runs the same three versions as priors on the seven multi-line spectra; its
photometric-redshift columns are drawn here.

Writes results/redshift_photoz.json, figures/redshift_photoz.png and, from results/redshift_degradation.json,
figures/redshift_photoz_noise.png.
Usage: python -m ladder_redshifts.photoz [--plot-only]     (a few seconds)
"""
from .paths import FIG, RESULTS, DATA  # noqa: E402
import argparse
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .ladder import guided_rule, xcorr_redshift, _rule_solutions, prior_spt, CONF_WIN, TOL_LN  # noqa: E402
from .batch_real import read_spectrum, to_grid_auto, read_truth  # noqa: E402
from .realspec import engine_likelihood, engine_posterior  # noqa: E402

INK, INK2, GRID_C, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
COL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300']
CONF = 0.95
WINDOW = 0.13                # half-width of the window rule in (1 + z_phot): the Sudoku code's dzUncertainty default
# Shift (in units of the published error) and widening that match the 2020 photometric redshifts to the
# spectroscopic ones for the 60 SPT sources outside this sample: results/redshift_photoz_calibration.json.
RECAL = dict(shift=-0.479, width=1.101)
STATES = {'photoz': 'photometric redshift, as published (2020)',
          'photoz_2024': 'photometric redshift, 2024 revision',
          'photoz_recal': 'photometric redshift, 2020 recalibrated'}


def read_photoz():
    """{source: (z_phot, error, z_phot_2024, error_2024)} from data/photoz.csv ({} if the file is absent)."""
    path = DATA / 'photoz.csv'
    if not path.exists():
        return {}
    with open(path, newline='') as fh:
        return {r['source']: tuple(float(r[k]) for k in ('z_phot', 'z_phot_err', 'z_phot_2024', 'z_phot_2024_err'))
                for r in csv.DictReader(fh)}


def gaussian(centre, width):
    """Redshift density (a function of z) for a photometric redshift centre +/- width."""
    return lambda z: np.exp(-0.5 * ((np.asarray(z, float) - centre) / width) ** 2)


def window(centre, frac=WINDOW):
    """Top-hat redshift density of half-width frac * (1 + centre)."""
    return lambda z: (np.abs(np.asarray(z, float) - centre) <= frac * (1 + centre)).astype(float)


def photoz_states(src, table=None):
    """{state: (centre, width)} for one source; {} if it has no published photometric redshift."""
    table = read_photoz() if table is None else table
    if src not in table:
        return {}
    zp, ep, zr, er = table[src]
    return {'photoz': (zp, ep), 'photoz_2024': (zr, er),
            'photoz_recal': (zp + RECAL['shift'] * ep, RECAL['width'] * ep)}


def photoz_priors(src, table=None):
    """{state: redshift density} for one source."""
    return {k: gaussian(*v) for k, v in photoz_states(src, table).items()}


def _outcome(z, z_true):
    if z is None:
        return 'none'
    return 'correct' if abs(np.log1p(z) - np.log1p(z_true)) < TOL_LN else 'wrong'


def _answer(z, z_true):
    return dict(z=None if z is None else round(float(z), 4), outcome=_outcome(z, z_true))


def analyse():
    """Every method on every real spectrum that has both a published and a photometric redshift."""
    table, rows = read_photoz(), []
    for src, t in sorted(read_truth().items()):
        path = DATA / 'spectra' / f'{src}.csv'
        if not path.exists() or src not in table:
            continue
        z_true = float(t['z_published'])
        g = to_grid_auto(read_spectrum(path))
        grid, d, _, cov = g
        like = engine_likelihood(*g)
        yd, dets = like['yd'], like['dets']
        near = np.abs(like['e'].lnz - np.log1p(z_true)) < CONF_WIN
        row = dict(source=src, z_published=z_true, lines='one' if '(1 line)' in t['notes'] else 'several',
                   detections=len(dets), reference={}, states={})
        for name, prior in (('uniform', None), ('spt', prior_spt)):
            r = engine_posterior(like, prior)
            row['reference'][name] = dict(z=round(r['sled'][0], 4), confidence=round(r['sled'][1], 3),
                                          p_published=round(float(r['post'][near].sum()), 4))
            if prior is None:
                row['candidates_no_prior'] = [(round(zc, 4), round(pc, 3)) for zc, pc in r['candidates'] if pc >= 0.02]
        sols = {ci: [round(float(np.exp(s) - 1), 4) for s in _rule_solutions(dets, grid, cov, yd, ci_prunes=ci)]
                if dets else [] for ci in (True, False)}
        row['rule_solutions'] = dict(ci_counted=sols[True], co_only=sols[False])
        for state, (centre, width) in photoz_states(src, table).items():
            prior = gaussian(centre, width)
            r = engine_posterior(like, prior)
            ze, ce = r['sled']
            q = dict(centre=round(centre, 3), width=round(width, 3), pull=round((z_true - centre) / width, 2),
                     engine=dict(z=round(ze, 4), confidence=round(ce, 3), outcome=_outcome(ze, z_true),
                                 p_published=round(float(r['post'][near].sum()), 4),
                                 candidates=[(round(zc, 4), round(pc, 3)) for zc, pc in r['candidates'][:4]]),
                     xcorr=_answer(xcorr_redshift(d, grid, cov, prior=prior)[0], z_true))
            for ci, key in ((True, 'ci_counted'), (False, 'co_only')):
                q[f'guided_{key}'] = _answer(guided_rule(dets, grid, cov, yd, prior=prior, ci_prunes=ci)[0], z_true)
                if state != 'photoz_recal':                # the window is centred on a published value
                    q[f'window_{key}'] = _answer(guided_rule(dets, grid, cov, yd, prior=window(centre),
                                                             ci_prunes=ci)[0], z_true)
                w = np.array([float(prior(s)) * (1 + s) for s in sols[ci]])
                pick = sols[ci][int(np.argmax(w))] if len(w) and w.sum() > 0 else None
                q[f'pick_{key}'] = dict(_answer(pick, z_true),
                                        share=None if pick is None else round(float(w.max() / w.sum()), 3))
            row['states'][state] = q
        rows.append(row)
    return rows


def summarise(rows):
    """Counts per group of sources (one line in the scan, or several), prior version and method."""
    out = {}
    for group in ('one', 'several'):
        rr = [r for r in rows if r['lines'] == group]
        g = dict(n=len(rr), sources=[r['source'] for r in rr], reference={}, states={})
        for name in ('uniform', 'spt'):
            ok = [abs(np.log1p(r['reference'][name]['z']) - np.log1p(r['z_published'])) < TOL_LN for r in rr]
            conf = [r['reference'][name]['confidence'] >= CONF for r in rr]
            g['reference'][name] = dict(confident_correct=int(sum(a and b for a, b in zip(ok, conf))),
                                        confident_wrong=int(sum((not a) and b for a, b in zip(ok, conf))),
                                        best_is_published=int(sum(ok)))
        for state in STATES:
            qq = [r['states'][state] for r in rr]
            s = dict(engine=dict(
                confident_correct=sum(q['engine']['outcome'] == 'correct' and q['engine']['confidence'] >= CONF for q in qq),
                confident_wrong=sum(q['engine']['outcome'] == 'wrong' and q['engine']['confidence'] >= CONF for q in qq),
                best_is_published=sum(q['engine']['outcome'] == 'correct' for q in qq)))
            for key in [k for k in qq[0] if k.startswith(('guided_', 'window_', 'pick_')) or k == 'xcorr']:
                s[key] = {o: sum(q[key]['outcome'] == o for q in qq) for o in ('correct', 'wrong', 'none')}
            g['states'][state] = s
        out[group] = g
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--plot-only', action='store_true', help='re-plot from results/redshift_photoz.json')
    a = ap.parse_args()
    path = RESULTS / 'redshift_photoz.json'
    if a.plot_only:
        plot(json.loads(path.read_text()))
        plot_noise()
        print('re-plotted the saved photometric-redshift test')
        return
    rows = analyse()
    out = dict(confidence_threshold=CONF, window_half_width=WINDOW, recalibration=RECAL, states=STATES,
               note='prior = Gaussian in z at the photometric redshift; pull = (published z - centre) / width; '
                    'p_published = posterior mass within the confidence window of the published redshift; '
                    'guided, window and pick are given with a missing [CI] line counted against a redshift '
                    '(ci_counted) and with CO lines only (co_only)',
               summary=summarise(rows), sources=rows)
    path.write_text(json.dumps(out, indent=1) + '\n')
    for group, title in (('one', 'One line in the scan'), ('several', 'Several lines')):
        g = out['summary'][group]
        print(f"\n{title}: {g['n']} sources. Engine with the survey prior: confident and correct "
              f"{g['reference']['spt']['confident_correct']}, best guess right {g['reference']['spt']['best_is_published']}")
        for state, label in STATES.items():
            s = g['states'][state]
            line = (f"  {label}: engine confident and correct {s['engine']['confident_correct']}, confident and "
                    f"wrong {s['engine']['confident_wrong']}, best guess right {s['engine']['best_is_published']}")
            for key in ('guided_ci_counted', 'guided_co_only', 'window_ci_counted', 'window_co_only', 'xcorr'):
                if key in s:
                    line += f" | {key} {s[key]['correct']}/{s[key]['wrong']}"
            print(line)
    print(f'\n(rules and cross-correlation: correct/wrong)  wrote {path.name}')
    plot(out)
    plot_noise()


def _style(ax):
    ax.set_facecolor(SURF); ax.set_axisbelow(True)
    ax.tick_params(colors=INK2, labelsize=8.5)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)


def plot(o):
    """The one-line sources. Left: the redshifts the scan allows and where the photometric redshift falls.
    Right: the engine's probability for the published redshift under each prior."""
    rows = sorted((r for r in o['sources'] if r['lines'] == 'one'), key=lambda r: r['z_published'])
    y = np.arange(len(rows))[::-1].astype(float)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(12.4, 6.0), facecolor=SURF, sharey=True,
                                 gridspec_kw=dict(width_ratios=[1.25, 1.0], wspace=0.06))
    blue, orange = COL[0], COL[1]
    for yy, r in zip(y, rows):
        for state, colr, dy in (('photoz', blue, 0.17), ('photoz_2024', orange, -0.17)):
            c, w = r['states'][state]['centre'], r['states'][state]['width']
            ax.plot([c - w, c + w], [yy + dy, yy + dy], color=colr, lw=5, alpha=0.30, solid_capstyle='butt', zorder=1)
            ax.plot([c], [yy + dy], '|', color=colr, ms=9, mew=2, zorder=2)
        others = [zc for zc, _ in r['candidates_no_prior'] if abs(np.log1p(zc) - np.log1p(r['z_published'])) >= TOL_LN]
        ax.plot(others, [yy] * len(others), 'o', mfc=SURF, mec=INK2, ms=7, mew=1.3, zorder=3)
        ax.plot([r['z_published']], [yy], 'o', color=INK, ms=7.5, mec=SURF, mew=1.2, zorder=4)
    ax.plot([], [], 'o', color=INK, ms=7.5, mec=SURF, mew=1.2, label='published redshift')
    ax.plot([], [], 'o', mfc=SURF, mec=INK2, ms=7, mew=1.3, label='other redshifts the scan allows (engine, no prior, p ≥ 0.02)')
    ax.plot([], [], color=blue, lw=5, alpha=0.30, label='photometric redshift ± error, 2020 paper')
    ax.plot([], [], color=orange, lw=5, alpha=0.30, label='the same, 2024 revision')
    ax.set_xlim(0.9, 6.6)
    ax.set_xlabel('redshift', fontsize=9, color=INK2)
    ax.set_yticks(y)
    ax.set_yticklabels([r['source'] for r in rows], fontsize=8.5, color=INK)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.grid(axis='x', color=GRID_C, lw=0.8)
    ax.set_title('A. What one line allows, and where the photometric redshift falls', fontsize=10, color=INK, loc='left')
    ax.legend(fontsize=8.2, frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.10), ncol=2, labelcolor=INK,
              handletextpad=0.5, columnspacing=1.6)
    _style(ax)

    logit = lambda p: np.log10(np.clip(p, 2e-3, 0.9995) / (1 - np.clip(p, 2e-3, 0.9995)))
    series = [('survey redshift distribution', lambda r: r['reference']['spt']['p_published'], 'o', INK2, False),
              ('photometric redshift, 2020 paper', lambda r: r['states']['photoz']['engine']['p_published'], 'o', blue, True),
              ('2024 revision', lambda r: r['states']['photoz_2024']['engine']['p_published'], 'D', orange, True),
              ('2020, recalibrated', lambda r: r['states']['photoz_recal']['engine']['p_published'], 's', COL[2], True)]
    for yy in y:
        bx.axhline(yy, color=GRID_C, lw=0.8, zorder=0)
    for (label, f, mk, colr, filled), dy in zip(series, (0.0, 0.2, 0.0, -0.2)):
        bx.plot([logit(f(r)) for r in rows], y + dy, mk, color=colr, mfc=colr if filled else SURF, mec=SURF if filled else colr,
                ms=7.5 if mk == 'o' else 6.5, mew=1.2 if filled else 1.5, label=label, zorder=3 if filled else 2, clip_on=False)
    ticks = [0.01, 0.1, 0.5, 0.9, 0.99, 0.999]
    bx.set_xticks([logit(t) for t in ticks]); bx.set_xticklabels([f'{t:g}' for t in ticks])
    bx.set_xlim(logit(2e-3) - 0.1, logit(0.9995) + 0.1)
    bx.axvline(logit(CONF), color=INK2, lw=0.9, ls=':')
    bx.text(logit(CONF) - 0.06, len(rows) - 0.36, 'confident\nfrom 0.95', fontsize=8, color=INK2, ha='right', va='top')
    bx.set_xlabel("engine's probability for the published redshift (logit scale)", fontsize=9, color=INK2)
    bx.set_title('B. The engine, by the prior it is given', fontsize=10, color=INK, loc='left')
    bx.legend(fontsize=8.2, frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.10), ncol=2, labelcolor=INK,
              handletextpad=0.3, columnspacing=1.6, title='redshift prior', title_fontsize=8.2, alignment='left')
    _style(bx)
    one = o['summary']['one']
    counts = ', '.join(f"{one['states'][s]['engine']['confident_correct']} ({lab})" for s, lab in
                       (('photoz', '2020 paper'), ('photoz_2024', '2024 revision'), ('photoz_recal', 'recalibrated')))
    fig.suptitle(f"The {one['n']} sources with one line in the 3 mm scan, each with its published photometric redshift "
                 f"as the prior.\nConfident and correct: {counts}; with the survey distribution "
                 f"{one['reference']['spt']['confident_correct']}. Confident and wrong: none.",
                 fontsize=10.5, color=INK, x=0.012, y=0.99, ha='left')
    fig.subplots_adjust(left=0.085, right=0.985, top=0.855, bottom=0.2)
    fig.savefig(FIG / 'redshift_photoz.png', dpi=130, facecolor=SURF)
    plt.close(fig)


def plot_noise():
    """The noise-added test with each source's 2020 photometric redshift as the prior: answers given, right and
    wrong. Drawn from results/redshift_degradation.json; skipped if that file has no photometric columns."""
    path = RESULTS / 'redshift_degradation.json'
    o = json.loads(path.read_text()) if path.exists() else {}
    rows = o.get('summary', [])
    if not rows or 'guided_co_only' not in rows[0]:
        return
    x = [q['target'] for q in rows]
    blue, orange = COL[0], COL[1]
    series = [('engine, photometric prior (confident)', blue, 'o', '-', True, lambda q, k: q['photoz']['engine'][k[0]]),
              ('engine, no prior (confident)', INK2, 'o', '--', False, lambda q, k: q['uniform']['engine'][k[0]]),
              ('guided rule, photometric prior, [CI] counted', orange, 'D', '-', True,
               lambda q, k: q['guided_photoz']['photoz'][k[1]]),
              ('guided rule, photometric prior, CO lines only', orange, 'D', '--', False,
               lambda q, k: q['guided_co_only']['photoz'][k[1]])]
    panels = [('A. The published redshift, returned as an answer', ('confident_correct', 'correct'), (0, 1.03)),
              ('B. A different redshift, returned as an answer', ('confidently_wrong', 'wrong'), (0, 0.155))]
    fig, axs = plt.subplots(1, 2, figsize=(12.4, 5.3), facecolor=SURF)
    for ax, (title, keys, ylim) in zip(axs, panels):
        for name, colr, mk, ls, filled, f in series:
            ax.plot(x, [f(q, keys) for q in rows], color=colr, marker=mk, ms=7 if mk == 'o' else 6, lw=2, ls=ls,
                    mfc=colr if filled else SURF, mec=SURF if filled else colr, mew=1.4, label=name, clip_on=False)
        ax.set_title(title, fontsize=10, color=INK, loc='left')
        ax.set_xticks(x); ax.set_ylim(*ylim)
        ax.set_xlabel('matched-filter S/N of the second-brightest line after adding noise', fontsize=9, color=INK2)
        ax.grid(axis='y', color=GRID_C, lw=0.8)
        _style(ax)
    axs[0].set_ylabel('fraction of noise draws (mean of the sources)', fontsize=9, color=INK2)
    axs[0].legend(fontsize=8.5, frameon=False, loc='lower right', labelcolor=INK, handlelength=4.2)
    lv = [q for s_ in o['sources'].values() for q in s_['levels']]
    total = sum(q['draws'] for q in lv)
    count = lambda f: int(round(sum(f(q) * q['draws'] for q in lv)))
    wrong = {p: (count(lambda q: q[p]['engine']['confidently_wrong']), count(lambda q: q['guided_photoz'][p]['wrong']),
                 count(lambda q: q['guided_co_only'][p]['wrong'])) for p in STATES}
    n_src, n = rows[0]['n_sources'], rows[0]['draws_per_source']
    fig.suptitle(f'Real multi-line spectra with noise added, each source with its 2020 photometric redshift as the '
                 f'prior: {n_src} sources, {n} noise draws per source and level', fontsize=10.5, color=INK, x=0.012,
                 ha='left')
    fig.text(0.012, 0.012, f'A different redshift returned, of {total:,} draws in all (engine; guided rule; guided rule '
             'with CO lines only): ' + '; '.join(f"{lab} {w[0]}, {w[1]}, {w[2]}" for lab, w in
                                                 (('2020 values', wrong['photoz']), ('2024 revision', wrong['photoz_2024']),
                                                  ('recalibrated', wrong['photoz_recal']))) + '.',
             fontsize=8, color=INK2, ha='left', va='bottom')
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    fig.savefig(FIG / 'redshift_photoz_noise.png', dpi=130, facecolor=SURF)
    plt.close(fig)


if __name__ == '__main__':
    main()
