"""Do the real multi-line spectra still give the published redshift once the second line sinks below 5 sigma?

As observed, the seven sources of 2016.1.00672.S whose scan shows two or three lines are settled both by the
engine and by the detection rule: their second-brightest lines are at 5.6 sigma or more, above the rule's
threshold. So the real data, as they are, do not test the engine where it should gain.

Here Gaussian noise is added to each real spectrum (realspec.degrade) until its second-brightest line would have
a chosen matched-filter S/N (5.5 down to 2), and five methods are run on many noise draws:
  - the engine, with no redshift prior and with the survey prior;
  - the equal-weight sum;
  - the detection rule (every line at 5 sigma);
  - a guided rule that, after the 5 sigma lines, accepts a weaker line (3 or 4 sigma) at a position one of the
    surviving solutions predicts (ladder.guided_rule), with and without the survey prior. It stands for what an
    observer does by eye; the 5 sigma rule alone would be too easy to beat;
  - a template cross-correlation after Yun et al. 2015 (ladder.xcorr_redshift), with and without the prior.
Each source's published photometric redshift is run as a third kind of prior, in the three versions of
photoz.py, for the engine, the guided rule and the cross-correlation. The guided rule is also run in a form
that does not count a missing [CI] line against a redshift (guided_co_only), under every prior.
Two kinds of outcome are kept apart. The engine (when confident) and the rules give an answer or none, so they
are compared on answers given. The cross-correlation always has a best redshift, so it is compared with the
engine's best guess. The comparison without a prior is the like-for-like one: the rule and the
cross-correlation use only the spectrum. The line shapes, line ratios, continuum residuals and window layout
are the real ones; only the extra noise is simulated.

Two limits on what this shows:
  - each draw keeps the noise of the real spectrum and adds to it, so draws of one source are not independent,
    least of all at the first levels (a share 1/factor^2 of the noise variance is the original);
  - these seven spectra are the ones the engine was developed on.

Writes results/redshift_degradation.json, figures/redshift_degradation.png (means) and
figures/redshift_degradation_sources.png (per source).
Usage: python -m ladder_redshifts.degrade [N_draws] [--plot-only]     (default 100; ~10 min on 2 cores)
"""
from .paths import FIG, RESULTS, DATA  # noqa: E402
import argparse
import json
import multiprocessing as mp
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .ladder import (line_list, matched_filter, detection_index, detect, guided_rule, xcorr_redshift,  # noqa: E402
                     prior_spt, BANK, C_KMS, TOL_LN)
from .batch_real import read_spectrum, to_grid_auto, read_truth  # noqa: E402
from .realspec import run_engine, engine_likelihood, engine_posterior, degrade  # noqa: E402
from .photoz import photoz_priors, STATES as PHOTOZ_STATES  # noqa: E402

INK, INK2, GRID_C, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
COL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300']
TARGETS = [5.5, 5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0]          # matched-filter S/N of the second-brightest line
PRIORS = {'spt': prior_spt, 'uniform': None}
CONF = 0.95
GUIDED = (3.0, 4.0)                                         # confirmation thresholds of the guided rule
SEED = 20261113


def multi_line_sources():
    """(source, published z) for the sources whose published scan shows two or three lines (data/truth.csv)."""
    out = []
    for s, t in sorted(read_truth().items()):
        if ('(2 lines)' in t['notes'] or '(3 lines)' in t['notes']) and (DATA / 'spectra' / f'{s}.csv').exists():
            out.append((s, float(t['z_published'])))
    return out


def published_lines(detections, z):
    """The detections that are lines of redshift z, as (name, frequency, S/N), brightest first."""
    names, nus = line_list()
    out = []
    for f, s in detections:
        off = np.log(nus / (1 + z) / f)
        j = int(np.argmin(np.abs(off)))
        if abs(off[j]) * C_KMS < 300.0:
            out.append((names[j], float(f), float(s)))
    return sorted(out, key=lambda q: -q[2])


def correct(z, z_true):
    return z is not None and abs(np.log1p(z) - np.log1p(z_true)) < TOL_LN


def run_level(task):
    """All draws of one source at one noise level. Returns counts (see main)."""
    src, z_true, f2, target, factor, n, seed = task
    spec = read_spectrum(DATA / 'spectra' / f'{src}.csv')
    rng = np.random.default_rng(seed)
    priors = dict(PRIORS, **photoz_priors(src))             # the source's own photometric redshift, three versions
    photoz = [p for p in priors if p in PHOTOZ_STATES]
    c = dict(rule_ok=0, rule_wrong=0, rule_none=0, second_detected=0, y2=[])
    for t in GUIDED:
        c[f'guided{t}'] = dict(ok=0, wrong=0, none=0)
    for p in ['spt'] + photoz:
        c[f'guided_{p}'] = dict(ok=0, wrong=0, none=0)
    for p in priors:
        c[f'guidedco_{p}'] = dict(ok=0, wrong=0, none=0)
        c[f'{p}_xcorr'] = 0
        for m in ('engine', 'flat'):
            c[f'{p}_{m}'] = dict(cc=0, cw=0, top=0, conf=[])
    for _ in range(n):
        g = to_grid_auto(degrade(spec, factor, rng))
        grid, d, _, cov = g
        like = engine_likelihood(*g)                        # once per draw: it does not depend on the prior
        yd, dets = like['yd'], like['dets']
        c['y2'].append(float(yd[int(round(grid.pix(f2)))]))
        for t in GUIDED:
            zg = guided_rule(dets, grid, cov, yd, confirm=t)[0]
            q = c[f'guided{t}']
            q['ok'] += correct(zg, z_true); q['wrong'] += zg is not None and not correct(zg, z_true)
            q['none'] += zg is None
        for p in ['spt'] + photoz:
            zg = guided_rule(dets, grid, cov, yd, confirm=GUIDED[0], prior=priors[p])[0]
            q = c[f'guided_{p}']
            q['ok'] += correct(zg, z_true); q['wrong'] += zg is not None and not correct(zg, z_true)
            q['none'] += zg is None
        for p, prior in priors.items():
            # the same rule when a missing [CI] line is not counted against a redshift (ladder._rule_solutions)
            zg = guided_rule(dets, grid, cov, yd, confirm=GUIDED[0], prior=prior, ci_prunes=False)[0]
            q = c[f'guidedco_{p}']
            q['ok'] += correct(zg, z_true); q['wrong'] += zg is not None and not correct(zg, z_true)
            q['none'] += zg is None
            c[f'{p}_xcorr'] += correct(xcorr_redshift(d, grid, cov, prior=prior)[0], z_true)
            r = engine_posterior(like, prior)
            for m, key in (('engine', 'sled'), ('flat', 'flat')):
                z, conf = r[key]
                ok = correct(z, z_true)
                q = c[f'{p}_{m}']
                q['cc'] += ok and conf >= CONF; q['cw'] += (not ok) and conf >= CONF; q['top'] += ok
                q['conf'].append(conf)
        zr = r['sudoku'][0]                                   # the rule does not use the prior
        c['rule_ok'] += correct(zr, z_true); c['rule_wrong'] += zr is not None and not correct(zr, z_true)
        c['rule_none'] += zr is None
        c['second_detected'] += len(published_lines(r['detections'], z_true)) >= 2
    out = dict(target=target, factor=round(factor, 3), draws=n,
               second_line_snr_measured=round(float(np.median(c['y2'])), 2),
               second_line_detected=c['second_detected'] / n,
               rule=dict(correct=c['rule_ok'] / n, wrong=c['rule_wrong'] / n, no_answer=c['rule_none'] / n),
               guided={f'{t:g}': dict(correct=c[f'guided{t}']['ok'] / n, wrong=c[f'guided{t}']['wrong'] / n,
                                      no_answer=c[f'guided{t}']['none'] / n) for t in GUIDED},
               guided_spt=dict(correct=c['guided_spt']['ok'] / n, wrong=c['guided_spt']['wrong'] / n,
                               no_answer=c['guided_spt']['none'] / n),
               guided_photoz={p: dict(correct=c[f'guided_{p}']['ok'] / n, wrong=c[f'guided_{p}']['wrong'] / n,
                                      no_answer=c[f'guided_{p}']['none'] / n) for p in photoz},
               guided_co_only={p: dict(correct=c[f'guidedco_{p}']['ok'] / n, wrong=c[f'guidedco_{p}']['wrong'] / n,
                                       no_answer=c[f'guidedco_{p}']['none'] / n) for p in priors},
               xcorr={p: dict(best_is_published=c[f'{p}_xcorr'] / n) for p in priors})
    for p in priors:
        out[p] = {m: dict(confident_correct=c[f'{p}_{m}']['cc'] / n, confidently_wrong=c[f'{p}_{m}']['cw'] / n,
                          best_is_published=c[f'{p}_{m}']['top'] / n,
                          median_confidence=round(float(np.median(c[f'{p}_{m}']['conf'])), 3))
                  for m in ('engine', 'flat')}
    return src, out


def as_observed():
    """Every method on every real spectrum that has a published redshift, with no noise added."""
    kdet = detection_index(BANK)
    rows = []
    for src, t in sorted(read_truth().items()):
        path = DATA / 'spectra' / f'{src}.csv'
        if not path.exists():
            continue
        z = float(t['z_published'])
        g = to_grid_auto(read_spectrum(path))
        grid, d, _, cov = g
        yd = matched_filter(d, BANK, cov)[kdet]
        dets = detect(yd, cov)
        outcome = lambda zz: 'none' if zz is None else 'correct' if correct(zz, z) else 'wrong'
        row = dict(source=src, z_published=z, detections=len(dets))
        for t_ in GUIDED:
            zg = guided_rule(dets, grid, cov, yd, confirm=t_)[0]
            row[f'guided_{t_:g}'] = dict(z=None if zg is None else round(zg, 4), outcome=outcome(zg))
        zg = guided_rule(dets, grid, cov, yd, confirm=GUIDED[0], prior=prior_spt)[0]
        row['guided_3_spt'] = dict(z=None if zg is None else round(zg, 4), outcome=outcome(zg))
        row['engine'], row['xcorr'] = {}, {}
        for p, prior in PRIORS.items():
            ze, ce = run_engine(*g, prior=prior)['sled']
            row['engine'][p] = dict(z=round(ze, 4), confidence=round(ce, 3), best_is_published=bool(correct(ze, z)))
            zx, sx = xcorr_redshift(d, grid, cov, prior=prior)
            row['xcorr'][p] = dict(z=round(zx, 4), snr=round(sx, 1), best_is_published=bool(correct(zx, z)))
        rows.append(row)
    return rows


def summarise(sources):
    """Mean over the sources at each level (each source counts equally)."""
    rows = []
    for k, target in enumerate(TARGETS):
        lv = [s['levels'][k] for s in sources.values()]
        row = dict(target=target, n_sources=len(lv), draws_per_source=lv[0]['draws'],
                   second_line_detected=round(float(np.mean([q['second_line_detected'] for q in lv])), 4),
                   rule={kk: round(float(np.mean([q['rule'][kk] for q in lv])), 4) for kk in lv[0]['rule']},
                   guided={t: {kk: round(float(np.mean([q['guided'][t][kk] for q in lv])), 4)
                               for kk in lv[0]['guided'][t]} for t in lv[0]['guided']},
                   guided_spt={kk: round(float(np.mean([q['guided_spt'][kk] for q in lv])), 4)
                               for kk in lv[0]['guided_spt']},
                   guided_photoz={p: {kk: round(float(np.mean([q['guided_photoz'][p][kk] for q in lv])), 4)
                                      for kk in lv[0]['guided_photoz'][p]} for p in lv[0]['guided_photoz']},
                   guided_co_only={p: {kk: round(float(np.mean([q['guided_co_only'][p][kk] for q in lv])), 4)
                                       for kk in lv[0]['guided_co_only'][p]} for p in lv[0]['guided_co_only']},
                   xcorr={p: dict(best_is_published=round(float(np.mean([q['xcorr'][p]['best_is_published']
                                                                        for q in lv])), 4)) for p in lv[0]['xcorr']})
        for p in lv[0]['xcorr']:
            row[p] = {m: {kk: round(float(np.mean([q[p][m][kk] for q in lv])), 4)
                          for kk in ('confident_correct', 'confidently_wrong', 'best_is_published')}
                      for m in ('engine', 'flat')}
        rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('n', nargs='?', type=int, default=100)
    ap.add_argument('--plot-only', action='store_true', help='re-plot from results/redshift_degradation.json')
    a = ap.parse_args()
    path = RESULTS / 'redshift_degradation.json'
    if a.plot_only:
        plot(json.loads(path.read_text()))
        print('re-plotted the saved degradation test')
        return
    t0 = time.time()
    sources, tasks = {}, []
    for i, (src, z) in enumerate(multi_line_sources()):
        r = run_engine(*to_grid_auto(read_spectrum(DATA / 'spectra' / f'{src}.csv')), prior=prior_spt)
        lines = published_lines(r['detections'], z)
        name2, f2, snr2 = lines[1]
        sources[src] = dict(z_published=z, lines_observed=[(n, round(f, 3), round(s, 1)) for n, f, s in lines],
                            second_line=dict(line=name2, freq_ghz=round(f2, 3), snr=round(snr2, 2)),
                            as_observed=dict(engine=dict(z=round(r['sled'][0], 4), confidence=round(r['sled'][1], 3)),
                                             rule_z=None if r['sudoku'][0] is None else round(r['sudoku'][0], 4)),
                            levels=[])
        print(f'{src}: lines as observed ' + ', '.join(f'{n} {s:.1f}' for n, _, s in lines))
        for k, target in enumerate(TARGETS):
            tasks.append((src, z, f2, target, snr2 / target, a.n, [SEED, i, k]))
    assert all(t[4] >= 1.0 for t in tasks), 'a target S/N above the observed one cannot be reached by adding noise'
    observed = as_observed()
    with mp.get_context('fork').Pool(2) as pool:
        for src, out in pool.imap(run_level, tasks):
            sources[src]['levels'].append(out)
    out = dict(n_draws=a.n, seed=SEED, targets=TARGETS, confidence_threshold=CONF,
               note='noise added to the real spectra; second line = second-brightest detected line of the '
                    'published redshift; fractions are over draws, the summary is the mean over sources',
               sources=sources, summary=summarise(sources), as_observed=observed)
    path.write_text(json.dumps(out, indent=1) + '\n')
    print(f'\n{a.n} draws x {len(sources)} sources. Columns: second-line S/N | no prior: engine confident+correct, '
          'guided rule correct, engine best guess, cross-correlation best guess | survey prior: the same four '
          '| 5 sigma rule')
    for q in out['summary']:
        u, sp = q['uniform'], q['spt']
        print(f"  {q['target']:4.1f} | {u['engine']['confident_correct']:.3f} {q['guided']['3']['correct']:.3f} "
              f"{u['engine']['best_is_published']:.3f} {q['xcorr']['uniform']['best_is_published']:.3f} | "
              f"{sp['engine']['confident_correct']:.3f} {q['guided_spt']['correct']:.3f} "
              f"{sp['engine']['best_is_published']:.3f} {q['xcorr']['spt']['best_is_published']:.3f} | "
              f"{q['rule']['correct']:.3f}")
    print('With each source\'s photometric redshift as the prior. Columns: second-line S/N | then, for the 2020 '
          'values, the 2024 revision and the recalibrated 2020 values: engine confident+correct, engine confident '
          'and wrong, guided rule correct, guided rule wrong, the same two with CO lines only, engine best guess, '
          'cross-correlation best guess')
    for q in out['summary']:
        print(f"  {q['target']:4.1f} | " + ' | '.join(
            f"{q[p]['engine']['confident_correct']:.3f} {q[p]['engine']['confidently_wrong']:.3f} "
            f"{q['guided_photoz'][p]['correct']:.3f} {q['guided_photoz'][p]['wrong']:.3f} "
            f"{q['guided_co_only'][p]['correct']:.3f} {q['guided_co_only'][p]['wrong']:.3f} "
            f"{q[p]['engine']['best_is_published']:.3f} {q['xcorr'][p]['best_is_published']:.3f}"
            for p in q['guided_photoz']))
    print(f'[{time.time() - t0:.0f}s] wrote {path.name}')
    plot(out)


def wrong_counts(o):
    """Draws, over all sources and levels, in which each method returned a different redshift with confidence."""
    lv = [q for s in o['sources'].values() for q in s['levels']]
    n = sum(q['draws'] for q in lv)
    count = lambda f: int(round(sum(f(q) * q['draws'] for q in lv)))
    return n, {'engine, no prior': count(lambda q: q['uniform']['engine']['confidently_wrong']),
               'engine, survey prior': count(lambda q: q['spt']['engine']['confidently_wrong']),
               'guided rule': count(lambda q: q['guided']['3']['wrong']),
               'guided rule, survey prior': count(lambda q: q['guided_spt']['wrong']),
               'detection rule': count(lambda q: q['rule']['wrong'])}


def wrong_counts_photoz(o):
    """The same count with each source's photometric redshift as the prior:
    {version: (engine, guided rule, guided rule with CO lines only)}."""
    lv = [q for s in o['sources'].values() for q in s['levels']]
    count = lambda f: int(round(sum(f(q) * q['draws'] for q in lv)))
    return {p: (count(lambda q: q[p]['engine']['confidently_wrong']), count(lambda q: q['guided_photoz'][p]['wrong']),
                count(lambda q: q['guided_co_only'][p]['wrong'])) for p in lv[0]['guided_photoz']}


def _style(ax):
    ax.set_facecolor(SURF); ax.set_axisbelow(True)
    ax.tick_params(colors=INK2, labelsize=8.5)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)


def plot(o):
    plot_curves(o)
    plot_sources(o)


def plot_curves(o):
    """Mean over the sources against the S/N left to the second line. Solid: no prior. Dashed: survey prior."""
    rows = o['summary']
    x = [q['target'] for q in rows]
    eng, rule3, rule5, xc = COL[0], COL[1], COL[2], COL[3]
    panels = [
        ('A. Answers given: the published redshift, returned', [
            ('engine, no prior (confident)', eng, 'o', '-', lambda q: q['uniform']['engine']['confident_correct']),
            ('engine, survey prior (confident)', eng, 'o', '--', lambda q: q['spt']['engine']['confident_correct']),
            ('guided rule: 5σ lines, then one at 3σ', rule3, 'D', '-', lambda q: q['guided']['3']['correct']),
            ('detection rule: every line at 5σ', rule5, 's', '-', lambda q: q['rule']['correct'])]),
        ('B. Best guesses: the published redshift ranked first', [
            ('engine, no prior', eng, 'o', '-', lambda q: q['uniform']['engine']['best_is_published']),
            ('engine, survey prior', eng, 'o', '--', lambda q: q['spt']['engine']['best_is_published']),
            ('cross-correlation, no prior', xc, '^', '-', lambda q: q['xcorr']['uniform']['best_is_published']),
            ('cross-correlation × survey prior', xc, '^', '--', lambda q: q['xcorr']['spt']['best_is_published'])])]
    fig, axs = plt.subplots(1, 2, figsize=(12.4, 5.3), facecolor=SURF, sharey=True)
    for ax, (title, series), where in zip(axs, panels, ('upper left', 'lower right')):
        for name, colr, mk, ls, f in series:
            dashed = ls == '--'
            ax.plot(x, [f(q) for q in rows], color=colr, marker=mk, ms=7, lw=2, ls=ls, mec=colr if dashed else SURF,
                    mfc=SURF if dashed else colr, mew=1.4, label=name, clip_on=False)
        ax.set_title(title, fontsize=10, color=INK, loc='left')
        ax.set_xticks(x); ax.set_ylim(0, 1.03)
        ax.set_xlabel('matched-filter S/N of the second-brightest line after adding noise', fontsize=9, color=INK2)
        ax.grid(axis='y', color=GRID_C, lw=0.8)
        ax.legend(fontsize=8.5, frameon=False, loc=where, labelcolor=INK, handlelength=4.2)
        _style(ax)
    axs[0].set_ylabel('fraction of noise draws (mean of the sources)', fontsize=9, color=INK2)
    axs[0].axvline(5.0, color=INK2, lw=0.8, ls=':')
    axs[0].text(4.97, 0.03, "the detection\nrule's threshold", fontsize=8, color=INK2, va='bottom', ha='right')
    n_src, n = rows[0]['n_sources'], rows[0]['draws_per_source']
    fig.suptitle(f'Real multi-line spectra with noise added: {n_src} sources of 2016.1.00672.S, {n} noise draws '
                 'per source and level', fontsize=10.5, color=INK, x=0.012, ha='left')
    total, wrong = wrong_counts(o)
    fig.text(0.012, 0.012, f'A different redshift returned as an answer, of {total:,} draws in all: '
             + '; '.join(f'{k} {v}' for k, v in wrong.items()) + '.', fontsize=8, color=INK2, ha='left', va='bottom')
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    fig.savefig(FIG / 'redshift_degradation.png', dpi=130, facecolor=SURF)
    plt.close(fig)


def plot_sources(o, levels=(4.5, 4.0, 3.5, 3.0)):
    """One row per source, one panel per S/N of the second line: where the engine is ahead and where it is not."""
    names = list(o['sources'])
    fig, axs = plt.subplots(1, len(levels), figsize=(12.4, 4.6), facecolor=SURF, sharey=True)
    y = np.arange(len(names))[::-1]
    for ax, t in zip(axs, levels):
        k = o['targets'].index(t)
        lv = [o['sources'][s]['levels'][k] for s in names]
        e = np.array([q['uniform']['engine']['confident_correct'] for q in lv])
        ep = np.array([q['spt']['engine']['confident_correct'] for q in lv])
        g = np.array([q['guided']['3']['correct'] for q in lv])
        r = np.array([q['rule']['correct'] for q in lv])
        ax.hlines(y, np.minimum(e, g), np.maximum(e, g), color='#c3c2b7', lw=2, zorder=1)
        ax.plot(r, y, 's', color=COL[2], ms=7, mec=SURF, mew=1.2, zorder=2, label='detection rule: every line at 5σ')
        ax.plot(g, y, 'D', color=COL[1], ms=7, mec=SURF, mew=1.2, zorder=3,
                label='guided rule: 5σ lines, then one at 3σ')
        ax.plot(ep, y, 'o', mfc=SURF, mec=COL[0], ms=7.5, mew=1.5, zorder=4, label='engine, survey prior')
        ax.plot(e, y, 'o', color=COL[0], ms=7.5, mec=SURF, mew=1.2, zorder=5, label='engine, no prior')
        ahead, behind = int((e > g).sum()), int((e < g).sum())
        ax.set_title(f'second line at {t:g}σ', fontsize=10, color=INK, loc='left')
        ax.set_xlabel(f'engine (no prior) ahead of the guided rule\nfor {ahead} of {len(names)}, behind for {behind}',
                      fontsize=8.5, color=INK2)
        ax.set_xlim(-0.04, 1.04); ax.set_xticks([0, 0.5, 1]); ax.set_xticklabels(['0', '0.5', '1'])
        ax.set_ylim(-0.6, len(names) - 0.4)
        ax.grid(axis='x', color=GRID_C, lw=0.8)
        _style(ax)
    axs[0].set_yticks(y)
    axs[0].set_yticklabels([f"{s}\n{o['sources'][s]['second_line']['line'].replace('-', '–')}" for s in names],
                           fontsize=8, color=INK)
    h, lab = axs[0].get_legend_handles_labels()
    fig.legend(h[::-1], lab[::-1], fontsize=8.5, frameon=False, loc='upper left', ncol=4, labelcolor=INK,
               bbox_to_anchor=(0.006, 0.955), handletextpad=0.3, columnspacing=1.4)
    n = o['summary'][0]['draws_per_source']
    fig.suptitle(f'Per source: fraction of {n} noise draws in which the published redshift is returned '
                 '(engine: with confidence). Each source is labelled with its second line.',
                 fontsize=10.5, color=INK, x=0.012, y=0.995, ha='left')
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(FIG / 'redshift_degradation_sources.png', dpi=130, facecolor=SURF)
    plt.close(fig)


if __name__ == '__main__':
    main()
