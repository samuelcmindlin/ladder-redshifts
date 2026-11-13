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

A comparison with one hand-built rule shows little: its thresholds decide the outcome. So a family of settings
of the guided rule (FAMILY: 72 with no prior, 216 with a photometric prior) is run on the same draws, and the
engine's answers are counted at several confidence thresholds (ENGINE_THRESHOLDS). The rule the engine is
compared with is chosen by a criterion fixed beforehand (see WRONG_CAP), and the two are also compared at equal
numbers of wrong answers.

Two kinds of outcome are kept apart. The engine (when confident) and the rules give an answer or none, so they
are compared on answers given. The cross-correlation always has a best redshift, so it is compared with the
engine's best guess. The comparison without a prior is the like-for-like one: the rule and the
cross-correlation use only the spectrum. The line shapes, line ratios, continuum residuals and window layout
are the real ones; only the extra noise is simulated.

Two limits on what this shows:
  - each draw keeps the noise of the real spectrum and adds to it, so draws of one source are not independent,
    least of all at the first levels (a share 1/factor^2 of the noise variance is the original);
  - these seven spectra are the ones the engine was developed on.

Writes results/redshift_degradation.json, results/redshift_rule_family.csv (one row per rule setting and prior),
figures/redshift_degradation.png (means), figures/redshift_degradation_sources.png (per source) and
figures/redshift_operating_points.png (answers against wrong answers).
Usage: python -m ladder_redshifts.degrade [N_draws] [--plot-only]     (default 100; ~16 min on 2 cores)
"""
from .paths import FIG, RESULTS, DATA  # noqa: E402
import argparse
import csv
import itertools
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
# A family of rule settings, so that the comparison does not rest on one hand-built rule. Each is the guided
# rule with: the S/N at which a line counts as detected; the S/N at which a weaker line confirms a solution;
# whether a missing [CI] line counts against a redshift; and the expected S/N above which a missing line counts
# against one. Every setting is run with no prior and with the source's 2020 photometric redshift.
FAMILY = dict(detect=(5.0, 4.5, 4.0), confirm=(4.0, 3.5, 3.0, 2.5), ci_counted=(True, False),
              prune_expect=(5.0, 8.0, 12.0))
FAMILY_SHARE = (0.95, 0.90, 0.80)       # with a prior: the share of prior weight at which one solution is taken
FAMILY_PRIORS = ('uniform', 'photoz')
# The engine's answers are counted at these confidence thresholds as well as at CONF, so that the engine and the
# rules can be compared at the same number of wrong answers.
ENGINE_THRESHOLDS = (0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)
# How the rule to compare with is chosen, fixed before the family was run. Over all noise draws, the setting with
# the most correct answers among those that return a wrong redshift (a) no more often than the engine does with
# the same prior, or (b) in at most 1% of draws. Ties go to the first setting in the order of rule_family().
# The envelope takes, at each noise level separately, the best setting with at most 1% wrong at that level; it
# is chosen with knowledge of the outcome, so it is an upper limit on any one rule of the family.
WRONG_CAP = 0.01


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


def rule_family(with_prior=False):
    """The settings of FAMILY, in a fixed order. With a prior each is run at every value of FAMILY_SHARE."""
    base = [dict(zip(FAMILY, v)) for v in itertools.product(*FAMILY.values())]
    if not with_prior:
        return [dict(v, share=None) for v in base]
    return [dict(v, share=sh) for sh in FAMILY_SHARE for v in base]


def family_answers(yd, grid, cov, prior=None):
    """The answer (z or None) of every rule setting for one matched-filter trace."""
    dets = {t: detect(yd, cov, t) for t in FAMILY['detect']}
    return [guided_rule(dets[v['detect']], grid, cov, yd, confirm=v['confirm'], prior=prior,
                        share=v['share'] or 0.95, ci_prunes=v['ci_counted'], prune_expect=v['prune_expect'])[0]
            for v in rule_family(prior is not None)]


def run_level(task):
    """All draws of one source at one noise level. Returns counts (see main)."""
    src, z_true, f2, target, factor, n, seed = task
    spec = read_spectrum(DATA / 'spectra' / f'{src}.csv')
    rng = np.random.default_rng(seed)
    priors = dict(PRIORS, **photoz_priors(src))             # the source's own photometric redshift, three versions
    photoz = [p for p in priors if p in PHOTOZ_STATES]
    c = dict(rule_ok=0, rule_wrong=0, rule_none=0, second_detected=0, y2=[])
    fam = {p: np.zeros((len(rule_family(priors[p] is not None)), 2), int) for p in FAMILY_PRIORS}   # correct, wrong
    curve = {p: np.zeros((len(ENGINE_THRESHOLDS), 2), int) for p in FAMILY_PRIORS}  # the engine, per threshold
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
        for p in FAMILY_PRIORS:
            for k, zg in enumerate(family_answers(yd, grid, cov, priors[p])):
                if zg is not None:
                    fam[p][k, 0 if correct(zg, z_true) else 1] += 1
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
                if m == 'engine' and p in curve:
                    curve[p][:, 0 if ok else 1] += conf >= np.array(ENGINE_THRESHOLDS)
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
    # the next two are moved out of the saved levels in main
    out['rule_family'] = {p: fam[p].tolist() for p in FAMILY_PRIORS}
    out['engine_curve'] = {p: curve[p].tolist() for p in FAMILY_PRIORS}
    return src, out


def family_summary(counts, curves, sources, observed):
    """Choose the rule settings to compare with, by the criterion stated at WRONG_CAP.

    counts[source][level][prior] = [[correct, wrong] per setting] over the draws of that source and level;
    curves[source][level][prior] = [[correct, wrong] per engine threshold].
    Returns (summary for the results file, rows for results/redshift_rule_family.csv)."""
    names = list(sources)
    n = sources[names[0]]['levels'][0]['draws']
    total = n * len(names) * len(TARGETS)
    out, rows = dict(settings=dict({k: list(v) for k, v in FAMILY.items()}, share_with_a_prior=list(FAMILY_SHARE)),
                     wrong_cap=WRONG_CAP, draws=total, engine_thresholds=list(ENGINE_THRESHOLDS), priors={}), []
    for p in FAMILY_PRIORS:
        fam = rule_family(p != 'uniform')
        a = np.array([[counts[s][k][p] for k in range(len(TARGETS))] for s in names])    # source, level, setting, 2
        e = np.array([[curves[s][k][p] for k in range(len(TARGETS))] for s in names])    # source, level, threshold, 2
        ok, wrong = a[..., 0].sum((0, 1)), a[..., 1].sum((0, 1))
        rate = a / n                                                                    # per source and level
        real = [[sum(r['rule_family'][p][j] == code for r in observed) for code in 'cwn'] for j in range(len(fam))]
        engine_wrong = int(round(sum(q[p]['engine']['confidently_wrong'] * q['draws']
                                     for s_ in sources.values() for q in s_['levels'])))
        chosen = {}
        for cap_name, cap in (('no_more_wrong_than_engine', engine_wrong), ('one_percent_wrong', WRONG_CAP * total)):
            eligible = np.where(wrong <= cap)[0]
            if not len(eligible):
                chosen[cap_name] = None
                continue
            j = int(eligible[np.argmax(ok[eligible])])
            chosen[cap_name] = dict(
                setting=fam[j], index=j, wrong_allowed=int(cap), total_correct=int(ok[j]), total_wrong=int(wrong[j]),
                levels=[dict(target=t, correct=round(float(rate[:, k, j, 0].mean()), 4),
                             wrong=round(float(rate[:, k, j, 1].mean()), 4)) for k, t in enumerate(TARGETS)],
                by_source={s: [round(float(x), 4) for x in rate[i, :, j, 0]] for i, s in enumerate(names)},
                real_spectra=dict(zip(('correct', 'wrong', 'no_answer'), real[j])))
        envelope = []
        for k, t in enumerate(TARGETS):
            okk = np.where(a[:, k, :, 1].sum(0) <= WRONG_CAP * n * len(names))[0]
            j = int(okk[np.argmax(rate[:, k, okk, 0].mean(0))])
            envelope.append(dict(target=t, correct=round(float(rate[:, k, j, 0].mean()), 4), index=j))
        # the engine at each threshold, and the rule with the most correct answers among those with no more wrong
        matched = []
        for i, thr in enumerate(ENGINE_THRESHOLDS):
            ec, ew = int(e[:, :, i, 0].sum()), int(e[:, :, i, 1].sum())
            eligible = np.where(wrong <= ew)[0]
            j = int(eligible[np.argmax(ok[eligible])]) if len(eligible) else None
            matched.append(dict(threshold=thr, engine_correct=ec, engine_wrong=ew,
                                engine_levels=[round(float((e[:, k, i, 0] / n).mean()), 4) for k in range(len(TARGETS))],
                                rule_index=j, rule_correct=None if j is None else int(ok[j]),
                                rule_wrong=None if j is None else int(wrong[j]),
                                rule_levels=None if j is None else
                                [round(float(rate[:, k, j, 0].mean()), 4) for k in range(len(TARGETS))],
                                rule_setting=None if j is None else fam[j]))
        out['priors'][p] = dict(n_settings=len(fam), engine_confident_wrong=engine_wrong, chosen=chosen,
                                envelope=envelope, engine_by_threshold=matched,
                                settings_with_no_more_wrong_than_engine=int((wrong <= engine_wrong).sum()),
                                settings_with_at_most_one_percent_wrong=int((wrong <= WRONG_CAP * total).sum()))
        for j, v in enumerate(fam):
            row = dict(prior=p, index=j, **v, total_correct=int(ok[j]), total_wrong=int(wrong[j]))
            for k, t in enumerate(TARGETS):
                row[f'correct_{t:g}'] = round(float(rate[:, k, j, 0].mean()), 4)
            for k, t in enumerate(TARGETS):
                row[f'wrong_{t:g}'] = round(float(rate[:, k, j, 1].mean()), 4)
            row.update(real_correct=real[j][0], real_wrong=real[j][1], real_no_answer=real[j][2])
            rows.append(row)
    return out, rows


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
        # every setting of the rule family: one letter each, c(orrect), w(rong) or n(o answer)
        fp = dict(uniform=None, photoz=photoz_priors(src).get('photoz'))
        row['rule_family'] = {p: ''.join(outcome(zz)[0] for zz in family_answers(yd, grid, cov, fp[p]))
                              for p in FAMILY_PRIORS if p == 'uniform' or fp[p] is not None}
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
    counts, curves = {src: [] for src in sources}, {src: [] for src in sources}
    with mp.get_context('fork').Pool(2) as pool:
        for src, out in pool.imap(run_level, tasks):
            counts[src].append(out.pop('rule_family'))
            curves[src].append(out.pop('engine_curve'))
            sources[src]['levels'].append(out)
    family, rows = family_summary(counts, curves, sources,
                                  [r for r in observed if len(r['rule_family']) == len(FAMILY_PRIORS)])
    out = dict(n_draws=a.n, seed=SEED, targets=TARGETS, confidence_threshold=CONF,
               note='noise added to the real spectra; second line = second-brightest detected line of the '
                    'published redshift; fractions are over draws, the summary is the mean over sources',
               sources=sources, summary=summarise(sources), rule_family=family, as_observed=observed)
    path.write_text(json.dumps(out, indent=1) + '\n')
    with open(RESULTS / 'redshift_rule_family.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
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
    print(f'Rule family, {family["draws"]} draws.')
    for p, f in family['priors'].items():
        print(f"  prior {p}, {f['n_settings']} settings: engine confident and wrong in {f['engine_confident_wrong']} draws; "
              f"{f['settings_with_no_more_wrong_than_engine']} settings are wrong no more often, "
              f"{f['settings_with_at_most_one_percent_wrong']} in at most 1%")
        for name, ch in f['chosen'].items():
            if ch:
                print(f"    {name}: {ch['setting']} correct {ch['total_correct']} wrong {ch['total_wrong']} | "
                      + ' '.join(f"{q['correct']:.3f}" for q in ch['levels']) + f" | real spectra {ch['real_spectra']}")
        print('    envelope: ' + ' '.join(f"{q['correct']:.3f}" for q in f['envelope']))
        print('    engine:   ' + ' '.join(f"{q[p]['engine']['confident_correct']:.3f}" for q in out['summary']))
        for q in f['engine_by_threshold']:
            print(f"    engine at {q['threshold']}: correct {q['engine_correct']} wrong {q['engine_wrong']} | best rule "
                  f"with no more wrong: correct {q['rule_correct']} wrong {q['rule_wrong']} {q['rule_setting']}")
    print(f'[{time.time() - t0:.0f}s] wrote {path.name} and redshift_rule_family.csv')
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
    plot_operating(o)


def _chosen(o, prior):
    """The rule setting the engine is compared with: the most correct answers among the settings of the family
    that return a wrong redshift no more often than the engine does with the same prior."""
    return o['rule_family']['priors'][prior]['chosen']['no_more_wrong_than_engine']


def plot_curves(o):
    """Mean over the sources against the S/N left to the second line. Solid: no prior. Dashed: survey prior
    (engine, cross-correlation) or the single rule setting used before the family was run."""
    rows = o['summary']
    x = [q['target'] for q in rows]
    eng, rule3, rule5, xc = COL[0], COL[1], COL[2], COL[3]
    best = _chosen(o, 'uniform')
    best_of = {q['target']: q['correct'] for q in best['levels']}
    panels = [
        ('A. Answers given: the published redshift, returned (engine: with confidence)', [
            ('engine, no prior', eng, 'o', '-', lambda q: q['uniform']['engine']['confident_correct']),
            ('engine, survey prior', eng, 'o', '--', lambda q: q['spt']['engine']['confident_correct']),
            (f"best of {o['rule_family']['priors']['uniform']['n_settings']} rule settings", rule3, 'D', '-',
             lambda q: best_of[q['target']]),
            ('one rule setting (5σ, then 3σ)', rule3, 'D', '--', lambda q: q['guided']['3']['correct']),
            ('every line at 5σ', rule5, 's', '-', lambda q: q['rule']['correct'])]),
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
        ax.legend(fontsize=8.3, frameon=False, loc=where, labelcolor=INK, handlelength=3.4)
        _style(ax)
    axs[0].set_ylabel('fraction of noise draws (mean of the sources)', fontsize=9, color=INK2)
    axs[0].axvline(5.0, color=INK2, lw=0.8, ls=':')
    axs[0].text(4.97, 0.03, "the detection\nrule's threshold", fontsize=8, color=INK2, va='bottom', ha='right')
    n_src, n = rows[0]['n_sources'], rows[0]['draws_per_source']
    fig.suptitle(f'Real multi-line spectra with noise added: {n_src} sources of 2016.1.00672.S, {n} noise draws '
                 'per source and level', fontsize=10.5, color=INK, x=0.012, ha='left')
    total, wrong = wrong_counts(o)
    wrong = dict(list(wrong.items())[:2] + [('best rule setting', best['total_wrong'])] + list(wrong.items())[2:])
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
        g = np.array([_chosen(o, 'uniform')['by_source'][s][k] for s in names])
        r = np.array([q['rule']['correct'] for q in lv])
        ax.hlines(y, np.minimum(e, g), np.maximum(e, g), color='#c3c2b7', lw=2, zorder=1)
        ax.plot(r, y, 's', color=COL[2], ms=7, mec=SURF, mew=1.2, zorder=2, label='detection rule: every line at 5σ')
        ax.plot(g, y, 'D', color=COL[1], ms=7, mec=SURF, mew=1.2, zorder=3,
                label=f"best of {o['rule_family']['priors']['uniform']['n_settings']} rule settings, no prior")
        ax.plot(ep, y, 'o', mfc=SURF, mec=COL[0], ms=7.5, mew=1.5, zorder=4, label='engine, survey prior')
        ax.plot(e, y, 'o', color=COL[0], ms=7.5, mec=SURF, mew=1.2, zorder=5, label='engine, no prior')
        ahead, behind = int((e > g).sum()), int((e < g).sum())
        ax.set_title(f'second line at {t:g}σ', fontsize=10, color=INK, loc='left')
        ax.set_xlabel(f'engine (no prior) ahead of the best rule\nfor {ahead} of {len(names)}, behind for {behind}',
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


def plot_operating(o, shown=30):
    """Answers against wrong answers over all draws: the engine at each confidence threshold and every rule
    setting, with no prior and with each source's 2020 photometric redshift."""
    path = RESULTS / 'redshift_rule_family.csv'
    if not path.exists():
        return
    with open(path, newline='') as fh:
        rows = list(csv.DictReader(fh))
    total = o['rule_family']['draws']
    fig, axs = plt.subplots(1, 2, figsize=(12.4, 5.4), facecolor=SURF, sharey=True)
    titles = dict(uniform='A. No redshift prior', photoz="B. Each source's 2020 photometric redshift as the prior")
    for ax, p in zip(axs, FAMILY_PRIORS):
        f = o['rule_family']['priors'][p]
        rr = [(int(r['total_wrong']), int(r['total_correct']) / total) for r in rows if r['prior'] == p]
        inside = [q for q in rr if q[0] <= shown]
        ax.plot([q[0] for q in inside], [q[1] for q in inside], 'D', color=COL[1], ms=5, alpha=0.55, mec='none',
                label=f"rule settings ({len(inside)} of {len(rr)}; the others give more than {shown} wrong answers)",
                zorder=2)
        ch = f['chosen']['no_more_wrong_than_engine']
        ax.plot([ch['total_wrong']], [ch['total_correct'] / total], 'D', mfc='none', mec=INK, ms=10, mew=1.5,
                label='the rule setting the engine is compared with', zorder=4)
        e = [q for q in f['engine_by_threshold'] if q['engine_wrong'] <= shown]
        ax.plot([q['engine_wrong'] for q in e], [q['engine_correct'] / total for q in e], '-o', color=COL[0], ms=7.5,
                lw=2, mec=SURF, mew=1.2, label='engine, at the confidence thresholds marked', zorder=3)
        for q in e:
            ax.annotate(f"{q['threshold']:g}", (q['engine_wrong'], q['engine_correct'] / total),
                        textcoords='offset points', xytext=(7, -11), fontsize=8.5, color=INK2)
        ax.set_title(titles[p], fontsize=10, color=INK, loc='left')
        ax.set_xlim(-1, shown + 1.5)
        ax.set_xlabel(f'draws in which a different redshift is returned, of {total:,}', fontsize=9, color=INK2)
        ax.grid(color=GRID_C, lw=0.8)
        ax.legend(fontsize=8.5, frameon=False, loc='lower right', labelcolor=INK)
        _style(ax)
    axs[0].set_ylabel('fraction of all draws in which the published redshift is returned', fontsize=9, color=INK2)
    fig.suptitle('Answers against wrong answers, over all noise levels: the engine and every rule setting',
                 fontsize=10.5, color=INK, x=0.012, ha='left')
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIG / 'redshift_operating_points.png', dpi=130, facecolor=SURF)
    plt.close(fig)


if __name__ == '__main__':
    main()
