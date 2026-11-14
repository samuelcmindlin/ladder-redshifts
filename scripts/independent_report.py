"""Report on the independent test: the outcomes fixed in docs/independent-test-plan.md, from the results of the
single run in results/independent/.

    LADDER_REDSHIFTS_SAMPLE=independent python scripts/independent_report.py

Reads results/independent/redshift_real_batch.json, redshift_photoz.json, redshift_degradation.json and
redshift_rule_family.csv (the last two from scripts/independent_noise_test.py), and the development sample's
results/redshift_degradation.json. Writes results/independent/report.json and figures/independent/summary.png.

The plan quotes two rule settings chosen on the development draws. The noise-added results keep per-source
figures only for settings chosen on the sample at hand, so the per-source figures of those two settings are
recomputed here on the same noise draws (same seeds), and checked against the stored means over sources.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import binomtest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ladder_redshifts import degrade as D  # noqa: E402
from ladder_redshifts.paths import ROOT, DATA, RESULTS, FIG, SAMPLE  # noqa: E402
from ladder_redshifts.batch_real import read_spectrum, to_grid_auto  # noqa: E402
from ladder_redshifts.realspec import run_engine, degrade  # noqa: E402
from ladder_redshifts.ladder import matched_filter, detection_index, prior_spt, BANK, TOL_LN  # noqa: E402
from ladder_redshifts.photoz import photoz_priors, upper_bound, CONF  # noqa: E402

# the rule settings chosen on the development draws (docs/independent-test-plan.md)
CHOSEN = dict(uniform=dict(detect=4.5, confirm=2.5, ci_counted=False, prune_expect=5.0, share=None),
              photoz=dict(detect=4.5, confirm=3.0, ci_counted=False, prune_expect=8.0, share=0.95))
LEVELS = (4.0, 3.5, 3.0)
INK, INK2, GRID_C, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
BLUE = '#2a78d6'


def rule_per_source(o):
    """Fraction of draws in which each chosen rule setting returns the published redshift, per source and level,
    on the noise draws of results/independent/redshift_degradation.json (regenerated from their seeds)."""
    kdet = detection_index(BANK)
    out = {p: {} for p in CHOSEN}
    for i, (src, s) in enumerate(o['sources'].items()):
        z, spec = s['z_published'], read_spectrum(DATA / 'spectra' / f'{src}.csv')
        r = run_engine(*to_grid_auto(spec), prior=prior_spt)
        snr2 = D.published_lines(r['detections'], z)[1][2]
        priors = dict(uniform=None, photoz=photoz_priors(src)['photoz'])
        for p in CHOSEN:
            out[p][src] = []
        for k, target in enumerate(o['targets']):
            rng = np.random.default_rng([o['seed'], i, k])
            ok = dict.fromkeys(CHOSEN, 0)
            for _ in range(o['n_draws']):
                grid, d, _, cov = to_grid_auto(degrade(spec, snr2 / target, rng))
                yd = matched_filter(d, BANK, cov)[kdet]
                for p, v in CHOSEN.items():
                    dets = D.detect(yd, cov, v['detect'])
                    zg = D.guided_rule(dets, grid, cov, yd, confirm=v['confirm'], prior=priors[p],
                                       share=v['share'] or 0.95, ci_prunes=v['ci_counted'],
                                       prune_expect=v['prune_expect'])[0]
                    ok[p] += D.correct(zg, z)
            for p in CHOSEN:
                out[p][src].append(ok[p] / o['n_draws'])
    return out


def sign_test(ahead, behind):
    n = ahead + behind
    return None if n == 0 else round(float(binomtest(min(ahead, behind), n, 0.5).pvalue), 3)


def main():
    assert SAMPLE == 'independent', 'run with LADDER_REDSHIFTS_SAMPLE=independent'
    batch = {r['source']: r for r in json.loads((RESULTS / 'redshift_real_batch.json').read_text())}
    ph = json.loads((RESULTS / 'redshift_photoz.json').read_text())
    o = json.loads((RESULTS / 'redshift_degradation.json').read_text())
    dev = json.loads((ROOT / 'results' / 'redshift_degradation.json').read_text())
    with open(RESULTS / 'redshift_rule_family.csv', newline='') as fh:
        fam_rows = list(csv.DictReader(fh))
    obs = {r['source']: r for r in o['as_observed']}
    idx = {p: D.rule_family(p != 'uniform').index(v) for p, v in CHOSEN.items()}
    letters = {'c': 'correct', 'w': 'wrong', 'n': 'none'}

    # per source, as observed
    rows = []
    for r in ph['sources']:
        s, st = r['source'], r['states']
        rows.append(dict(
            source=s, z_published=r['z_published'], group=r['lines'], status=batch[s]['verdict'],
            detections_5sigma=r['detections'],
            engine={k: dict(z=e['z'], confidence=e['confidence'], p_published=e['p_published'],
                            correct=bool(abs(np.log1p(e['z']) - np.log1p(r['z_published'])) < TOL_LN))
                    for k, e in (('uniform', r['reference']['uniform']), ('spt', r['reference']['spt']),
                                 ('photoz', st['photoz']['engine']), ('photoz_2024', st['photoz_2024']['engine']))},
            photoz=dict(centre=st['photoz']['centre'], width=st['photoz']['width'], offset=-st['photoz']['pull']),
            rule={p: letters[obs[s]['rule_family'][p][idx[p]]] for p in CHOSEN},
            rule_5sigma='none' if batch[s]['z_sudoku'] is None else
            'correct' if abs(np.log1p(batch[s]['z_sudoku']) - np.log1p(r['z_published'])) < TOL_LN else 'wrong',
            xcorr={p: ('correct' if obs[s]['xcorr'][p]['best_is_published'] else 'wrong') for p in ('uniform', 'spt')},
            xcorr_photoz=st['photoz']['xcorr']['outcome']))

    def count(group, f):
        return int(sum(bool(f(r)) for r in rows if group is None or r['group'] == group))
    conf_ok = lambda p: (lambda r: r['engine'][p]['confidence'] >= CONF and r['engine'][p]['correct'])
    conf_bad = lambda p: (lambda r: r['engine'][p]['confidence'] >= CONF and not r['engine'][p]['correct'])
    dev_conf = dict(spt=8, photoz=15)               # confident answers on the development spectra, none wrong
    primary = {}
    for p in ('spt', 'photoz'):
        n, k = count(None, conf_ok(p)) + count(None, conf_bad(p)), count(None, conf_bad(p))
        primary[p] = dict(confident_answers=n, wrong=k, wrong_rate_upper_95=round(upper_bound(k, n), 3),
                          pooled_confident_answers=n + dev_conf[p], pooled_wrong=k,
                          pooled_wrong_rate_upper_95=round(upper_bound(k, n + dev_conf[p]), 3))
    groups = {}
    for g in ('several', 'one', 'none'):
        groups[g] = dict(
            n=count(g, lambda r: True),
            engine={p: dict(confident_correct=count(g, conf_ok(p)), confident_wrong=count(g, conf_bad(p)),
                            best_is_published=count(g, lambda r: r['engine'][p]['correct']))
                    for p in ('uniform', 'spt', 'photoz', 'photoz_2024')},
            chosen_rule={p: {v: count(g, lambda r: r['rule'][p] == v) for v in ('correct', 'wrong', 'none')}
                         for p in CHOSEN},
            rule_5sigma={v: count(g, lambda r: r['rule_5sigma'] == v) for v in ('correct', 'wrong', 'none')},
            xcorr_best_is_published=dict(uniform=count(g, lambda r: r['xcorr']['uniform'] == 'correct'),
                                         photoz=count(g, lambda r: r['xcorr_photoz'] == 'correct')))
    # sources where the engine (confident and correct) and the chosen rule (correct) differ, as observed
    differ = {}
    for p, e in (('uniform', 'uniform'), ('photoz', 'photoz')):
        a = [r['source'] for r in rows if conf_ok(e)(r) and r['rule'][p] != 'correct']
        b = [r['source'] for r in rows if r['rule'][p] == 'correct' and not conf_ok(e)(r)]
        differ[p] = dict(engine_only=a, rule_only=b, sign_test_p=sign_test(len(a), len(b)),
                         rule_wrong=[r['source'] for r in rows if r['rule'][p] == 'wrong'])

    # noise added, the settings chosen on the development draws
    per = rule_per_source(o)
    T = o['targets']
    noise = dict(sources=list(o['sources']), draws=o['rule_family']['draws'], priors={})
    for p in CHOSEN:
        row = next(r for r in fam_rows if r['prior'] == p and int(r['index']) == idx[p])
        mean = [float(np.mean([per[p][s][k] for s in o['sources']])) for k in range(len(T))]
        assert all(abs(mean[k] - float(row[f'correct_{t:g}'])) < 1e-4 for k, t in enumerate(T)), 'draws differ'
        eng = {s: [q[p]['engine']['confident_correct'] for q in v['levels']] for s, v in o['sources'].items()}
        dch = dev['rule_family']['priors'][p]['chosen']['no_more_wrong_than_engine']
        assert dch['setting'] == CHOSEN[p]
        lev = []
        for t in LEVELS:
            k, kd = T.index(t), dev['targets'].index(t)
            ahead = int(sum(eng[s][k] > per[p][s][k] for s in eng))
            behind = int(sum(eng[s][k] < per[p][s][k] for s in eng))
            da = int(sum(v['levels'][kd][p]['engine']['confident_correct'] > dch['by_source'][s][kd]
                         for s, v in dev['sources'].items()))
            db = int(sum(v['levels'][kd][p]['engine']['confident_correct'] < dch['by_source'][s][kd]
                         for s, v in dev['sources'].items()))
            lev.append(dict(target=t, engine=round(float(np.mean([eng[s][k] for s in eng])), 3), rule=round(mean[k], 3),
                            engine_ahead=ahead, engine_behind=behind, development_ahead=da, development_behind=db,
                            pooled_sign_test_p=sign_test(ahead + da, behind + db)))
        fo = o['rule_family']['priors'][p]
        noise['priors'][p] = dict(
            setting=CHOSEN[p], levels=lev, rule_wrong=int(row['total_wrong']), rule_correct=int(row['total_correct']),
            engine_wrong=fo['engine_confident_wrong'],
            engine_correct=next(q['engine_correct'] for q in fo['engine_by_threshold'] if q['threshold'] == CONF),
            by_source={s: dict(engine=[round(float(eng[s][T.index(t)]), 2) for t in LEVELS],
                               rule=[round(float(per[p][s][T.index(t)]), 2) for t in LEVELS]) for s in eng},
            engine_by_threshold=[{k: q[k] for k in ('threshold', 'engine_correct', 'engine_wrong', 'rule_correct',
                                                    'rule_wrong')} for q in fo['engine_by_threshold']])
    report = dict(plan='docs/independent-test-plan.md', code='tag v0.2.0', sample=len(rows), primary=primary,
                  groups=groups, as_observed_engine_against_chosen_rule=differ, noise_added=noise, sources=rows)
    (RESULTS / 'report.json').write_text(json.dumps(report, indent=1) + '\n')
    plot(rows, groups)
    print(json.dumps(dict(primary=primary, groups=groups, differ=differ), indent=1))
    for p, v in noise['priors'].items():
        print(p, json.dumps(v['levels']), 'engine wrong', v['engine_wrong'], 'rule wrong', v['rule_wrong'],
              'engine correct', v['engine_correct'], 'rule correct', v['rule_correct'])
        print('   ', json.dumps(v['by_source']))


def plot(rows, groups):
    """One row per source: the engine's probability for the published redshift under two priors, and what the
    rule chosen beforehand answers."""
    order = [r for g in ('several', 'one', 'none') for r in sorted((x for x in rows if x['group'] == g),
                                                                    key=lambda x: x['z_published'])]
    y, gaps, last = [], 0, None
    for r in order:
        if last is not None and r['group'] != last:
            gaps += 1
        y.append(len(y) + gaps * 0.8); last = r['group']
    y = -np.array(y)
    logit = lambda p: np.log10(np.clip(p, 2e-3, 0.9995) / (1 - np.clip(p, 2e-3, 0.9995)))
    fig, ax = plt.subplots(figsize=(11.6, 8.2), facecolor=SURF)
    for yy in y:
        ax.axhline(yy, color=GRID_C, lw=0.8, zorder=0)
    ax.plot([logit(r['engine']['spt']['p_published']) for r in order], y + 0.13, 'o', mfc=SURF, mec=INK2, ms=7.5, mew=1.5,
            label='engine, survey redshift distribution as prior', zorder=2, clip_on=False)
    ax.plot([logit(r['engine']['photoz']['p_published']) for r in order], y - 0.13, 'o', color=BLUE, ms=7.5, mec=SURF,
            mew=1.2, label="engine, the source's 2020 photometric redshift as prior", zorder=3, clip_on=False)
    ticks = [0.01, 0.1, 0.5, 0.9, 0.99, 0.999]
    ax.set_xticks([logit(t) for t in ticks]); ax.set_xticklabels([f'{t:g}' for t in ticks])
    x1 = logit(0.9995) + 0.15
    ax.set_xlim(logit(2e-3) - 0.1, x1 + 2.6)
    ax.axvline(logit(CONF), color=INK2, lw=0.9, ls=':')
    ax.text(logit(CONF) - 0.06, y[0] + 0.95, 'confident from 0.95', fontsize=8, color=INK2, ha='right', va='center')
    word = dict(correct='correct', wrong='WRONG', none='no answer')
    for k, (p, head) in enumerate((('uniform', 'rule, no prior'), ('photoz', 'rule, photometric prior'))):
        xc = x1 + 0.55 + 1.25 * k
        ax.text(xc, y[0] + 0.95, head, fontsize=8, color=INK2, ha='center', va='center')
        for yy, r in zip(y, order):
            v = r['rule'][p]
            ax.text(xc, yy, word[v], fontsize=8, ha='center', va='center', color=INK if v != 'none' else '#a3a29a',
                    fontweight='bold' if v == 'wrong' else 'normal')
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['source']}  z = {r['z_published']:.3f}" for r in order], fontsize=8.5, color=INK)
    names = dict(several='two or three lines published from the 3 mm scan',
                 one='one line published from the 3 mm scan', none='no line in the scan used')
    last = None
    for yy, r in zip(y, order):
        if r['group'] != last:
            ax.text(logit(2e-3) - 0.1, yy + 0.62, names[r['group']], fontsize=8.5, color=INK2, ha='left', va='center',
                    style='italic')
            last = r['group']
    ax.set_ylim(y[-1] - 0.7, y[0] + 1.3)
    ax.set_xlabel("engine's probability for the published redshift (logit scale)", fontsize=9, color=INK2, x=0.36)
    ax.tick_params(colors=INK2, labelsize=8.5)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.set_facecolor(SURF)
    ax.legend(fontsize=8.5, frameon=False, loc='upper center', bbox_to_anchor=(0.36, -0.07), ncol=2, labelcolor=INK)
    g = groups
    wrong = sum(v['engine'][p]['confident_wrong'] for v in g.values() for p in ('spt', 'photoz'))
    fig.suptitle('23 sources the engine was not developed on (ALMA 2015.1.00504.S), run once with the code of tag v0.2.0.\n'
                 f"Confident and correct: {g['several']['engine']['spt']['confident_correct']} of {g['several']['n']} "
                 f"multi-line sources with the survey prior and {g['several']['engine']['photoz']['confident_correct']} "
                 f"with the photometric prior;\n{g['one']['engine']['photoz']['confident_correct']} of {g['one']['n']} "
                 f"one-line sources with the photometric prior. Confident and wrong: {wrong or 'none'}.",
                 fontsize=10.2, color=INK, x=0.012, y=0.985, ha='left')
    fig.subplots_adjust(left=0.19, right=0.985, top=0.875, bottom=0.11)
    fig.savefig(FIG / 'summary.png', dpi=130, facecolor=SURF)
    plt.close(fig)


if __name__ == '__main__':
    main()
