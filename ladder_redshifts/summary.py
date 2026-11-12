"""Headline figure for the real-data test: the 18 sources of ALMA project 2016.1.00672.S with published redshifts.

For each source: the confidence the SLED likelihood (this engine) and the equal-weight sum give their best redshift,
and whether that redshift is the published one; plus whether the simplified Sudoku-style detection rule finds a
unique solution. Sources are grouped by how many lines the published 3 mm scan showed (Reuter et al. 2020).

Reads results/redshift_real_batch.json (written by batch_real) and data/truth.csv.
Writes figures/summary_real.png.
Usage: python -m ladder_redshifts.summary
"""
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from .paths import DATA, FIG, RESULTS
from .ladder import TOL_LN

INK, INK2, GRID_C, SURF, SHADE = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb', '#efeee9'
ENGINE, FLAT = '#2a78d6', '#eb6834'
THRESH = 0.95


def load():
    with open(DATA / 'truth.csv', newline='') as fh:
        truth = {r['source']: r for r in csv.DictReader(fh)}
    rows = json.loads((RESULTS / 'redshift_real_batch.json').read_text())
    groups = {'multi': [], 'single': []}
    for r in rows:
        t = truth.get(r['source'])
        if not t or r.get('z_engine') is None:
            continue
        notes = t['notes']
        g = 'multi' if ('(2 lines)' in notes or '(3 lines)' in notes) else \
            'single' if ('(1 line)' in notes or ' only;' in notes) else None
        if g is None:
            continue
        zt = float(t['z_published'])
        ok = lambda z: z is not None and abs(np.log1p(z) - np.log1p(zt)) < TOL_LN
        groups[g].append(dict(src=r['source'], tentative=t['status'] == 'tentative',
                              c_eng=r['confidence'], ok_eng=ok(r['z_engine']),
                              c_flat=r['conf_flat'], ok_flat=ok(r['z_flat']),
                              sudoku=None if r['z_sudoku'] is None else ok(r['z_sudoku'])))
    for g in groups.values():
        g.sort(key=lambda d: (-d['c_eng'], d['src']))
    return groups


def panel(ax, items, title):
    n = len(items)
    ax.axvspan(THRESH, 1.04, color=SHADE, lw=0, zorder=0)
    ax.axvline(THRESH, color=INK2, lw=0.8, ls=(0, (3, 3)), zorder=1)
    for i, d in enumerate(items):
        y = n - 1 - i
        ax.axhline(y, color=GRID_C, lw=0.6, zorder=1)
        for c, ok, col, dy, mk in ((d['c_flat'], d['ok_flat'], FLAT, -0.17, 'D'),
                                   (d['c_eng'], d['ok_eng'], ENGINE, 0.17, 'o')):
            if ok:
                ax.plot(c, y + dy, marker=mk, ms=8 if mk == 'o' else 7, color=col, mec=SURF, mew=1.2, ls='',
                        zorder=3)
            else:
                ax.plot(c, y + dy, marker='X', ms=9, color=col, mec=SURF, mew=1.0, ls='', zorder=3)
        lab = d['src'] + (' *' if d['tentative'] else '')
        ax.text(-0.02, y, lab, ha='right', va='center', fontsize=8.5, color=INK, transform=ax.get_yaxis_transform())
        s = d['sudoku']
        ax.text(1.02, y, 'solved' if s else ('wrong' if s is False else '–'), ha='left', va='center', fontsize=8,
                color=INK if s else INK2, transform=ax.get_yaxis_transform())
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlim(-0.02, 1.04)
    ax.set_yticks([])
    ax.set_title(title, fontsize=9.5, color=INK, loc='left', pad=6)
    ax.set_facecolor(SURF)
    for sp in ('top', 'right', 'left'):
        ax.spines[sp].set_visible(False)
    ax.spines['bottom'].set_color(INK2)
    ax.tick_params(axis='x', colors=INK2, labelsize=8)


def main():
    g = load()
    nm, ns = len(g['multi']), len(g['single'])
    fig, axs = plt.subplots(2, 1, figsize=(8.2, 0.34 * (nm + ns) + 2.4), facecolor=SURF, sharex=True,
                            gridspec_kw=dict(height_ratios=[nm + 1.2, ns + 1.2]))
    n_ok = sum(d['ok_eng'] and d['c_eng'] >= THRESH for d in g['multi'])
    n_top = sum(d['ok_eng'] for d in g['single'])
    panel(axs[0], g['multi'], f'Two or three lines in the scan ({nm} sources): engine confident and correct for {n_ok}')
    panel(axs[1], g['single'], f'One line in the scan ({ns} sources): engine never confident; '
                               f'published z ranked first for {n_top}')
    axs[0].text(1.02, 1.0, 'detection\nrule', ha='left', va='bottom', fontsize=7.5, color=INK2,
                transform=axs[0].transAxes)
    axs[0].text((THRESH + 1.04) / 2, 1.0, 'confident\n(≥ 0.95)', ha='center', va='bottom', fontsize=7.5, color=INK2,
                transform=axs[0].get_xaxis_transform())
    axs[1].set_xlabel('confidence score for the best redshift', fontsize=8.5, color=INK2)
    handles = [Line2D([], [], marker='o', ms=8, color=ENGINE, mec=SURF, ls='', label='SLED likelihood (this engine)'),
               Line2D([], [], marker='D', ms=7, color=FLAT, mec=SURF, ls='', label='equal-weight sum'),
               Line2D([], [], marker='o', ms=8, color=INK2, mec=SURF, ls='', label='best z = published'),
               Line2D([], [], marker='X', ms=9, color=INK2, mec=SURF, ls='', label='best z ≠ published')]
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.01, 0.995), ncol=4, frameon=False, fontsize=8,
               labelcolor=INK, handletextpad=0.3, columnspacing=1.2)
    star = '; * tentative published identification' if any(d['tentative'] for grp in g.values() for d in grp) else ''
    fig.text(0.01, 0.012, 'ALMA project 2016.1.00672.S, Band 3 scan, 84–115 GHz; published redshifts from Reuter et al. '
             f'2020{star}.\nBaselines are simplified re-implementations: equal-weight '
             'sum of line S/N; Sudoku-style detection rule (unique solution from lines ≥ 5σ), not the authors\' code.', fontsize=6.8, color=INK2, ha='left',
             va='bottom', linespacing=1.4)
    fig.subplots_adjust(left=0.17, right=0.89, top=0.9, bottom=0.115, hspace=0.26)
    out = FIG / 'summary_real.png'
    fig.savefig(out, dpi=150, facecolor=SURF)
    print(f'wrote {out.relative_to(FIG.parent)}')


if __name__ == '__main__':
    main()
