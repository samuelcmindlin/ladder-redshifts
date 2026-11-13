#!/usr/bin/env python3
"""What the published High-z Sudoku code says about the same sources and setups.

The code of Bakx & Dannerbauer (2022), https://github.com/tjlcbakx/redshift-search-graphs, is not redistributed
here (it carries no licence file). Clone it, then run from the repository root:

    python scripts/compare_sudoku.py /path/to/redshift-search-graphs/RSG.py

Its function RSGquality(windows, redshifts) does not read a spectrum. For each redshift it counts the lines of
its list (CO up to J = 20 and the far-infrared fine-structure lines; [CI] only with includeCI=True) that fall
inside the windows, and classes the redshift as
  - no line;
  - one line: robust if a photometric redshift good to dz/(1+z) = 0.13 excludes the neighbouring CO
    identifications at 5 sigma, or if those would have shown more than one line;
  - two lines: robust unless a neighbouring solution would also show exactly two CO lines;
  - three or more lines.
The authors' "robust redshift coverage" (their tuningGraph.py) is two + more - non-robust two + robust one.
Every covered line is assumed to be detected: the metric describes a scan of unlimited depth.

This script calls RSGquality unmodified
  1. at the published redshift of every source that has a spectrum in data/spectra/, with that spectrum's windows;
  2. on the redshifts of the synthetic atlas (results/redshift_atlas_raw.npz), for its six setups,
and writes its answers next to this engine's in results/redshift_sudoku_comparison.json.

Two details of the call:
  - RSGquality returns fractions over the redshifts it is given, and builds its redshift grid up to 1.5 times
    the largest of them; neighbouring solutions beyond the grid are read at its last point. A sentinel redshift
    (Z_TOP = 12) is therefore added to every call, so that the grid reaches every neighbouring solution, and its
    own class is subtracted from the result.
  - the grid has 40,000 points (the function's default is 10,000 in its README, 1,000 in the authors' scripts),
    which places a redshift to 0.0005.
"""
import csv
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ladder_redshifts.batch_real import read_spectrum  # noqa: E402
from ladder_redshifts.ladder import TOL_LN  # noqa: E402
from ladder_redshifts.atlas import SETUPS  # noqa: E402

Z_TOP = 12.0
GRID = 40000
BRIGHT = 16.0            # atlas sources with CO(4-3) S/N at least this are compared with the unlimited-depth metric
NAMES = ['no_line', 'one_line', 'two_lines', 'more_lines', 'robust_one_line', 'non_robust_two_lines']


def load_rsg(path):
    spec = importlib.util.spec_from_file_location('RSG', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def merge(windows):
    """Union of overlapping or touching windows (RSGquality loops over windows for every grid point)."""
    out = []
    for a, b in sorted(windows):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [(a, b) for a, b in out]


class Sudoku:
    def __init__(self, rsg):
        self.rsg = rsg
        self._top = {}

    def fractions(self, windows, z, **kw):
        """RSGquality's six fractions over the redshifts z (the sentinel's own class removed)."""
        z = np.atleast_1d(np.asarray(z, float))
        assert z.max() < Z_TOP
        lo, hi = [a for a, b in windows], [b for a, b in windows]
        key = (tuple(windows), tuple(sorted(kw.items())))
        if key not in self._top:
            self._top[key] = self.rsg.RSGquality(lo, hi, np.array([Z_TOP]), lin_arr_size=GRID, **kw)
        both = self.rsg.RSGquality(lo, hi, np.append(z, Z_TOP), lin_arr_size=GRID, **kw)
        return (both * (len(z) + 1) - self._top[key]) / len(z)

    def classify(self, windows, z, **kw):
        f = self.fractions(windows, [z], **kw)
        assert np.all(np.abs(f - np.round(f)) < 1e-9), f
        no, one, two, more, rob1, bad2 = (int(round(x)) for x in f)
        if no:
            label, n = 'no line', 0
        elif one:
            label, n = ('one line, robust with a photometric redshift' if rob1 else 'one line, not robust'), 1
        elif two:
            label, n = ('two lines, degenerate' if bad2 else 'two lines, robust'), 2
        else:
            label, n = 'three or more lines', 3
        return dict(lines=n, label=label, robust=bool(two + more - bad2 + rob1))


def robust_coverage(f):
    return float(f[2] + f[3] - f[5] + f[4])


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    rsg_path = Path(sys.argv[1])
    sud = Sudoku(load_rsg(rsg_path))
    out = dict(code='https://github.com/tjlcbakx/redshift-search-graphs', rsg_sha256=hashlib.sha256(
        rsg_path.read_bytes()).hexdigest(), z_top=Z_TOP, grid_points=GRID,
        settings=dict(sigma_threshold=5, dzUncertainty=0.13, nr_of_CO_lines=20))

    # 1. the real sources
    with open(ROOT / 'data' / 'truth.csv', newline='') as fh:
        truth = {r['source']: r for r in csv.DictReader(fh)}
    batch = {r['source']: r for r in json.loads((ROOT / 'results' / 'redshift_real_batch.json').read_text())}
    rows = []
    for src in sorted(truth):
        path = ROOT / 'data' / 'spectra' / f'{src}.csv'
        if not path.exists() or src not in batch:
            continue
        z = float(truth[src]['z_published'])
        wins = merge([(float(w['freq'].min()), float(w['freq'].max())) for w in read_spectrum(path)])
        b = batch[src]
        rows.append(dict(
            source=src, z_published=z, windows_ghz=[(round(a, 3), round(b_, 3)) for a, b_ in wins],
            sudoku_co=sud.classify(wins, z, includeCI=False), sudoku_co_ci=sud.classify(wins, z, includeCI=True),
            engine=dict(z=b['z_engine'], confidence=b['confidence'], verdict=b['verdict']),
            repo_rule_solved=b['z_sudoku'] is not None, detections_5sigma=len(b['detections'])))
        r = rows[-1]
        print(f"{src:11s} z={z:.4f}  CO only: {r['sudoku_co']['label']:45s} with [CI]: {r['sudoku_co_ci']['label']:45s}"
              f" engine: {b['confidence']:.2f} {b['verdict']}")
    out['sources'] = rows

    # 2. the synthetic atlas
    d = np.load(ROOT / 'results' / 'redshift_atlas_raw.npz', allow_pickle=True)
    z, rho = d['z'], d['rho']
    bright = rho >= BRIGHT

    def engine(setup, est, sel):
        zz, cc = d[f'z|{setup}|{est}'], d[f'conf|{setup}|{est}']
        ok = np.abs(np.log1p(np.nan_to_num(zz, nan=-0.5)) - np.log1p(z)) < TOL_LN
        return round(float((ok & (cc >= 0.95))[sel].mean()), 4)

    atlas = {}
    for name, wins in SETUPS.items():
        wins = merge(wins)
        e = {}
        for label, sel in (('bright', bright), ('all', np.ones(len(z), bool))):
            s = {}
            for ci in (False, True):
                for dz in ((0.13, 0.20) if label == 'bright' else (0.13,)):
                    f = sud.fractions(wins, z[sel], includeCI=ci, dzUncertainty=dz)
                    s[f"{'co_ci' if ci else 'co'}_dz{dz:.2f}"] = dict(
                        {k: round(float(v), 4) for k, v in zip(NAMES, f)},
                        robust_lines_only=round(float(f[2] + f[3] - f[5]), 4),
                        robust_with_photoz=round(robust_coverage(f), 4))
            e[label] = dict(n=int(sel.sum()), sudoku=s, engine=dict(
                sled=engine(name, 'sled', sel), sled_photoz=engine(name, 'sled + photo-z', sel),
                repo_rule=engine(name, 'sudoku rule', sel)))
        atlas[name] = e
        s, g = e['bright']['sudoku'], e['bright']['engine']
        print(f"{name:18s} bright: Sudoku robust, lines only {s['co_dz0.13']['robust_lines_only']:.3f} "
              f"(with [CI] {s['co_ci_dz0.13']['robust_lines_only']:.3f}); with photo-z "
              f"{s['co_dz0.13']['robust_with_photoz']:.3f} ({s['co_ci_dz0.13']['robust_with_photoz']:.3f}) | engine "
              f"{g['sled']:.3f}, with photo-z {g['sled_photoz']:.3f}; repo rule {g['repo_rule']:.3f}")
    out['atlas'] = dict(bright_snr=BRIGHT, setups=atlas)
    (ROOT / 'results' / 'redshift_sudoku_comparison.json').write_text(json.dumps(out, indent=1) + '\n')
    print('wrote results/redshift_sudoku_comparison.json')


if __name__ == '__main__':
    main()
