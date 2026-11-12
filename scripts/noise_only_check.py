#!/usr/bin/env python3
"""How does the engine behave on spectra that contain only noise?

Allowing each line a velocity offset (ladder.DV_MAX) lets noise fit a little better, because the largest of
several neighbouring filter outputs is positive on average. This script measures the price: the confidence the
engine reports on pure-noise spectra, with the tolerance and with exact line positions (dv_max = 0), for a
Band 3 scan and for Band 3 + Band 4 scans, under a uniform and an SPT-like redshift prior.

Writes results/redshift_noise_only.json. Usage (from the repository root):
    python scripts/noise_only_check.py [N_spectra]      (default 600; about 5 minutes)
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402

from ladder_redshifts.ladder import DV_MAX, Engine, Grid, prior_spt, templates_grid  # noqa: E402
from ladder_redshifts.paths import RESULTS  # noqa: E402

SETUPS = {'B3 scan': [(84.2, 114.9)], 'B3 + B4 scans': [(84.2, 114.9), (125.0, 163.0)]}
PRIORS = {'uniform': None, 'SPT-like': prior_spt}


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 600
    seed = 99
    rng = np.random.default_rng(seed)
    grid = Grid(66.0, 164.0)
    T, lw = templates_grid()
    out = dict(n_spectra=n, seed=seed, dv_max_kms=DV_MAX, setups=SETUPS, results={})
    for name, wins in SETUPS.items():
        cov = grid.mask(wins)
        engines = {dv: Engine(grid, cov, T, logw=lw, dv_max=dv) for dv in (DV_MAX, 0.0)}
        conf = {(dv, p): [] for dv in engines for p in PRIORS}
        for _ in range(n):
            d = rng.normal(0, 1, grid.npix)
            for dv, e in engines.items():
                L = e.loglike(e.filter(d))
                for p, pr in PRIORS.items():
                    conf[(dv, p)].append(e.posterior(L, prior=pr)[1])
        out['results'][name] = {}
        for (dv, p), v in conf.items():
            v = np.array(v)
            row = dict(median_confidence=float(np.median(v)), frac_above_0p5=float(np.mean(v > 0.5)),
                       frac_above_0p8=float(np.mean(v > 0.8)), frac_confident=float(np.mean(v >= 0.95)),
                       max_confidence=float(v.max()))
            key = ('with tolerance' if dv > 0 else 'exact positions') + f', {p} prior'
            out['results'][name][key] = row
            print(f"{name:14s} {key:34s} median {row['median_confidence']:.3f}  >0.5: {row['frac_above_0p5']:.3f}  "
                  f">0.8: {row['frac_above_0p8']:.3f}  >=0.95: {row['frac_confident']:.4f}  max {row['max_confidence']:.3f}")
    (RESULTS / 'redshift_noise_only.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
