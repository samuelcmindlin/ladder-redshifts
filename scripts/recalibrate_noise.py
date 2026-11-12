#!/usr/bin/env python3
"""Apply the noise calibration to cube-route spectra that were extracted before it existed.

The spectra under data/spectra/ and the two SPT0027-50 Band 2 files were written with the extractor's annulus
noise, which is too high at the source position (docs/development-notes.md, item 20). realspec.spectrum_point()
and realspec.spectrum() now rescale the noise themselves (realspec.calibrate_noise), so a spectrum rebuilt from
its cut-outs needs nothing more. The cut-outs are not shipped, so this script applies the same calibration to
the saved files once: only the noise_mjy column changes, and the factors go to data/noise_calibration.csv
(noise_mjy = extractor noise x scale_applied).

It refuses to run on files that are already listed in data/noise_calibration.csv.

Usage (from the repository root):  python scripts/recalibrate_noise.py [--dry-run]
"""
import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402

from ladder_redshifts.batch_real import read_spectrum, write_spectrum  # noqa: E402
from ladder_redshifts.paths import DATA  # noqa: E402
from ladder_redshifts.realspec import calibrate_noise, record_noise_scales  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--dry-run', action='store_true', help='report the factors without writing anything')
    a = ap.parse_args()
    table = DATA / 'noise_calibration.csv'
    done = set()
    if table.exists():
        with open(table, newline='') as fh:
            done = {r['spectrum'] for r in csv.DictReader(fh)}
    files = sorted((DATA / 'spectra').glob('*.csv')) + sorted(DATA.glob('*_band2_spectrum*.csv'))
    seen = {}
    for f in files:
        label = f.relative_to(DATA).as_posix()
        if label in done:
            print(f'{label:52s} already calibrated; skipped')
            continue
        spec = read_spectrum(f)
        if any(w.get('kind', 'measured') != 'measured' for w in spec):
            print(f'{label:52s} not a cube-route spectrum; skipped')
            continue
        cal = calibrate_noise(spec)
        applied = np.array([w['noise_scale'] for w in cal])
        print(f'{label:52s} {len(cal):2d} windows, noise x {applied.min():.2f}-{applied.max():.2f} '
              f'(median {np.median(applied):.2f})')
        seen[label] = applied
        if not a.dry_run:
            write_spectrum(cal, f)
            record_noise_scales(label, cal, table)
    if not seen:
        print('nothing to do')


if __name__ == '__main__':
    main()
