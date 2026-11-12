#!/usr/bin/env python3
"""Check data/truth.csv and the cut-out centres against Reuter et al. 2020 (ApJ 902, 78), row by row.

The paper is not redistributed here. Download its PDF, then run from the repository root:
    python scripts/verify_truth_table.py Reuter_2020_ApJ_902_78.pdf
Needs `pdftotext` (poppler). Compares, for every source whose reference is Reuter et al. 2020:
  - Table E1: the redshift (to the last published digit) and the lines seen in the 3 mm scan;
  - Table A1: the ALMA position, against data/cutouts/2016.1.00672.S_targets.csv.
Exits with status 1 on any mismatch.
"""
import csv
import math
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def hms(h):
    a = [float(x) for x in h.split(':')]
    return 15 * (a[0] + a[1] / 60 + a[2] / 3600)


def dms(d):
    d = d.replace('−', '-')
    sgn = -1 if d.startswith('-') else 1
    a = [abs(float(x)) for x in d.lstrip('+-').split(':')]
    return sgn * (a[0] + a[1] / 60 + a[2] / 3600)


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    text = subprocess.run(['pdftotext', '-layout', sys.argv[1], '-'], capture_output=True, text=True).stdout
    with open(ROOT / 'data' / 'truth.csv', newline='') as fh:
        truth = {r['source']: r for r in csv.DictReader(fh) if 'Reuter' in r['reference']}
    with open(ROOT / 'data' / 'cutouts' / '2016.1.00672.S_targets.csv', newline='') as fh:
        targets = {r['target']: r for r in csv.DictReader(fh)}
    e1 = {m.group(1): (m.group(2), m.group(3).strip()) for m in re.finditer(
        r'(SPT\d{4}-\d{2})\s+\d\s+(\d\.\d+)\s*\(\d\)\s+[BC]\d(?:,\s*[BC]\d)*\s+[\d.]+\s+(.*)', text)}
    a1 = {m.group(1): (hms(m.group(2)), dms(m.group(3))) for m in re.finditer(
        r'(SPT\d{4}-\d{2})\s+SPT-?S\s*J[\d.\-]+\s+(\d+:\d+:[\d.]+)\s+([−+-]?\d+:\d+:[\d.]+)', text)}
    bad = 0
    for s, t in sorted(truth.items()):
        if s not in e1:
            print(f'{s}: not found in Table E1'); bad += 1; continue
        z, lines = e1[s]
        ok_z = abs(float(z) - float(t['z_published'])) < 1e-9
        seen = re.findall(r'(CO\(\d.\d\)|\[C I\]\(\d.\d\))', lines.split('  ')[0])
        print(f"{s}: z {z} vs {t['z_published']} {'ok' if ok_z else 'MISMATCH'}; 3 mm lines {', '.join(seen)}")
        bad += not ok_z
        if s in targets and s in a1:
            ra, dec = a1[s]
            r = targets[s]
            off = math.hypot((ra - float(r['ra_deg'])) * math.cos(math.radians(dec)), dec - float(r['dec_deg'])) * 3600
            print(f'   position offset {off:.2f}"' + ('' if off < 0.1 else '  MISMATCH'))
            bad += off >= 0.1
    print(f'{len(truth)} sources checked, {bad} mismatches')
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
