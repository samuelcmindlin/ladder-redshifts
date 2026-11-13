"""Define the independent test sample (docs/independent-test-plan.md) and write its input tables.

The sample is every science target of the 3 mm scans of ALMA project 2015.1.00504.S (Cycle 3) that has a redshift
in Reuter et al. 2020 (ApJ 902, 78), Table F.1, and is not in data/truth.csv: 23 sources. None of their spectra
had been downloaded when the plan and this script were written.

    python scripts/independent_sample.py spt_dndz_tab_a1.txt spt_dndz_tab_f1.txt spt_dndz_tab_f1_kappa_wd03.txt

The three arguments are the tables the SPT-SMG collaboration posts at https://github.com/spt-smg/publicdata
(directory reuter2020_final_dndz; not copied here). Written:
    data/cutouts/2015.1.00504.S_targets.csv   positions (Table A.1) for the cut-outs and the extraction
    data/independent/truth.csv                published redshifts (Table F.1) and the lines of the 3 mm scan
    data/independent/photoz.csv               photometric redshifts, 2020 and 2024 (as data/photoz.csv)

LINES is Table E.1, column "Lines from 3 mm scans", read from the arXiv copy (v3). The script checks that every
listed line falls inside the scan at the published redshift. Standard library only.
"""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PROJECT = '2015.1.00504.S'
REST = {'CO(3-2)': 345.796, 'CO(4-3)': 461.041, 'CO(5-4)': 576.268, 'CO(6-5)': 691.473, 'CO(7-6)': 806.652,
        '[CI](1-0)': 492.161, '[CI](2-1)': 809.342}
SCAN = (84.0, 115.2)                    # GHz; the scans cover 84.2-114.9
LINES = {
    'SPT0020-51': ('CO(4-3)', 'CO(5-4)', '[CI](1-0)'),
    'SPT0106-64': ('CO(5-4)',),
    'SPT0109-47': ('CO(4-3)', '[CI](1-0)'),
    'SPT0112-55': ('CO(4-3)', '[CI](1-0)'),
    'SPT0136-63': ('CO(4-3)', 'CO(5-4)', '[CI](1-0)'),
    'SPT0147-64': ('CO(5-4)', '[CI](1-0)'),
    'SPT0150-59': ('CO(3-2)',),
    'SPT0155-62': ('CO(4-3)', 'CO(5-4)', '[CI](1-0)'),
    'SPT0202-61': ('CO(5-4)', 'CO(6-5)'),
    'SPT0245-63': ('CO(5-4)', 'CO(6-5)'),
    'SPT0311-58': ('CO(6-5)', 'CO(7-6)', '[CI](2-1)'),
    'SPT0348-62': ('CO(5-4)', 'CO(6-5)'),
    'SPT0516-59': ('CO(4-3)',),
    'SPT0544-40': ('CO(4-3)', 'CO(5-4)', '[CI](1-0)'),
    'SPT0552-42': ('CO(4-3)', 'CO(5-4)'),
    'SPT0553-50': ('CO(5-4)', 'CO(6-5)'),
    'SPT0555-62': ('CO(5-4)',),
    'SPT0604-64': ('CO(3-2)',),
    'SPT0611-55': ('CO(3-2)',),
    'SPT0625-58': ('CO(3-2)',),
    'SPT2037-65': ('CO(4-3)',),
    'SPT2048-55': ('CO(4-3)', 'CO(5-4)', '[CI](1-0)'),
    'SPT2101-60': ('CO(4-3)',),
}
# Two sources need a word. The table's lines for SPT0112-55 are from a deeper Cycle 7 scan (2019.1.00486.S), which
# is not part of this sample: the Cycle 3 scan showed none. SPT0150-59 had one line and no second observation in
# 2020: its redshift rested on the photometric redshift, as SPT0314-44's did.
SPECIAL = {'SPT0112-55': ('(0 lines in the Cycle 3 scan used here)', 'confirmed'),
           'SPT0150-59': (None, 'tentative')}


def table(path):
    return [line.split('\t') for line in Path(path).read_text().splitlines() if line.startswith('SPT')]


def degrees(text, hours=False):
    sign = -1 if text.strip().startswith('-') else 1
    a, b, c = (abs(float(x)) for x in text.strip().split(':'))
    return sign * (a + b / 60 + c / 3600) * (15 if hours else 1)


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    pos = {r[0].strip(): [x for x in r[1:] if x.strip()] for r in table(sys.argv[1])}
    f1 = {r.split()[0]: r.split() for r in Path(sys.argv[2]).read_text().splitlines() if r.startswith('SPT')}
    f24 = {r.split()[0]: r.split() for r in Path(sys.argv[3]).read_text().splitlines() if r.startswith('SPT')}
    with open(ROOT / 'data' / 'truth.csv', newline='') as fh:
        seen = {r['source'] for r in csv.DictReader(fh)}
    assert len(LINES) == 23 and not (set(LINES) & seen) and set(LINES) <= set(f1) and set(LINES) <= set(pos)
    targets, truth, photoz = [], [], []
    for s in sorted(LINES):
        z = float(f1[s][1])
        for line in LINES[s]:
            f = REST[line] / (1 + z)
            assert SCAN[0] <= f <= SCAN[1], f'{s}: {line} at {f:.1f} GHz is outside the scan'
        ra, dec = pos[s][-2], pos[s][-1]
        targets.append([s, f'{degrees(ra, True):.7f}', f'{degrees(dec):.7f}', 'Reuter+20 Table A1', 1, 20, 8])
        n = len(LINES[s])
        tag, status = SPECIAL.get(s, (None, 'confirmed'))
        tag = tag or ('(1 line)' if n == 1 else f'({n} lines)')
        truth.append([s, f1[s][1], 'spt', 'Reuter et al. 2020, Tables E.1 and F.1',
                      f"Table E.1, lines from 3 mm scans: {', '.join(LINES[s])} {tag}", status])
        photoz.append([s, f'{float(f1[s][2]):g}', f'{float(f1[s][3]):g}', f'{float(f24[s][2]):g}', f'{float(f24[s][3]):g}'])
    out = ROOT / 'data' / 'independent'
    out.mkdir(parents=True, exist_ok=True)
    for path, head, rows in (
            (ROOT / 'data' / 'cutouts' / f'{PROJECT}_targets.csv',
             ['target', 'ra_deg', 'dec_deg', 'position_from', 'priority', 'radius_arcsec', 'box_arcsec'], targets),
            (out / 'truth.csv', ['source', 'z_published', 'prior', 'reference', 'notes', 'status'], truth),
            (out / 'photoz.csv', ['source', 'z_phot', 'z_phot_err', 'z_phot_2024', 'z_phot_2024_err'], photoz)):
        with open(path, 'w', newline='') as fh:
            w = csv.writer(fh, lineterminator='\n')
            w.writerow(head)
            w.writerows(rows)
        print(f'wrote {path.relative_to(ROOT)} ({len(rows)} sources)')
    groups = {}
    for r in truth:
        groups.setdefault(r[4][r[4].rindex('('):], []).append(r[0])
    for k, v in sorted(groups.items()):
        print(f'  {k}: {len(v)}')


if __name__ == '__main__':
    main()
