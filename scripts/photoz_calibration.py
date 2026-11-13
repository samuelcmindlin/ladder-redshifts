"""Published photometric redshifts of the SPT sample: copy the rows this repository uses and measure how
well the published values and their error bars match the spectroscopic redshifts.

The values are Table F.1 of Reuter et al. 2020 (ApJ 902, 78). The collaboration posts the table in
machine-readable form at https://github.com/spt-smg/publicdata (directory reuter2020_final_dndz; no licence
file, so the tables are not copied into this repository). Download two files from there and pass their paths:

    spt_dndz_tab_f1.txt             the table as published, except for the SPT0452-50 row (see below)
    spt_dndz_tab_f1_kappa_wd03.txt  the 2024 revision, with photometric redshifts refitted after the redshift
                                    of SPT0452-50 was corrected

    python scripts/photoz_calibration.py spt_dndz_tab_f1.txt spt_dndz_tab_f1_kappa_wd03.txt [--write]

Without --write the script checks that data/photoz.csv agrees with the two tables; with --write it creates
that file. Either way it writes results/redshift_photoz_calibration.json.

SPT0452-50: the journal table lists z_spec = 2.010 and z_phot = 4.4 +/- 0.7. The posted file keeps that row
as a comment and adds a recomputed row for the corrected z_spec = 5.016. The journal photometric redshift is
taken from the commented row; the spectroscopic redshift is the corrected one.

Standard library only, so it runs without the package installed.
"""
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PHOTOZ = ROOT / 'data' / 'photoz.csv'
TRUTH = ROOT / 'data' / 'truth.csv'
OUT = ROOT / 'results' / 'redshift_photoz_calibration.json'
COLUMNS = ['source', 'z_phot', 'z_phot_err', 'z_phot_2024', 'z_phot_2024_err']
JOURNAL_ROW_KEPT_AS_COMMENT = 'SPT0452-50'


def read_table(path):
    """{source: (z_spec, z_phot, err)} from a posted table, plus the commented-out journal rows."""
    rows, commented = {}, {}
    for line in Path(path).read_text().splitlines():
        t = line.split()
        if len(t) < 4:
            continue
        name = t[0].lstrip('#')
        if not name.startswith('SPT'):
            continue
        try:
            vals = tuple(float(x) for x in t[1:4])
        except ValueError:
            continue
        (commented if t[0].startswith('#') else rows)[name] = vals
    return rows, commented


def stats(rows):
    """rows: list of (z_spec, z_phot, err). Pull = (z_spec - z_phot) / err."""
    n = len(rows)
    pull = [(zs - zp) / e for zs, zp, e in rows]
    dz = [(zp - zs) / (1 + zs) for zs, zp, _ in rows]
    mean = sum(pull) / n
    sd = math.sqrt(sum((p - mean) ** 2 for p in pull) / (n - 1))
    med = sorted(dz)[n // 2] if n % 2 else 0.5 * (sorted(dz)[n // 2 - 1] + sorted(dz)[n // 2])
    dev = sorted(abs(d - med) for d in dz)
    mad = dev[n // 2] if n % 2 else 0.5 * (dev[n // 2 - 1] + dev[n // 2])
    mdz = sum(dz) / n
    return dict(n=n, pull_mean=round(mean, 3), pull_sd=round(sd, 3),
                pull_rms=round(math.sqrt(sum(p * p for p in pull) / n), 3),
                within_1_sigma=round(sum(abs(p) <= 1 for p in pull) / n, 3),
                within_2_sigma=round(sum(abs(p) <= 2 for p in pull) / n, 3),
                largest_pull=round(max(pull, key=abs), 2),
                dz_over_1pz_mean=round(mdz, 3), dz_over_1pz_median=round(med, 3),
                dz_over_1pz_sd=round(math.sqrt(sum((d - mdz) ** 2 for d in dz) / (n - 1)), 3),
                dz_over_1pz_nmad=round(1.4826 * mad, 3))


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    if len(args) != 2:
        sys.exit(__doc__)
    posted, commented = read_table(args[0])
    revised, _ = read_table(args[1])
    if len(posted) != 81 or set(posted) != set(revised) or JOURNAL_ROW_KEPT_AS_COMMENT not in commented:
        sys.exit('The tables do not look like the two files named in the docstring.')
    journal = dict(posted)
    zs_now = posted[JOURNAL_ROW_KEPT_AS_COMMENT][0]
    journal[JOURNAL_ROW_KEPT_AS_COMMENT] = (zs_now,) + commented[JOURNAL_ROW_KEPT_AS_COMMENT][1:]

    with open(TRUTH, newline='') as fh:
        sample = sorted(r['source'] for r in csv.DictReader(fh) if r['source'] in journal)
    want = [[s, f'{journal[s][1]:g}', f'{journal[s][2]:g}', f'{revised[s][1]:g}', f'{revised[s][2]:g}'] for s in sample]
    if '--write' in sys.argv:
        with open(PHOTOZ, 'w', newline='') as fh:
            w = csv.writer(fh, lineterminator='\n')
            w.writerow(COLUMNS)
            w.writerows(want)
        print(f'wrote {PHOTOZ} ({len(want)} sources)')
    else:
        with open(PHOTOZ, newline='') as fh:
            have = [[r[c] for c in COLUMNS] for r in csv.DictReader(fh)]
        if have != want:
            sys.exit(f'{PHOTOZ} does not agree with the tables')
        print(f'{PHOTOZ} agrees with the tables ({len(want)} sources)')

    outside = [s for s in journal if s not in sample]
    pick = lambda table, names: [table[s] for s in names]
    cal = stats(pick(journal, outside))
    wrong = commented[JOURNAL_ROW_KEPT_AS_COMMENT]
    out = dict(
        paper='Reuter et al. 2020 (ApJ 902, 78), Table F.1',
        tables='https://github.com/spt-smg/publicdata, directory reuter2020_final_dndz',
        files={Path(a).name: hashlib.sha256(Path(a).read_bytes()).hexdigest() for a in args},
        definitions='pull = (z_spec - z_phot) / published error; dz_over_1pz = (z_phot - z_spec) / (1 + z_spec); '
                    'z_spec of SPT0452-50 is the corrected 5.016 throughout',
        sample=sample,
        published_2020=dict(all=stats(list(journal.values())), outside_sample=cal,
                            sample=stats(pick(journal, sample))),
        revised_2024=dict(all=stats(list(revised.values())), outside_sample=stats(pick(revised, outside)),
                          sample=stats(pick(revised, sample)),
                          rows_changed=sum(journal[s][1:] != revised[s][1:] for s in journal)),
        recalibration=dict(from_sources='the %d sources of the table that are not in data/truth.csv' % len(outside),
                           shift_in_published_errors=cal['pull_mean'], width_factor=cal['pull_sd'],
                           meaning='recalibrated centre = z_phot + shift * error; recalibrated width = factor * error'),
        spt0452_50=dict(z_spec_in_journal_table=wrong[0], z_phot=wrong[1], z_phot_err=wrong[2],
                        pull_at_journal_z_spec=round((wrong[0] - wrong[1]) / wrong[2], 2),
                        pull_at_corrected_z_spec=round((zs_now - wrong[1]) / wrong[2], 2),
                        larger_pulls_elsewhere_in_table=sum(
                            abs((zs - zp) / e) >= abs((wrong[0] - wrong[1]) / wrong[2])
                            for s, (zs, zp, e) in journal.items() if s != JOURNAL_ROW_KEPT_AS_COMMENT)))
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print(f'wrote {OUT}')
    for k in ('all', 'outside_sample', 'sample'):
        print('published', k, out['published_2020'][k])
    print('revised all', out['revised_2024']['all'])
    print(out['recalibration'], out['spt0452_50'])


if __name__ == '__main__':
    main()
