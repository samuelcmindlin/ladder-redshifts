#!/usr/bin/env python3
"""Build the cut-out manifest for get_cutouts.sh: one row per pipeline cube to cut, centred on its source.

Inputs
  --files    list of archive file names that name the cubes. ADMIT names (*.cube.I.pbcor.admit.native.tgz) are
             turned into the cube names (*.cube.I.pbcor.fits); ARI-L re-processed products ('.ari_l.') are skipped.
  --ari-l    use the ARI-L cubes instead (*.ari_l.<target>_sci.spwN_<freq>MHz.12m.cube.I.pbcor.fits). Projects
             from before Cycle 4 have no pipeline cubes per window; ARI-L re-imaged them.
  --targets  CSV: target,ra_deg,dec_deg,position_from,priority,radius_arcsec,box_arcsec
             (priority 1 = fetch first; radius = cut-out radius; box = half-size used later by the extractor).
  --cont     targets (comma-separated) whose source position is unknown: their continuum images are listed too
             (whole files, ~0.3 MB each), so the source can be located before extraction.

Output: a tab-separated manifest, sorted by priority:
  file  kind(cutout|file)  target  ra_deg  dec_deg  radius_arcsec  priority
"""
import argparse
import csv
import re
import sys

NAME = re.compile(r'^(member\.uid___\w+?)\.(.+?)_sci\.spw(\d+)\.cube\.I\.pbcor\.(?:admit\.native\.tgz|fits)$')
ARI_L = re.compile(r'^(member\.uid___\w+?)\.ari_l\.(.+?)_sci\.spw(\d+)_\d+MHz\.12m\.cube\.I\.pbcor\.fits$')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--files', required=True)
    ap.add_argument('--targets', required=True)
    ap.add_argument('--cont', default='')
    ap.add_argument('--ari-l', action='store_true')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    tg = {r['target']: r for r in csv.DictReader(open(a.targets))}
    rows, mous_of, skipped = [], {}, set()
    for line in open(a.files):
        name = line.strip()
        if not name or ('.ari_l.' in name) != a.ari_l:
            continue
        m = (ARI_L if a.ari_l else NAME).match(name)
        if not m:
            continue
        member, target, spw = m.group(1), m.group(2), m.group(3)
        if target not in tg:
            skipped.add(target)
            continue
        t = tg[target]
        cube = name if a.ari_l else f'{member}.{target}_sci.spw{spw}.cube.I.pbcor.fits'
        rows.append([cube, 'cutout', target, t['ra_deg'], t['dec_deg'], t['radius_arcsec'], t['priority']])
        mous_of.setdefault(target, {}).setdefault(member, []).append(int(spw))
    for target in [x for x in a.cont.split(',') if x]:
        t = tg[target]
        for member, spws in sorted(mous_of.get(target, {}).items()):
            tag = '_'.join(str(s) for s in sorted(spws))
            rows.append([f'{member}.{target}_sci.spw{tag}.cont.I.tt0.pbcor.fits', 'file', target, t['ra_deg'],
                         t['dec_deg'], t['radius_arcsec'], t['priority']])
    if skipped:
        print('not in the targets file, skipped: ' + ', '.join(sorted(skipped)), file=sys.stderr)
    rows.sort(key=lambda r: (int(r[6]), r[2], r[1] != 'file', r[0]))
    with open(a.out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['file', 'kind', 'target', 'ra_deg', 'dec_deg', 'radius_arcsec', 'priority'])
        w.writerows(rows)
    by = {}
    for r in rows:
        by.setdefault(r[6], set()).add(r[2])
    print(f'{len(rows)} rows: ' + '; '.join(f'priority {p}: {len(v)} targets' for p, v in sorted(by.items())))


if __name__ == '__main__':
    main()
