#!/usr/bin/env python3
"""Turn ALMA ADMIT summary archives (*.admit.native.tgz) into engine-ready spectra, without downloading cubes.

The ALMA archive ships an ADMIT summary next to every pipeline cube: a few MB instead of hundreds. Per spectral
window it holds
  - the cube's per-channel noise (CubeStats "sigma", robust over the whole field) and per-channel peak;
  - the source finder's detections on the line-moment map (SFind2D, 5 sigma);
  - spectra at up to two positions: the cube's brightest voxel and the moment-map peak (CubeSpectrum);
  - the line finder's verdict (LineSegment: >= 4 channels above 5 sigma, gaps <= 3).
Those two spectra sit on the target only when the window contains a line bright enough for ADMIT to find.
Elsewhere they land on noise spikes anywhere in the field.

So the spectrum written here is a pseudo-spectrum:
  measured     window where a CubeSpectrum position lies within --match x beam FWHM of the target. That spectrum
               is used, binned to ~--dv km/s, with a linear baseline refitted on line-free channels (the pipeline's
               continuum subtraction can eat into broad lines), the noise scaled to the line-free scatter, and a
               point-source beam correction for the small position offset.
  zero-filled  window with no spectrum on the target and no sign of emission there: flux 0, noise = the cube's
               per-channel sigma x --nondet-scale (default 10). ADMIT misses lines that peak below ~5 sigma per
               native channel; for a 500 km/s line that is an integrated S/N of up to ~20, so x10 turns a missed
               line into a ~2 sigma deficit, and only lines far brighter than ADMIT's limit count as absent.
               (The first default, x2.5, under-weighted the true redshift in three of the 2016.1.00672.S sources
               where ADMIT had missed a line; see docs/development-notes.md.)
  dropped      window where ADMIT saw emission but no spectrum lies on the target (we cannot tell where the line
               is), or where a higher-priority window already covers the channels (overlaps).
Why a spectrum near the target counts even when ADMIT reports no line: ADMIT's line finder needs >= 4 channels
above 5 sigma, and real lines miss that (SPT0403-58's CO(4-3), a published detection, makes the cube's brightest
voxel, 5.2 sigma, right on the source, yet ADMIT reports no line segment and no source). When the cube's brightest
voxel or moment peak falls within half a beam of the target, that is itself the detection: a noise peak lands there
by chance in well under 1% of windows. The position was chosen for its peak, so the peak channel is biased upwards,
but a line spans many channels. --require-detection restores the strict rule (measured only where ADMIT itself
reports emission).
The target position is the S/N-weighted mean of the source-finder detections within --max-offset of the phase
centre, across all windows (sky offsets, since pixel grids differ between tunings); the phase centre otherwise.

Tar members are read into memory only; nothing is extracted or executed. Parsed windows are cached as JSON, so a
large set can be processed over several short runs (each run continues where the last stopped).

Usage:
  python3 scripts/alma/admit_spectra.py <dir|file.tgz> [...] --out data/admit/spectra [--nondet-scale 10] [--dv 50]
Writes <out>/<TARGET>.csv (spw,freq_ghz,flux_mjy,noise_mjy,continuum_mjy,admit) and <out>/<TARGET>_admit.json.
Then score them with  python -m ladder_redshifts.batch_real --admit
"""
import argparse
import ast
import json
import math
import re
import sys
import tarfile
import time
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import numpy as np

C_KMS = 299792.458
NAME_RE = re.compile(r'member\.uid___(\w+?)_(\w+?)_(\w+?)\.(.+?)_sci\.spw(\d+)\.cube\.I\.pbcor\.admit\.'
                     r'(native|binned16)\.tgz$')
NEED = ('admit.log', 'x.cst.bdp')


# ---------------------------------------------------------------------------------------------------- parsing
def parse_log(log):
    out = {}
    m = re.search(r'RA\s+Axis 1:\s+(\S+)\s+(\S+)\s+(\S+)', log)
    n = re.search(r'DEC\s+Axis 2:\s+(\S+)\s+(\S+)\s+(\S+)', log)
    if m and n:
        out['ra0_deg'], out['dx_arcsec'], out['x0'] = map(float, m.groups())
        out['dec0_deg'], out['dy_arcsec'], out['y0'] = map(float, n.groups())
    b = re.search(r"COMMONBEAM\[\d+\] \{'major': \{'unit': 'arcsec', 'value': ([\d.eE+-]+)\}, "
                  r"'minor': \{'unit': 'arcsec', 'value': ([\d.eE+-]+)\}", log)
    if b:
        out['bmaj'], out['bmin'] = map(float, b.groups())
    else:
        b = re.search(r'Restoring Beam: Major axis:\s+([\d.]+) arcsec , Minor axis:\s+([\d.]+) arcsec', log)
        if b:
            out['bmaj'], out['bmin'] = map(float, b.groups())
    pos = []
    c = re.search(r'CubeStats::maxpos,val=\[\s*(\d+),\s*(\d+)', log)
    if c:
        pos.append([int(c.group(1)), int(c.group(2))])
    s = re.search(r'CubeSum::maxpos,val=\[\s*(\d+),\s*(\d+)\]', log)
    if s and [int(s.group(1)), int(s.group(2))] not in pos:
        pos.append([int(s.group(1)), int(s.group(2))])
    out['spec_pos'] = pos
    w = re.search(r'Writing (\d+) testCubeSpectrum tables', log)
    out['n_spec_tables'] = int(w.group(1)) if w else None
    src = []
    for line in log.splitlines():
        if 'SFind2D_AT.py :' not in line:
            continue
        t = line.split('SFind2D_AT.py :')[1].split()
        if len(t) == 10 and re.match(r'\d{2}:\d{2}:\d{2}', t[0]):
            try:
                src.append(dict(x=float(t[2]), y=float(t[3]), peak=float(t[4]), flux=float(t[5]),
                                snr=float(t[9])))
            except ValueError:
                pass
    out['sfind'] = src
    seg = re.findall(r'Segment Coverage (\d+) / (\d+)', log)
    out['seg_channels'] = sum(int(a) for a, _ in seg) if seg else None
    return out


def parse_cst(text):
    root = ET.fromstring(text)
    cols = ast.literal_eval(root.find('table/columns').text)
    data = np.array(ast.literal_eval(root.find('table/data').text), dtype=float)
    return {c: data[:, i] for i, c in enumerate(cols)}


def parse_tgz(path):
    m = NAME_RE.search(path.name)
    if not m:
        return None
    rec = dict(file=path.name, target=m.group(4), mous=m.group(3), spw=int(m.group(5)), kind=m.group(6))
    got, tabs = {}, {}
    with tarfile.open(path, 'r:gz') as t:
        for mem in t:
            if not mem.isfile():
                continue
            base = mem.name.rsplit('/', 1)[-1]
            if base in NEED:
                got[base] = t.extractfile(mem).read().decode('utf-8', 'replace')
            elif re.fullmatch(r'testCubeSpectrum_(\d+)\.tab', base):
                tabs[int(base.split('_')[1].split('.')[0])] = t.extractfile(mem).read().decode()
    if 'admit.log' not in got or 'x.cst.bdp' not in got:
        rec['error'] = 'missing admit.log or x.cst.bdp'
        return rec
    rec.update(parse_log(got['admit.log']))
    cst = parse_cst(got['x.cst.bdp'])
    rec['freq_ghz'] = cst['frequency'].tolist()
    rec['sigma_jyb'] = cst['sigma'].tolist()
    rec['peak_jyb'] = cst['max'].tolist()
    spectra = []
    for k in sorted(tabs):
        a = np.loadtxt(tabs[k].splitlines(), comments='#', ndmin=2)
        spectra.append(dict(index=k, pos=rec['spec_pos'][k] if k < len(rec['spec_pos']) else None,
                            freq_ghz=a[:, 0].tolist(), flux_jyb=a[:, 1].tolist()))
    rec['spectra'] = spectra
    return rec


# ---------------------------------------------------------------------------------------------------- geometry
def offset_arcsec(rec, x, y):
    """Sky offset of pixel (x, y) from this window's phase centre, in arcsec (east, north)."""
    return (x - rec['x0']) * rec['dx_arcsec'], (y - rec['y0']) * rec['dy_arcsec']


def target_position(recs, max_offset):
    pts, wts = [], []
    for r in recs:
        for s in r.get('sfind', []):
            dx, dy = offset_arcsec(r, s['x'], s['y'])
            if math.hypot(dx, dy) <= max_offset:
                pts.append((dx, dy)); wts.append(s['snr'])
    if not pts:
        return (0.0, 0.0), 'phase centre (no source-finder detection within %.0f")' % max_offset
    p = np.average(np.array(pts), axis=0, weights=np.array(wts))
    return (float(p[0]), float(p[1])), 'source finder: %d detection(s), S/N-weighted' % len(pts)


# ---------------------------------------------------------------------------------------------------- spectra
def binned(freq, flux, sigma, k):
    n = len(freq) // k * k
    if k <= 1:
        return freq, flux, sigma
    f = freq[:n].reshape(-1, k).mean(1)
    s = flux[:n].reshape(-1, k).mean(1)
    e = np.sqrt((sigma[:n].reshape(-1, k) ** 2).sum(1)) / k
    return f, s, e


def refit_baseline(f, s, e, nbox_kms=300.0):
    """Linear baseline on line-free channels; returns (residual, baseline, noise scale, line mask)."""
    dv = abs(np.median(np.diff(f))) / np.median(f) * C_KMS
    nb = max(3, int(round(nbox_kms / dv)))
    line = np.zeros(len(f), bool)
    base = np.zeros(len(f))
    x = (f - f.mean()) / max(np.ptp(f), 1e-9)
    for _ in range(4):
        r = s - base
        sm = np.convolve(r / e, np.ones(nb) / nb, mode='same') * math.sqrt(nb)
        new = np.abs(sm) > 3.0
        new = np.convolve(new, np.ones(2 * (nb // 2) + 1), mode='same') > 0      # grow by ~half a box
        if new.sum() > 0.6 * len(f):                                             # never mask most of a window
            new = np.abs(sm) > 5.0
        line = new
        good = ~line
        if good.sum() < 10:
            break
        c = np.polyfit(x[good], s[good], 1, w=1 / e[good])
        base = np.polyval(c, x)
    r = s - base
    good = ~line
    mad = 1.4826 * np.median(np.abs(r[good] / e[good] - np.median(r[good] / e[good]))) if good.sum() > 10 else 1.0
    return r, base, float(np.clip(mad, 0.7, 3.0)), line


def build_target(recs, args):
    (tx, ty), how = target_position(recs, args.max_offset)
    wins = []
    for r in sorted(recs, key=lambda r: min(r['freq_ghz'])):
        f = np.array(r['freq_ghz']); sig = np.array(r['sigma_jyb'])
        o = np.argsort(f); f = f[o]; sig = sig[o]
        dvn = abs(np.median(np.diff(f))) / np.median(f) * C_KMS
        k = max(1, int(round(args.dv / dvn)))
        bmaj, bmin = r.get('bmaj', 4.0), r.get('bmin', 3.0)
        best = None
        for sp in r['spectra']:
            if sp['pos'] is None:
                continue
            dx, dy = offset_arcsec(r, *sp['pos'])
            d = math.hypot(dx - tx, dy - ty)
            if d <= args.match * bmaj and (best is None or d < best[0]):
                best = (d, sp)
        near_src = any(math.hypot(*(np.subtract(offset_arcsec(r, s['x'], s['y']), (tx, ty)))) <= bmaj
                       for s in r.get('sfind', []))
        info = dict(window=f"{r['target']}_{r['mous']}_spw{r['spw']}", file=r['file'],
                    freq_ghz=[round(float(f.min()), 4), round(float(f.max()), 4)],
                    native_kms=round(dvn, 1), bin=k, sigma_mjy=round(1e3 * float(np.median(sig)), 3),
                    seg_channels=r.get('seg_channels'), sfind_near_target=near_src)
        emission = (r.get('seg_channels') or 0) > 0 or near_src
        if best is not None and (emission or not args.require_detection):
            d, sp = best
            fs = np.array(sp['freq_ghz']); ys = np.array(sp['flux_jyb'])
            oo = np.argsort(fs); ys = np.interp(f, fs[oo], ys[oo])
            att = math.exp(-4 * math.log(2) * d ** 2 / (bmaj * bmin))
            fb, sb, eb = binned(f, ys / att, sig / att, k)
            res, base, scale, line = refit_baseline(fb, sb, eb)
            eb = eb * scale
            info.update(status='measured', admit_detection=bool(emission), offset_arcsec=round(d, 2),
                        beam_factor=round(att, 3),
                        spectrum='CubeSpectrum_%d' % sp['index'], noise_scale=round(scale, 2),
                        line_channels=int(line.sum()), peak_snr=round(float((res / eb).max()), 1))
            wins.append(dict(info=info, f=fb, S=1e3 * res, e=1e3 * eb, kind='measured'))
        elif emission:
            info.update(status='dropped', reason='ADMIT found emission but no spectrum on the target')
            wins.append(dict(info=info, f=None, kind='dropped'))
        else:
            fb, _, eb = binned(f, np.zeros_like(f), sig, k)
            info.update(status='zero-filled', nondet_scale=args.nondet_scale)
            wins.append(dict(info=info, f=fb, S=np.zeros_like(fb), e=1e3 * eb * args.nondet_scale,
                             kind='zero-filled'))
    # overlaps: measured windows first, then quieter zero-filled windows; later ones lose overlapping channels
    order = sorted([w for w in wins if w['f'] is not None],
                   key=lambda w: (w['kind'] != 'measured', float(np.median(w['e']))))
    taken = []
    for w in order:
        keep = np.ones(len(w['f']), bool)
        for lo, hi in taken:
            keep &= ~((w['f'] >= lo) & (w['f'] <= hi))
        if keep.sum() < len(w['f']):
            w['info']['overlap_channels_dropped'] = int((~keep).sum())
        for key in ('f', 'S', 'e'):
            w[key] = w[key][keep]
        if len(w['f']):
            taken.append((w['f'].min() - 1e-6, w['f'].max() + 1e-6))
    return (tx, ty), how, wins


def write_target(name, pos, how, wins, out, args):
    rows = ['spw,freq_ghz,flux_mjy,noise_mjy,continuum_mjy,admit']
    for w in sorted([w for w in wins if w['f'] is not None and len(w['f'])], key=lambda w: w['f'].min()):
        for f, s, e in zip(w['f'], w['S'], w['e']):
            rows.append(f"{w['info']['window']},{f:.6f},{s:.5f},{e:.5f},0.00000,{w['kind']}")
    (out / f'{name}.csv').write_text('\n'.join(rows) + '\n')
    summ = dict(target=name, target_offset_arcsec=[round(pos[0], 2), round(pos[1], 2)], position_from=how,
                settings=dict(nondet_scale=args.nondet_scale, dv_kms=args.dv, match_beams=args.match,
                              max_offset_arcsec=args.max_offset, require_detection=args.require_detection),
                counts={k: sum(w['kind'] == k for w in wins) for k in ('measured', 'zero-filled', 'dropped')},
                windows=[w['info'] for w in wins])
    (out / f'{name}_admit.json').write_text(json.dumps(summ, indent=1))
    return summ


# ---------------------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('inputs', nargs='+')
    ap.add_argument('--out', default='spectra_admit')
    ap.add_argument('--cache', default=None, help='parsed-window cache (default <out>/.admit_cache)')
    ap.add_argument('--nondet-scale', type=float, default=10.0)
    ap.add_argument('--dv', type=float, default=50.0, help='target channel width, km/s')
    ap.add_argument('--match', type=float, default=0.5, help='max spectrum-to-target offset, in beam FWHM')
    ap.add_argument('--max-offset', type=float, default=10.0, help='max target offset from phase centre, arcsec')
    ap.add_argument('--targets', nargs='*', default=None)
    ap.add_argument('--require-detection', action='store_true',
                    help='use a spectrum only where ADMIT itself reports emission (stricter; misses faint lines)')
    ap.add_argument('--budget', type=float, default=1e9, help='stop parsing after this many seconds')
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    cache = Path(a.cache) if a.cache else out / '.admit_cache'; cache.mkdir(parents=True, exist_ok=True)
    files = []
    for p in map(Path, a.inputs):
        files += sorted(p.rglob('*.admit.*.tgz')) if p.is_dir() else [p]
    keyed = {}
    for f in files:
        m = NAME_RE.search(f.name)
        if not m or (a.targets and m.group(4) not in a.targets):
            continue
        key = (m.group(4), m.group(3), int(m.group(5)))
        if m.group(6) == 'native' or key not in keyed:          # prefer native over binned16
            keyed[key] = f
    t0, todo = time.time(), 0
    recs = defaultdict(list)
    failed = defaultdict(list)
    for key, f in sorted(keyed.items()):
        cf = cache / (f.name + '.json')
        if cf.exists():
            rec = json.loads(cf.read_text())
        elif time.time() - t0 < a.budget:
            try:
                rec = parse_tgz(f)
            except (tarfile.TarError, OSError, EOFError, ValueError, SyntaxError) as exc:
                rec = dict(file=f.name, error=f'{type(exc).__name__}: {exc} (incomplete download?)')
            if rec and 'error' not in rec:          # failures are not cached, so they are retried next run
                cf.write_text(json.dumps(rec))
        else:
            todo += 1
            continue
        if rec and 'error' not in rec:
            recs[key[0]].append(rec)
        elif rec:
            failed[key[0]].append(f.name)
            print(f"  unreadable {f.name}: {rec['error']}")
    if todo:
        print(f'time budget reached: {todo} archive(s) still to parse; run the same command again to continue')
        return 2
    for name in sorted(recs):
        if name in failed:
            print(f"{name:14s} NOT WRITTEN: {len(failed[name])} archive(s) unreadable")
            continue
        pos, how, wins = build_target(recs[name], a)
        s = write_target(name, pos, how, wins, out, a)
        meas = [w for w in s['windows'] if w['status'] == 'measured']
        nm = len({r['mous'] for r in recs[name]})
        print(f"{name:14s} {len(recs[name])} windows from {nm} tunings: {s['counts']}  "
              f"target at ({pos[0]:+.1f}\", {pos[1]:+.1f}\") [{how}]")
        for w in meas:
            print(f"    measured {w['window']:28s} {w['freq_ghz'][0]:.3f}-{w['freq_ghz'][1]:.3f} GHz  "
                  f"offset {w['offset_arcsec']:.2f}\"  peak S/N {w['peak_snr']}  noise x{w['noise_scale']}"
                  f"{'' if w['admit_detection'] else '  (no ADMIT detection)'}")
        for w in s['windows']:
            if w['status'] == 'dropped':
                print(f"    DROPPED  {w['window']:28s} ({w['reason']})")
    print(f'wrote {len([n for n in recs if n not in failed])} target(s) to {out}/')
    if failed:
        print(f'INCOMPLETE: {sum(map(len, failed.values()))} archive(s) unreadable for {len(failed)} target(s) '
              f'(still downloading, or corrupt). Rerun when the download has finished.')
        return 3
    print('ALL DONE')
    return 0


if __name__ == '__main__':
    sys.exit(main())
