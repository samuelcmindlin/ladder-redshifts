#!/usr/bin/env python3
"""Keep only small cut-outs of ALMA science cubes, without ever unpacking a whole product bundle.

Usage:  python3 extract_cube_spectra.py INPUT [INPUT ...] --out OUTROOT [--box 12] [--dv 50] [--scratch DIR]
                                       [--targets targets.csv] [--tag spw|freq]
        (old form, still accepted:  python3 extract_cube_spectra.py PRODUCT.tgz OUTDIR)

INPUT can be a pipeline product bundle (.tgz / .tar.gz, streamed; .tar, streamed), a single cube
(.fits, e.g. from the archive's cut-out service), or a directory searched recursively for cubes.
Results go to OUTROOT/<target>/<target>_spwNN.npz (one folder per science target).

For every member matching '*_sci.spw*.cube*.I.pbcor.fits' (one cube per spectral window):
  * copy it to a scratch file (one cube at a time, deleted afterwards - nothing big is kept),
  * cut a box of +/- BOX arcsec around the image reference pixel (the phase centre), or around the target's
    position when --targets gives one (CSV with columns target,ra_deg,dec_deg and optionally box_arcsec;
    sources can sit 10-20" off the phase centre),
  * average channels to ~DV km/s,
  * measure the per-channel noise (robust MAD) in an annulus 1.3*BOX .. 2.6*BOX arcsec. In primary-beam-corrected
    cubes the noise rises away from the pointing centre, so this is higher than the noise at the source (by
    about 1.3 for BOX = 8", more for larger boxes and at higher frequency). The shape of the noise along the
    spectrum is right; ladder_redshifts.realspec.calibrate_noise() sets its scale from the line-free channels
    of the extracted spectrum,
  * save OUTDIR/<target>_spwNN.npz  (freq_hz, cutout [chan,y,x] in Jy/beam, rms_jyb, beam, pixel scale, header).
    With --tag freq the file is <target>_fNNN.NNN.npz (centre frequency in GHz) instead: needed when several
    tunings reuse the same spw numbers, as in multi-tuning spectral scans.
The small continuum (mfs) images of the science target are copied as they are.
Needs numpy and astropy.
"""
import argparse, os, re, shutil, sys, tarfile, tempfile, time, json
import numpy as np
from astropy.io import fits

C_KMS = 299792.458


def log(msg, fh):
    line = time.strftime('%H:%M:%S ') + msg
    print(line, flush=True); fh.write(line + '\n'); fh.flush()


def spectral_axis(h):
    for i in range(1, h['NAXIS'] + 1):
        if str(h.get(f'CTYPE{i}', '')).upper().startswith('FREQ'):
            return i
    raise ValueError('no FREQ axis')


def centre_freq_ghz(path):
    h = fits.getheader(path)
    ia = spectral_axis(h)
    n = h[f'NAXIS{ia}']
    return (h[f'CRVAL{ia}'] + ((n + 1) / 2 - h[f'CRPIX{ia}']) * h[f'CDELT{ia}']) / 1e9


def process_cube(path, box_as, dv_kms, out_npz, fh, centre=None):
    with fits.open(path, memmap=True) as hd:
        h = hd[0].header
        data = hd[0].data
        ia = spectral_axis(h)                         # FITS axis number of frequency
        nchan = h[f'NAXIS{ia}']
        f = h[f'CRVAL{ia}'] + (np.arange(nchan) + 1 - h[f'CRPIX{ia}']) * h[f'CDELT{ia}']
        pix_as = abs(h['CDELT2']) * 3600.0
        cx, cy = int(round(h['CRPIX1'] - 1)), int(round(h['CRPIX2'] - 1))
        if centre is not None:                        # target position (deg) -> pixel, via the celestial WCS
            from astropy.wcs import WCS
            px, py = WCS(h).celestial.world_to_pixel_values(centre[0], centre[1])
            cx, cy = int(round(float(px))), int(round(float(py)))
        nx, ny = h['NAXIS1'], h['NAXIS2']
        nb = int(round(box_as / pix_as)); na = int(round(2.6 * box_as / pix_as))
        # numpy axis order is reversed FITS order: (..., chan, y, x); drop degenerate Stokes axis
        arr = data
        while arr.ndim > 3:
            arr = arr[0]
        if ia != 3 and data.ndim == 4:                # e.g. RA, DEC, STOKES, FREQ
            arr = data[:, 0]
        dv_nat = C_KMS * abs(h[f'CDELT{ia}']) / np.median(f)
        k = max(1, int(round(dv_kms / dv_nat)))
        nout = nchan // k
        y0, y1 = max(cy - na, 0), min(cy + na + 1, ny)
        x0, x1 = max(cx - na, 0), min(cx + na + 1, nx)
        yy, xx = np.mgrid[y0:y1, x0:x1]
        rr = np.hypot((yy - cy) * pix_as, (xx - cx) * pix_as)
        ann = (rr > 1.3 * box_as) & (rr < 2.6 * box_as)
        if not (nb <= cx < nx - nb and nb <= cy < ny - nb):
            raise ValueError(f'box of +/-{box_as}" around pixel ({cx}, {cy}) does not fit in the {nx}x{ny} image')
        cut = np.empty((nout, 2 * nb + 1, 2 * nb + 1), np.float32)
        rms = np.empty(nout, np.float32)
        fo = np.empty(nout)
        for j in range(nout):
            blk = np.asarray(arr[j * k:(j + 1) * k, y0:y1, x0:x1], dtype=np.float64)
            m = np.nanmean(blk, axis=0)
            cut[j] = m[cy - nb - y0:cy + nb + 1 - y0, cx - nb - x0:cx + nb + 1 - x0]
            v = m[ann]; v = v[np.isfinite(v)]
            rms[j] = 1.4826 * np.median(np.abs(v - np.median(v))) if v.size else np.nan
            fo[j] = f[j * k:(j + 1) * k].mean()
        beams = None
        for ext in hd[1:]:
            if ext.name == 'BEAMS':
                b = ext.data
                beams = np.array([b['BMAJ'], b['BMIN'], b['BPA']], float).T
                units = [ext.header.get('TUNIT1', ''), ext.header.get('TUNIT2', ''), ext.header.get('TUNIT3', '')]
        if beams is None:
            beams = np.array([[h.get('BMAJ', np.nan) * 3600, h.get('BMIN', np.nan) * 3600, h.get('BPA', np.nan)]])
            units = ['arcsec', 'arcsec', 'deg']
        hdr = {kk: str(h.get(kk)) for kk in ('OBJECT', 'BUNIT', 'SPECSYS', 'RESTFRQ', 'CRVAL1', 'CRVAL2', 'CDELT1',
                                              'CDELT2', 'NAXIS1', 'NAXIS2', f'NAXIS{ia}', f'CDELT{ia}', 'DATE-OBS')}
        hdr['centre_pix'] = f'{cx},{cy}'
        hdr['centre_from'] = 'target position' if centre is not None else 'reference pixel'
        hdr['annulus_pixels'] = str(int(ann.sum()))
        np.savez_compressed(out_npz, freq_hz=fo, cutout=cut, rms_jyb=rms, beams=beams, beam_units=np.array(units),
                            pix_arcsec=pix_as, chan_bin=k, dv_native_kms=dv_nat, header=json.dumps(hdr))
        log(f'  saved {os.path.basename(out_npz)}: {nout} chans of {dv_nat * k:.1f} km/s, '
            f'{fo.min() / 1e9:.3f}-{fo.max() / 1e9:.3f} GHz, box {2 * nb + 1} px ({pix_as:.3f}"/px), '
            f'median rms {np.nanmedian(rms) * 1e3:.3f} mJy/beam', fh)


def target_and_tag(name, path=None):
    """Target name and spectral-window tag from a pipeline file name, or from the FITS header."""
    m = re.search(r'\.([^.]+)_sci\.spw(\d+)\.', name)
    if m:
        return m.group(1), f'spw{m.group(2)}'
    if path is not None:
        h = fits.getheader(path)
        obj = str(h.get('OBJECT', 'target')).strip().replace(' ', '_') or 'target'
        ia = spectral_axis(h)
        f0 = h[f'CRVAL{ia}'] + (1 - h[f'CRPIX{ia}']) * h[f'CDELT{ia}']
        return obj, f'f{f0 / 1e9:.3f}'
    return 'target', 'cube'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('inputs', nargs='+')
    ap.add_argument('--out', default=None, help='output root; results go to OUT/<target>/')
    ap.add_argument('--box', type=float, default=12.0, help='half-size of the cut-out, arcsec')
    ap.add_argument('--dv', type=float, default=50.0, help='target channel width, km/s')
    ap.add_argument('--scratch', default=os.path.join(os.path.expanduser('~'), 'scratch'))
    ap.add_argument('--targets', default=None,
                    help='CSV target,ra_deg,dec_deg[,box_arcsec]: centre each box on its target (optional box size)')
    ap.add_argument('--tag', choices=['spw', 'freq'], default='spw',
                    help='name outputs by spw number (default) or by centre frequency (multi-tuning scans)')
    ap.add_argument('--allow-failures', action='store_true',
                    help='exit with status 0 even if some cubes could not be processed (a deliberately partial run)')
    a = ap.parse_args()
    failed, done = [], []
    centres, boxes = {}, {}
    if a.targets:
        import csv
        for r in csv.DictReader(open(a.targets)):
            centres[r['target']] = (float(r['ra_deg']), float(r['dec_deg']))
            if r.get('box_arcsec'):
                boxes[r['target']] = float(r['box_arcsec'])
    flat_out = None
    if a.out is None:                      # old form: last positional argument is the output folder
        if len(a.inputs) < 2:
            ap.error('give --out OUTROOT')
        flat_out = a.inputs[-1]; a.inputs = a.inputs[:-1]; a.out = flat_out
    os.makedirs(a.out, exist_ok=True); os.makedirs(a.scratch, exist_ok=True)
    fh = open(os.path.join(a.out, 'extract_log.txt'), 'a')
    pat_cube = re.compile(r'\.cube.*\.I\.pbcor\.fits$|_sci\.spw\d+\.cube.*\.fits$')
    pat_mfs = re.compile(r'_sci\.spw\d+\.mfs.*\.I\.pbcor\.fits$')

    def outpath(target, tag):
        d = a.out if flat_out else os.path.join(a.out, target)
        os.makedirs(d, exist_ok=True)
        return os.path.join(d, f'{target}_{tag}.npz'), d

    def handle_fits(path, name):
        target, tag = target_and_tag(name, path)
        if a.tag == 'freq':
            try:
                tag = f'f{centre_freq_ghz(path):.3f}'
            except Exception as e:
                log(f'  FAILED {name}: {e!r}', fh); failed.append(name); return
        out, _ = outpath(target, tag)
        if os.path.exists(out):
            log(f'skip (done) {name}', fh); return
        if a.targets and target not in centres:
            log(f'  note: {target} not in {a.targets}; using the reference pixel', fh)
        log(f'cube {name}', fh)
        try:
            process_cube(path, boxes.get(target, a.box), a.dv, out, fh, centre=centres.get(target))
            done.append(name)
        except Exception as e:
            log(f'  FAILED {name}: {e!r}', fh); failed.append(name)

    for inp in a.inputs:
        log(f'start {inp}', fh)
        if os.path.isdir(inp):
            for root, _, files in os.walk(inp):
                for fn in sorted(files):
                    if fn.endswith('.fits') and ('pbcor' in fn or 'cube' in fn or 'cutout' in fn.lower()):
                        if not fn.endswith('.pb.fits') and 'mask' not in fn and '.mfs.' not in fn \
                                and '.cont.' not in fn:
                            handle_fits(os.path.join(root, fn), fn)
        elif inp.endswith('.fits'):
            handle_fits(inp, os.path.basename(inp))
        else:
            mode = 'r|gz' if inp.endswith(('.tgz', '.gz')) else 'r|'
            with tarfile.open(inp, mode) as tf:
                for m in tf:
                    name = os.path.basename(m.name)
                    if pat_cube.search(name) and '_sci.' in name and not name.endswith('.pb.fits'):
                        target, tag = target_and_tag(name)
                        out, _ = outpath(target, tag)
                        if a.tag == 'spw' and os.path.exists(out):
                            log(f'skip (done) {name}', fh); continue
                        log(f'cube {name} ({m.size / 1e9:.2f} GB)', fh)
                        tmp = os.path.join(a.scratch, 'cube.fits')
                        with tf.extractfile(m) as src, open(tmp, 'wb') as dst:
                            shutil.copyfileobj(src, dst, 16 * 1024 * 1024)
                        try:
                            if a.tag == 'freq':
                                out, _ = outpath(target, f'f{centre_freq_ghz(tmp):.3f}')
                                if os.path.exists(out):
                                    log(f'skip (done) {name}', fh); continue
                            process_cube(tmp, boxes.get(target, a.box), a.dv, out, fh, centre=centres.get(target))
                            done.append(name)
                        except Exception as e:
                            log(f'  FAILED {name}: {e!r}', fh); failed.append(name)
                        finally:
                            os.remove(tmp)
                    elif pat_mfs.search(name):
                        target, _ = target_and_tag(name)
                        _, d = outpath(target, 'mfs')
                        with tf.extractfile(m) as src, open(os.path.join(d, name), 'wb') as dst:
                            shutil.copyfileobj(src, dst)
                        log(f'mfs  {name}', fh)
    if failed:
        log(f'INCOMPLETE: {len(done)} cubes extracted, {len(failed)} FAILED (listed above)', fh)
        if not a.allow_failures:
            sys.exit(1)
    else:
        log(f'DONE: {len(done)} cubes extracted', fh)


if __name__ == '__main__':
    main()
