"""Run the ladder-redshift engine on a REAL spectrum: ALMA Band 2 science verification of SPT0027-50.

Input: the small cut-outs written by scripts/alma/extract_cube_spectra.py (one .npz per spectral window,
from the ALMA pipeline product cubes of project 2011.0.00025.SV), staged into data/SPT0027-50/.

Steps:
  1. aperture: collapse all channels into a continuum(+line) image, keep the connected >3 sigma region
     at the centre, grown by one beam; also report a fixed 3" circle for comparison;
  2. spectrum: sum Jy/beam over the aperture / beam area -> mJy; per-channel noise from the off-source
     annulus rms x sqrt(N_pix / beam_pix); per-window continuum removed with a sigma-clipped linear fit;
     the noise is then put on the scale of the spectrum's own line-free channels (calibrate_noise);
  3. engine: put S/N on the 50 km/s log-frequency grid, with per-channel noise as the template weights;
     run the SLED likelihood (robust prior), flat sum, Sudoku-style rule, naive spacing; uniform and
     SPT-like redshift priors;
  4. stress tests on the real data: (a) add noise to lower the S/N; (b) keep only parts of the coverage
     (a Band 3-like subset, single tunings); (c) drop lines, to see when it becomes ambiguous.
Writes data/SPT0027-50_band2_spectrum.csv, results/redshift_real_SPT0027.json and
figures/redshift_real_SPT0027.png.
ALMA data: ADS/JAO.ALMA#2011.0.00025.SV (public science verification data).
"""
from pathlib import Path as _Path
from .paths import AUDIO, FIG, RESULTS, DATA  # noqa: E402
import glob
import json
import numpy as np
from scipy import ndimage
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
from .ladder import (Grid, Engine, line_list, templates_grid, matched_filter, detect, naive_spacing, BANK,  # noqa: E402
                    sudoku_rule, classify, prior_spt, C_KMS, DLN, PIX_KMS, TOL_LN, CONF_WIN, CO, CI10, CI21,
                    detection_index)

SRC = 'SPT0027-50'
Z_PUBLISHED = 3.444          # Band 3 + Band 2 literature value (De Breuck 2026 talk; Bakx 2026)
INK, INK2, GRID_C, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
COL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300']
# lines NOT in the engine's template (dense-gas tracers etc.), drawn for reference only
EXTRA = {'HCN(4-3)': 354.505, 'HCO+(4-3)': 356.734, 'HNC(4-3)': 362.630, 'H2O(211-202)': 752.033}


def load(src_dir):
    spws = []
    for f in sorted(glob.glob(str(src_dir / f'{SRC}_spw*.npz'))):
        a = np.load(f, allow_pickle=False)
        spws.append(dict(name=_Path(f).stem, freq=a['freq_hz'] / 1e9, cut=a['cutout'].astype(np.float64),
                         rms=a['rms_jyb'].astype(np.float64), beams=a['beams'], pix=float(a['pix_arcsec']),
                         k=int(a['chan_bin'])))
    spws.sort(key=lambda s: s['freq'].min())
    return spws


def beam_pix(s):
    b = np.nanmedian(s['beams'], axis=0)                 # (bmaj, bmin, bpa) in arcsec
    return np.pi * b[0] * b[1] / (4 * np.log(2)) / s['pix'] ** 2, b


def remap(mask, ref_pix, shape, pix):
    """Map a mask drawn on one pixel grid onto another (both centred on the phase centre), by angle."""
    if mask.shape == shape and abs(ref_pix - pix) < 1e-6:
        return mask
    ny, nx = shape
    yy, xx = np.mgrid[:ny, :nx]
    ry = np.rint((yy - ny // 2) * pix / ref_pix).astype(int) + mask.shape[0] // 2
    rx = np.rint((xx - nx // 2) * pix / ref_pix).astype(int) + mask.shape[1] // 2
    ok = (ry >= 0) & (ry < mask.shape[0]) & (rx >= 0) & (rx < mask.shape[1])
    out = np.zeros(shape, bool)
    out[ok] = mask[ry[ok], rx[ok]]
    return out


def aperture(spws):
    """Continuum(+line) >3 sigma region touching the central 3", grown by half a beam. Built on the
    pixel grid shared by most windows; other grids get it by remap(). Falls back to a one-beam circle
    when nothing is detected (faint continuum)."""
    from collections import Counter
    key = Counter((s['cut'].shape[1:], round(s['pix'], 6)) for s in spws).most_common(1)[0][0]
    ref = [s for s in spws if (s['cut'].shape[1:], round(s['pix'], 6)) == key]
    num = sum(np.nansum(s['cut'] / s['rms'][:, None, None] ** 2, axis=0) for s in ref)
    den = sum(np.nansum(1 / s['rms'] ** 2) for s in ref)
    img = num / den                                       # inverse-variance weighted mean over all channels
    ny, nx = img.shape
    pix = ref[0]['pix']
    yy, xx = np.mgrid[:ny, :nx]
    rr = np.hypot(yy - ny // 2, xx - nx // 2) * pix
    off = img[rr > 8]
    sig = 1.4826 * np.median(np.abs(off - np.median(off)))
    lab, n = ndimage.label(img > 3 * sig)
    keep = np.zeros_like(img, bool)
    for i in range(1, n + 1):
        reg = lab == i
        if rr[reg].min() < 3.0 and reg.sum() > 5:         # connected regions touching the central 3"
            keep |= reg
    bp, b = beam_pix(ref[0])
    if not keep.any():
        keep = rr < max(b[0], 1.0)                        # fallback: one-beam circle at the phase centre
    grow = max(1, int(round(b[0] / pix / 2)))
    mask = ndimage.binary_dilation(keep, iterations=grow)
    circ = rr < 3.0
    return img, sig, mask, circ


def spectrum(spws, mask):
    out = []
    from collections import Counter
    key = Counter((s['cut'].shape[1:], round(s['pix'], 6)) for s in spws).most_common(1)[0][0]
    ref_pix = key[1]
    for s in spws:
        bp, _ = beam_pix(s)
        m = remap(mask, ref_pix, s['cut'].shape[1:], s['pix'])
        npix = m.sum()
        S = np.nansum(s['cut'][:, m], axis=1) / bp * 1e3                # mJy
        sig = s['rms'] * np.sqrt(npix / bp) * 1e3                       # mJy
        # continuum: sigma-clipped linear fit per spectral window
        x = s['freq'] - s['freq'].mean()
        good = np.isfinite(S)
        for _ in range(5):
            p = np.polyfit(x[good], S[good], 1, w=1 / sig[good])
            res = (S - np.polyval(p, x)) / sig
            good = np.isfinite(S) & (np.abs(res) < 2.5)
        cont = np.polyval(p, x)
        out.append(dict(spw=s['name'], freq=s['freq'], S=S - cont, cont=cont, sig=sig))
    return calibrate_noise(out)


def spectrum_point(spws, radius_as=0.75):
    """Spectrum of an unresolved source at the centre of each cut-out (the cut-outs must be centred on the
    source): the mean surface brightness within radius_as (well inside the beam), in mJy (Jy/beam = Jy for a
    point source), with the cut-out's per-channel pixel noise, rescaled by calibrate_noise(). Same continuum
    handling as spectrum().
    Used when the beam (~4") is much larger than the source and the cubes are continuum-subtracted, where a
    detection-grown aperture would sum several beams of noise."""
    out = []
    for s in spws:
        ny, nx = s['cut'].shape[1:]
        yy, xx = np.mgrid[:ny, :nx]
        m = np.hypot(yy - ny // 2, xx - nx // 2) * s['pix'] <= radius_as
        S = np.nanmean(s['cut'][:, m], axis=1) * 1e3
        sig = s['rms'] * 1e3
        x = s['freq'] - s['freq'].mean()
        good = np.isfinite(S) & np.isfinite(sig)
        for _ in range(5):
            p = np.polyfit(x[good], S[good], 1, w=1 / sig[good])
            res = (S - np.polyval(p, x)) / sig
            good = np.isfinite(S) & np.isfinite(sig) & (np.abs(res) < 2.5)
        cont = np.polyval(p, x)
        out.append(dict(spw=s['name'], freq=s['freq'], S=S - cont, cont=cont, sig=sig))
    return calibrate_noise(out)


def noise_scale(S, sig, nmin=30):
    """Scatter of flux/noise over the line-free channels of one window: 1 when the stated noise is right, below 1
    when it is too large. Lines are removed by iterative 3-sigma clipping about the median (robust scale), and
    the standard deviation of the remaining channels is corrected for that clipping. NaN for windows with fewer
    than nmin usable channels."""
    with np.errstate(divide='ignore', invalid='ignore'):
        s = np.asarray(S, float) / np.asarray(sig, float)
    ok = np.isfinite(s)
    if ok.sum() < nmin:
        return np.nan
    m = ok.copy()
    for _ in range(5):
        med = np.median(s[m])
        sd = 1.4826 * np.median(np.abs(s[m] - med))
        if not sd > 0:
            return np.nan
        m = ok & (np.abs(s - med) < 3 * sd)
    return float(np.std(s[m]) / 0.9866)                     # 0.9866: std of a unit normal truncated at +/-3


def calibrate_noise(spec):
    """Put each window's per-channel noise on the scale of the spectrum's own line-free channels.

    The extractor measures noise in an annulus around the source in primary-beam-corrected cubes, where the
    noise rises away from the pointing centre, so the stated noise is too high at the source: by about 1.3 for
    the 8" boxes of 2016.1.00672.S, more for larger boxes and at higher frequency (docs/development-notes.md,
    item 20). The channel-to-channel shape of that noise is right; only its scale is off.

    Each window's scale is measured with noise_scale(). One window has about 120 channels, so its own value is
    uncertain by about 6%; the scale varies smoothly with frequency, so with six or more windows spanning at
    least 5% in frequency a straight line in frequency is fitted to the windows' values (one pass of 3-sigma
    outlier rejection) and each window gets the fitted value. With fewer windows they share the median.
    Windows that are not measured spectra (zero-filled ADMIT windows) are left alone.

    Returns new window dicts with 'sig' rescaled, 'noise_scale' (the factor applied) and 'noise_scatter' (the
    window's own measured value)."""
    meas = np.array([w.get('kind', 'measured') == 'measured' for w in spec])
    raw = np.array([noise_scale(w['S'], w['sig']) if m else np.nan for w, m in zip(spec, meas)])
    good = np.isfinite(raw)
    if not good.any():
        return [dict(w) for w in spec]
    fr = np.array([float(np.nanmean(w['freq'])) for w in spec])
    applied = np.full(len(spec), float(np.median(raw[good])))
    if good.sum() >= 6 and np.ptp(fr[good]) > 0.05 * fr[good].mean():
        x = fr - fr[good].mean()
        keep = good.copy()
        for _ in range(2):
            coef = np.polyfit(x[keep], raw[keep], 1)
            res = raw - np.polyval(coef, x)
            sd = 1.4826 * np.median(np.abs(res[keep] - np.median(res[keep])))
            keep = good & (np.abs(res) < 3 * max(sd, 1e-6))
            if keep.sum() < 6:
                break
        fit = np.polyval(coef, x)
        if np.all(fit[meas] > 0.5 * raw[good].min()) and np.all(fit[meas] < 2.0 * raw[good].max()):
            applied = fit
    out = []
    for w, m, a, r in zip(spec, meas, applied, raw):
        if not m:
            out.append(dict(w)); continue
        out.append(dict(w, sig=np.asarray(w['sig'], float) * a, noise_scale=float(a),
                        noise_scatter=float(r) if np.isfinite(r) else None))
    return out


NOISE_TABLE_COLUMNS = ['spectrum', 'spw', 'freq_ghz', 'channels', 'scatter_measured', 'scale_applied']


def record_noise_scales(label, spec, path=None):
    """Record the factors calibrate_noise() applied to `spec` in data/noise_calibration.csv (earlier rows for
    the same `label` are replaced). noise_mjy in the spectrum file = extractor noise x scale_applied."""
    import csv as _csv
    path = (DATA / 'noise_calibration.csv') if path is None else _Path(path)
    rows = []
    if path.exists():
        with open(path, newline='') as fh:
            rows = [r for r in _csv.DictReader(fh) if r['spectrum'] != label]
    for w in spec:
        if w.get('noise_scale') is None:
            continue
        sc = w.get('noise_scatter')
        rows.append(dict(spectrum=label, spw=w['spw'], freq_ghz=f"{float(np.nanmean(w['freq'])):.3f}",
                         channels=str(len(w['freq'])), scatter_measured='' if sc is None else f'{sc:.3f}',
                         scale_applied=f"{w['noise_scale']:.4f}"))
    rows.sort(key=lambda r: (r['spectrum'], float(r['freq_ghz'])))
    with open(path, 'w', newline='') as fh:
        wr = _csv.DictWriter(fh, fieldnames=NOISE_TABLE_COLUMNS)
        wr.writeheader(); wr.writerows(rows)


def grid_windows(spec, grid):
    """Put spectral windows on the engine's log-frequency grid, independently of their order.

    Each window's flux and noise are interpolated onto the grid pixels it spans (never across a gap of more than
    two pixels inside a window). Where windows overlap, the measurements are combined in flux units with
    inverse-variance weights, as independent measurements (different tunings or spectral windows have independent
    noise). Returns (grid, S/N per pixel, noise per pixel (inf where uncovered), coverage mask).

    The engine assumes channels of about one grid pixel (50 km/s). Wider channels make neighbouring pixels
    correlated after interpolation, which raises the matched filter's noise by about sqrt(width / 50 km/s);
    for runs of channels wider than 1.25 pixels the noise is inflated by that factor (with a warning).
    Channels much narrower than a pixel are only sampled, not averaged: bin them to ~50 km/s first."""
    num = np.zeros(grid.npix)
    den = np.zeros(grid.npix)
    key = lambda w: (float(np.nanmin(w['freq'])), float(np.nanmax(w['freq'])), str(w.get('spw', '')))
    for w in sorted(spec, key=key):
        f, S, sg = (np.asarray(w[k], float) for k in ('freq', 'S', 'sig'))
        ok = np.isfinite(f) & np.isfinite(S) & np.isfinite(sg) & (sg > 0)
        if not ok.any():
            continue
        o = np.argsort(f[ok])                               # np.interp needs ascending sample points
        f, S, sg = f[ok][o], S[ok][o], sg[ok][o]
        p = grid.pix(f)
        cuts = np.flatnonzero(np.diff(p) > 2.0) + 1          # gaps inside a window are not bridged
        for seg in np.split(np.arange(len(p)), cuts):
            ps, Ss, ns = p[seg], S[seg], sg[seg]
            if len(seg) > 1:
                width = float(np.median(np.diff(ps))) * PIX_KMS
                if width > 1.25 * PIX_KMS:
                    warnings.warn(f"window {w.get('spw', '?')}: {width:.0f} km/s channels; noise inflated by "
                                  f"{np.sqrt(width / PIX_KMS):.2f} for the 50 km/s grid")
                    ns = ns * np.sqrt(width / PIX_KMS)
                elif width < 0.5 * PIX_KMS:
                    warnings.warn(f"window {w.get('spw', '?')}: {width:.0f} km/s channels are only sampled on the "
                                  f"50 km/s grid; bin them to ~50 km/s for full sensitivity")
            idx = np.clip(np.rint(ps).astype(int), 0, grid.npix - 1)
            g = np.arange(idx.min(), idx.max() + 1)
            s_g = np.interp(g, ps, Ss)
            n_g = np.interp(g, ps, ns)
            num[g] += s_g / n_g ** 2
            den[g] += 1.0 / n_g ** 2
    cov = den > 0
    sig = np.full(grid.npix, np.inf)
    sig[cov] = 1.0 / np.sqrt(den[cov])
    d = np.zeros(grid.npix)
    d[cov] = num[cov] / np.sqrt(den[cov])                  # (combined flux) / (combined noise)
    return grid, d, sig, cov


def to_grid(spec, nu_lo=60.0, nu_hi=125.0):
    """Place the spectrum on the engine's 50 km/s log-frequency grid between nu_lo and nu_hi GHz."""
    return grid_windows(spec, Grid(nu_lo, nu_hi))


def degrade(spec, factor, rng):
    """Raise the noise by `factor` (>=1): add Gaussian noise of sig*sqrt(factor^2-1) and scale sig."""
    out = []
    for w in spec:
        extra = rng.normal(0, 1, len(w['S'])) * w['sig'] * np.sqrt(max(factor ** 2 - 1, 0))
        out.append(dict(w, S=w['S'] + extra, sig=w['sig'] * factor))
    return out


def engine_likelihood(grid, d, sig, cov, bank=BANK):
    """The part of run_engine that does not depend on the redshift prior: engines, likelihoods and detections.
    Several priors can then be applied to one spectrum (engine_posterior) without repeating it."""
    T, lw = templates_grid(robust=True)
    sig_e = np.where(np.isfinite(sig), sig / np.nanmin(sig[np.isfinite(sig)]), 1.0)
    e = Engine(grid, cov, T, sigma=sig_e, logw=lw, zmin=0.5, zmax=8.0)
    # comparison baseline: equal weights on the matched-filter S/N of every covered line (no noise weighting)
    ef = Engine(grid, cov, np.ones(T.shape[1]), zmin=0.5, zmax=8.0)
    Y = matched_filter(d, bank, cov)
    # solutions are chosen with a velocity offset of up to DV_MAX allowed per line; the redshift quoted for a
    # solution is the one at which its lines sit where predicted (the exact-position likelihood, dv_max = 0)
    yd = Y[detection_index(bank)]
    return dict(grid=grid, cov=cov, e=e, ef=ef, L=e.loglike(Y, bank), L0=e.loglike(Y, bank, dv_max=0),
                Lf=ef.loglike(Y, bank), Lf0=ef.loglike(Y, bank, dv_max=0), yd=yd, dets=detect(yd, cov))


def engine_posterior(s, prior=None):
    """The results of run_engine for one redshift prior, from the output of engine_likelihood."""
    grid, cov, e, ef, L0, yd, dets = (s[k] for k in ('grid', 'cov', 'e', 'ef', 'L0', 'yd', 'dets'))
    zs, c, p = e.posterior(s['L'], prior=prior, refine=L0)
    zf, cf, _ = ef.posterior(s['Lf'], prior=prior, refine=s['Lf0'])
    zn = naive_spacing(dets, grid)
    zr, nsol = sudoku_rule(dets, grid, cov, yd)
    # candidate list: posterior mass per cluster
    exact = L0 + e.log_prior(prior=prior)
    order = np.argsort(p, kind='stable')[::-1]
    cands, used = [], np.zeros(len(p), bool)
    for i in order:
        if used[i] or p[i] < 1e-4 * p.max():
            continue
        i = e.locate(i, exact)
        win = np.abs(e.lnz - e.lnz[i]) < CONF_WIN
        cands.append((float(np.exp(e.lnz[i]) - 1), float(p[win].sum())))
        used |= np.abs(e.lnz - e.lnz[i]) < 0.02
        if len(cands) >= 6:
            break
    return dict(sled=(float(zs), float(c)), flat=(float(zf), float(cf)), naive=zn, sudoku=(zr, nsol),
                detections=[(float(np.exp(grid.lnmin + q * DLN)), float(s_)) for q, s_ in dets],
                candidates=cands, lnz=e.lnz, post=p)


def run_engine(grid, d, sig, cov, prior=None, bank=BANK):
    return engine_posterior(engine_likelihood(grid, d, sig, cov, bank), prior)


def read_csv_spectrum(path):
    import csv as _csv
    with open(path, newline='') as fh:
        rows = list(_csv.DictReader(fh))
    out = []
    for name in dict.fromkeys(r['spw'] for r in rows):
        rr = [r for r in rows if r['spw'] == name]
        out.append(dict(spw=name, freq=np.array([float(r['freq_ghz']) for r in rr]),
                        S=np.array([float(r['flux_mjy']) for r in rr]),
                        sig=np.array([float(r['noise_mjy']) for r in rr]),
                        cont=np.array([float(r['continuum_mjy']) for r in rr])))
    return out


def ok(z):
    return z is not None and abs(np.log1p(z) - np.log1p(Z_PUBLISHED)) < TOL_LN


def main():
    src_dir = DATA / SRC
    spws = load(src_dir)
    csv = DATA / f'{SRC}_band2_spectrum.csv'
    imgf = DATA / f'{SRC}_band2_collapsed.npz'
    if spws:                                   # full route: from the cube cut-outs
        img, csig, mask, circ = aperture(spws)
        spec = spectrum(spws, mask)
        spec_c = spectrum(spws, circ)
        rows = ['spw,freq_ghz,flux_mjy,noise_mjy,continuum_mjy']
        for w in spec:
            for f, s, n, c in zip(w['freq'], w['S'], w['sig'], w['cont']):
                rows.append(f"{w['spw']},{f:.6f},{s:.5f},{n:.5f},{c:.5f}")
        csv.write_text('\n'.join(rows) + '\n')
        rows = ['spw,freq_ghz,flux_mjy,noise_mjy,continuum_mjy']
        for w in spec_c:
            for f, s, n, c in zip(w['freq'], w['S'], w['sig'], w['cont']):
                rows.append(f"{w['spw']},{f:.6f},{s:.5f},{n:.5f},{c:.5f}")
        (DATA / f'{SRC}_band2_spectrum_circle3arcsec.csv').write_text('\n'.join(rows) + '\n')
        record_noise_scales(csv.name, spec)                          # factors applied by calibrate_noise
        record_noise_scales(f'{SRC}_band2_spectrum_circle3arcsec.csv', spec_c)
        img = img.astype(np.float32)               # same precision as the saved copy, so both routes plot identically
        np.savez_compressed(imgf, img=img, csig=csig, mask=mask, circ=circ)
        spec = read_csv_spectrum(csv)              # use the saved (rounded) numbers, as the light route does
        spec_c = read_csv_spectrum(DATA / f'{SRC}_band2_spectrum_circle3arcsec.csv')
    elif csv.exists():                         # light route: from the saved 1-D spectra (cut-outs not shipped)
        print('cut-outs not found; using the saved 1-D spectra in data/')
        spec = read_csv_spectrum(csv); spec_c = read_csv_spectrum(DATA / f'{SRC}_band2_spectrum_circle3arcsec.csv')
        z = np.load(imgf); img, csig, mask = z['img'], float(z['csig']), z['mask']
        spws = [None] * len(spec)
    else:
        print(f'no data for {SRC}: run scripts/alma/extract_cube_spectra.py first (see data/README.md)'); return
    grid, d, sig, cov = to_grid(spec)
    res = {}
    for pname, pr in (('uniform', None), ('SPT-like', prior_spt)):
        res[pname] = run_engine(grid, d, sig, cov, prior=pr)
    resc = run_engine(*to_grid(spec_c))
    print(f'{SRC}: {len(spws)} spectral windows, {sum(len(w["freq"]) for w in spec)} channels, '
          f'{min(w["freq"].min() for w in spec):.2f}-{max(w["freq"].max() for w in spec):.2f} GHz')
    print(f'aperture: {mask.sum()} px (continuum >3 sigma, grown by half a beam)')
    for pname, r in res.items():
        print(f'\n[{pname} prior] detections >5 sigma: ' + ', '.join(f'{f:.3f} GHz ({s:.1f})' for f, s in r['detections']))
        print(f"  SLED likelihood: z = {r['sled'][0]:.4f}, confidence {r['sled'][1]:.3f}  -> {'CORRECT' if ok(r['sled'][0]) else 'wrong'}")
        print(f"  flat sum:        z = {r['flat'][0]:.4f}, confidence {r['flat'][1]:.3f}  -> {'CORRECT' if ok(r['flat'][0]) else 'wrong'}")
        print(f"  Sudoku rule:     {r['sudoku'][0]} ({r['sudoku'][1]} surviving solutions)")
        print(f"  naive spacing:   {r['naive']}")
        print('  candidates (z, posterior mass): ' + ', '.join(f'({z:.3f}, {m:.3f})' for z, m in r['candidates']))
    print(f"\n3\" circle aperture: SLED z = {resc['sled'][0]:.4f}, conf {resc['sled'][1]:.3f}")

    # stress tests on the real spectrum
    rng = np.random.default_rng(5)
    factors = [1, 1.5, 2, 3, 4, 6, 8]
    deg = []
    for fct in factors:
        n_ok = n_conf_ok = n_conf_bad = 0; confs = []
        nrep = 1 if fct == 1 else 30
        for _ in range(nrep):
            sp = spec if fct == 1 else degrade(spec, fct, rng)
            r = run_engine(*to_grid(sp), prior=prior_spt)
            z, c = r['sled']
            n_ok += ok(z); n_conf_ok += ok(z) and c >= 0.95; n_conf_bad += (not ok(z)) and c >= 0.95; confs.append(c)
        deg.append(dict(noise_factor=fct, reps=nrep, map_correct=n_ok / nrep, confident_correct=n_conf_ok / nrep,
                        confidently_wrong=n_conf_bad / nrep, median_conf=float(np.median(confs))))
    print('\nStress test 1 (noise raised by a factor, SPT-like prior):')
    for q in deg:
        print(f"  x{q['noise_factor']}: MAP correct {q['map_correct']:.2f}, confident&correct {q['confident_correct']:.2f}, "
              f"confidently wrong {q['confidently_wrong']:.2f}, median confidence {q['median_conf']:.2f}")
    subsets = {'full Band 2 SV coverage': None, 'Band 3-like (>= 84.2 GHz only)': (84.2, 200),
               'below 84 GHz only (CO(3-2) alone)': (0, 84.0), '100-108 GHz only (CO(4-3) alone)': (100.0, 108.0)}
    sub_res = {}
    print('\nStress test 2 (coverage subsets, real data, SPT-like prior):')
    for lab, rngf in subsets.items():
        sp = spec if rngf is None else [dict(w, S=w['S'][(w['freq'] >= rngf[0]) & (w['freq'] <= rngf[1])],
                                             sig=w['sig'][(w['freq'] >= rngf[0]) & (w['freq'] <= rngf[1])],
                                             freq=w['freq'][(w['freq'] >= rngf[0]) & (w['freq'] <= rngf[1])])
                                        for w in spec]
        sp = [w for w in sp if len(w['freq']) > 3]
        if not sp:
            continue
        r = run_engine(*to_grid(sp), prior=prior_spt)
        sub_res[lab] = dict(sled=r['sled'], candidates=r['candidates'][:4], sudoku=r['sudoku'],
                            correct=ok(r['sled'][0]))
        print(f"  {lab:38s} SLED z = {r['sled'][0]:.4f} conf {r['sled'][1]:.3f} ({'correct' if ok(r['sled'][0]) else 'WRONG'}); "
              f"top candidates {[(round(z, 3), round(m, 3)) for z, m in r['candidates'][:3]]}")

    out = dict(source=SRC, alma_project='2011.0.00025.SV (Band 2 science verification)', z_published=Z_PUBLISHED,
               n_spw=len(spws), aperture_pix=int(mask.sum()),
               coverage_ghz=[[float(w['freq'].min()), float(w['freq'].max())] for w in spec],
               results={k: {kk: vv for kk, vv in v.items() if kk not in ('lnz', 'post')} for k, v in res.items()},
               circle_aperture={kk: vv for kk, vv in resc.items() if kk not in ('lnz', 'post')},
               stress_noise=deg, stress_coverage=sub_res)
    (RESULTS / 'redshift_real_SPT0027.json').write_text(json.dumps(out, indent=1, default=float))
    plot(spec, img, csig, mask, res['SPT-like'], deg, sub_res, path=FIG / 'redshift_real_SPT0027.png')


def plot(spec, img, csig, mask, r, deg, sub_res, path):
    from matplotlib.ticker import NullFormatter
    fig = plt.figure(figsize=(15, 8.6), facecolor=SURF)
    gs = fig.add_gridspec(2, 3, height_ratios=[1.1, 1], width_ratios=[1, 1.3, 1.3])
    ax = fig.add_subplot(gs[0, 1:])
    for w in spec:
        ax.step(w['freq'], w['S'] / w['sig'], where='mid', color=COL[0], lw=0.8)
    zb = r['sled'][0]
    names, nus = line_list()
    ymax = max(np.nanmax(w['S'] / w['sig']) for w in spec)
    for n, nu in zip(names, nus):
        f = nu / (1 + zb)
        if 74 < f < 112:
            ax.axvline(f, color=COL[1], lw=1, ls='--')
            ax.text(f, ymax * 1.02, n, rotation=90, fontsize=7.5, color=COL[1], ha='center', va='bottom')
    for n, nu in EXTRA.items():
        f = nu / (1 + zb)
        if 74 < f < 112:
            ax.axvline(f, color=INK2, lw=0.8, ls=':')
            ax.text(f, ymax * 1.02, n, rotation=90, fontsize=7, color=INK2, ha='center', va='bottom')
    ax.set_xlim(74, 112); ax.set_ylim(-5, ymax * 1.35)
    ax.set_xlabel('observed frequency (GHz)', color=INK2, fontsize=9); ax.set_ylabel('S/N per 50 km/s channel', color=INK2, fontsize=9)
    ax.set_title(f'B. Real Band 2 spectrum of {SRC} (continuum removed); dashed: template lines at best z = {zb:.4f}; '
                 f'dotted: not in template', fontsize=9.5, color=INK, loc='left')

    ax = fig.add_subplot(gs[0, 0])
    ext = img.shape[0] * 0.14 / 2
    ax.imshow(img / csig, origin='lower', cmap='Greys', extent=[-ext, ext, -ext, ext], vmin=-2, vmax=10)
    ax.contour(mask.astype(float), levels=[0.5], colors=[COL[1]], extent=[-ext, ext, -ext, ext], linewidths=1.2)
    ax.set_xlim(-6, 6); ax.set_ylim(-6, 6)
    ax.set_xlabel('arcsec', color=INK2, fontsize=9)
    ax.set_title('A. All channels collapsed (S/N); orange = aperture', fontsize=9.5, color=INK, loc='left')

    ax = fig.add_subplot(gs[1, :2])
    labs = list(sub_res.keys())
    for i, lab in enumerate(labs):
        cands = sub_res[lab]['candidates']
        bottom = 0.0
        for z, m in cands[:4]:
            is_true = abs(np.log1p(z) - np.log1p(Z_PUBLISHED)) < CONF_WIN
            ax.barh(i, m, left=bottom, color=COL[2] if is_true else COL[1] if bottom == 0 else '#c9c8c3',
                    edgecolor=SURF, height=0.6)
            if m > 0.06:
                ax.text(bottom + m / 2, i, f'z={z:.2f}', ha='center', va='center', fontsize=7.5,
                        color='white' if (is_true or bottom == 0) else INK)
            bottom += m
    ax.set_yticks(range(len(labs))); ax.set_yticklabels(labs, fontsize=8); ax.invert_yaxis()
    ax.set_xlim(0, 1.0); ax.set_xlabel('posterior mass on each candidate redshift', color=INK2, fontsize=9)
    ax.set_title(f'C. Same spectrum, part of the coverage only: posterior mass per candidate (green = true z)',
                 fontsize=9.5, color=INK, loc='left')

    ax = fig.add_subplot(gs[1, 2])
    f = [q['noise_factor'] for q in deg]
    ax.plot(f, [q['confident_correct'] for q in deg], color=COL[2], marker='o', lw=2, label='confident & correct')
    ax.plot(f, [q['confidently_wrong'] for q in deg], color=COL[1], marker='s', lw=2, label='confidently wrong')
    ax.plot(f, [q['median_conf'] for q in deg], color=COL[0], marker='^', lw=1.5, ls='--', label='median confidence')
    ax.set_xscale('log'); ax.set_xticks(f); ax.set_xticklabels([str(x) for x in f]); ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_ylim(-0.02, 1.05); ax.set_xlabel('noise raised by factor (real spectrum + added noise)', color=INK2, fontsize=9)
    ax.set_title('D. Same spectrum + extra noise (30 draws)', fontsize=9.5, color=INK, loc='left'); ax.legend(fontsize=8, frameon=False)
    for a in fig.axes:
        a.set_facecolor(SURF)
        for sp_ in ('top', 'right'):
            a.spines[sp_].set_visible(False)
        a.tick_params(colors=INK2, labelsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=130, facecolor=SURF)


if __name__ == '__main__':
    main()
