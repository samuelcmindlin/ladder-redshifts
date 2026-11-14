"""A joint likelihood of the whole spectrum and the photometry, after Bing et al. (2024), as a comparison method.

Bing et al. 2024 (A&A 683, A232; arXiv:2403.00911, section 4.2) find the redshift of a dusty galaxy from a blind
spectral scan
without detecting lines one by one. At each trial redshift they predict the flux of every CO and [CI] line from
the galaxy's infrared luminosity, through published luminosity relations (their Table 7), draw each line as a
Gaussian of 500 km/s FWHM, and compare that model spectrum with the observed one by chi-square (their Eq. 5).
The likelihood exp(-chi2_spec / 2) is multiplied by the likelihood of the photometry at that redshift (Eq. 6).
It is the published method nearest to the engine of this repository when the engine is given a photometric
redshift. The differences: one set of line strengths with the amplitude fixed by the photometry, where the
engine has a bank of line-ratio templates with a free amplitude; one line width, where the engine has four;
lines exactly at their predicted frequencies, where the engine allows 150 km/s.

This is a re-implementation from the paper's text, not the authors' code (none is published), on SPT sources,
which the paper did not treat. Where the paper leaves a choice open, or these data force one, it is listed here:

  - Line luminosities: log L_IR = N log L'_line + A, with L' in K km/s pc2 and L_IR in solar luminosities. The
    paper prints "L_line = N x log(L_FIR) + A"; the coefficients are those of Greve et al. 2014 and Liu et al.
    2015, whose relations have the form used here. [CII] is left out: it is never inside 67-116 GHz below z = 8.
  - Infrared luminosity at each trial redshift: the paper fits galaxy templates to the photometry. Here a
    modified blackbody (emissivity index 2, optically thick below 100 micron rest, temperature free) is fitted
    to the published photometry of Reuter et al. 2020, Table D.1, at rest wavelengths of 50 micron or more, with
    10% added to every error for calibration, and integrated from 8 to 1000 micron. Non-detections enter as zero
    with their stated error. The luminosity is the apparent (lensed) one, as are the line fluxes.
  - The photometric term: the paper multiplies by exp(-chi2_SED / 2) of its template fit. Here that term is the
    published photometric redshift of the source (photoz.py, the 2020 values), the same prior the engine is
    given, so that the two methods differ only in how they read the spectrum. With no prior the method still
    uses the photometry, to set the line strengths; the engine with no prior uses none.
  - Confidence: the paper reports the maximum of the joint likelihood and defines no threshold. Here the
    likelihood is turned into a posterior on the engine's redshift grid with the engine's own code
    (Engine.posterior), so "confident" means the same for both: 95% of the posterior within 1800 km/s.
  - Spectra: continuum-subtracted, on the engine's 50 km/s grid, with the same noise per pixel.

Because a single reading of the paper could be a weak one, a small family of variants is run (VARIANTS):
  - the line strengths as predicted, and shifted by -1, -0.5, +0.5 and +1 times each relation's scatter (the
    test of the paper's Appendix D);
  - the scatter of each relation taken into the likelihood: every line's strength is integrated over a
    log-normal of the width in Table 7, centred on the prediction or 0.5 or 1 times the scatter below it. This is
    not in the paper; it is what its Table 7 allows, and it is the form that suits bright lines, for which a
    strength fixed in advance is wrong by many times the noise;
  - line widths of 300 and 500 km/s (the paper's own, and the one its Appendix A also tries).
The variant the engine is compared with is chosen on the noise-added draws of the development sample, by the
criterion used for the rule family (degrade.py), with both methods at the confidence threshold of 0.95: the most
confident and correct answers among the variants that are confidently wrong no more often than the engine with
the same prior; if none qualifies, the variant with the fewest wrong answers, and among those the most correct.
The choice is written to results/redshift_joint_likelihood.json, and a run on a second sample
(LADDER_REDSHIFTS_SAMPLE) reads it from there and does not choose again. The variant "as published" (predicted
strengths, 500 km/s) is always reported as well.

A threshold of 0.95 need not mean the same for two likelihoods, so both are also counted at every threshold of
degrade.ENGINE_THRESHOLDS: answers against wrong answers for the engine, for the chosen variant, and for the
best variant and threshold found with hindsight on the sample at hand.

The noise-added draws are those of degrade.py (same seeds), and the script stops if the engine's counts on them
differ from the stored ones.

Usage:
    python scripts/compare_joint_likelihood.py photometry spt_dndz_tab_d1.txt   copy the photometry rows used
    python scripts/compare_joint_likelihood.py run [N_draws]                    as observed and noise added
    python scripts/compare_joint_likelihood.py plot                             both samples together, and the figure
The table is posted at https://github.com/spt-smg/publicdata (directory reuter2020_final_dndz). With
LADDER_REDSHIFTS_SAMPLE=independent the first two act on the second sample.
Writes data/photometry.csv and results/redshift_joint_likelihood.json (per sample), and from both samples
results/redshift_joint_likelihood_pooled.json and figures/redshift_joint_likelihood.png.
"""
import csv
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np
from scipy.signal import fftconvolve
from scipy.special import logsumexp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ladder_redshifts import degrade as D  # noqa: E402
from ladder_redshifts.paths import DATA, RESULTS, ROOT, SAMPLE  # noqa: E402
from ladder_redshifts.batch_real import read_spectrum, to_grid_auto, read_truth  # noqa: E402
from ladder_redshifts.realspec import run_engine, engine_likelihood, engine_posterior, degrade  # noqa: E402
from ladder_redshifts.ladder import CO, CI10, CI21, DLN, PIX_KMS, TOL_LN, prior_spt  # noqa: E402
from ladder_redshifts.photoz import photoz_priors  # noqa: E402

# Table 7 of Bing et al. 2024: log L_IR = N log L'_line + A, and the scatter of each relation in dex.
# CO(1-0) to CO(3-2): Greve et al. 2014; CO(4-3) to CO(12-11): Liu et al. 2015; [CI]: Valentino et al. 2018.
TABLE7 = [('CO(1-0)', CO[0], 0.99, 1.90, 0.30), ('CO(2-1)', CO[1], 1.03, 1.60, 0.30),
          ('CO(3-2)', CO[2], 0.99, 2.10, 0.30), ('CO(4-3)', CO[3], 1.06, 1.49, 0.27),
          ('CO(5-4)', CO[4], 1.07, 1.71, 0.22), ('CO(6-5)', CO[5], 1.10, 1.79, 0.19),
          ('CO(7-6)', CO[6], 1.03, 2.62, 0.19), ('CO(8-7)', CO[7], 1.02, 2.82, 0.21),
          ('CO(9-8)', CO[8], 1.01, 3.10, 0.27), ('CO(10-9)', CO[9], 0.96, 3.67, 0.26),
          ('CO(11-10)', CO[10], 1.00, 3.51, 0.27), ('CO(12-11)', CO[11], 0.99, 3.83, 0.28),
          ('[CI](1-0)', CI10, 1.28, 0.00, 0.26), ('[CI](2-1)', CI21, 1.28, 0.61, 0.50)]
NAMES = [t[0] for t in TABLE7]
NU = np.array([t[1] for t in TABLE7])
SLOPE = np.array([t[2] for t in TABLE7])
ICPT = np.array([t[3] for t in TABLE7])
SCATTER = np.array([t[4] for t in TABLE7])
PAIR = (NAMES.index('CO(7-6)'), NAMES.index('[CI](2-1)'))      # the one pair close enough to overlap

WIDTHS = (500.0, 300.0)                         # km/s FWHM; the first is the paper's
OFFSETS = (0.0, -1.0, -0.5, 0.5, 1.0)           # shift of the predicted strengths, in units of each relation's scatter
SCATTER_OFFSETS = (0.0, -0.5, -1.0)             # centre of the log-normal when the scatter is in the likelihood
VARIANTS = ([dict(strength='fixed', offset=o, fwhm=w) for w in WIDTHS for o in OFFSETS]
            + [dict(strength='scatter', offset=o, fwhm=w) for w in WIDTHS for o in SCATTER_OFFSETS])
AS_PUBLISHED = 0                                # index in VARIANTS: predicted strengths, 500 km/s
PRIORS = ('uniform', 'photoz')                  # for the noise-added test; as observed, the survey prior too
THRESHOLDS = D.ENGINE_THRESHOLDS
CONF = D.CONF
KCONF = THRESHOLDS.index(CONF)
T_NODES = np.linspace(-6.0, 6.0, 61)            # nodes of the integral over a relation's scatter, in sigma
T_LOGW = -0.5 * T_NODES ** 2 - logsumexp(-0.5 * T_NODES ** 2)
U_NODES = np.linspace(-6.0, 6.0, 49)            # nodes across the data's own constraint on a bright line's strength
SNR_SWITCH = 4.0                                # above this measured S/N the second set of nodes is used

BANDS_UM = (3000.0, 2000.0, 1400.0, 870.0, 500.0, 350.0, 250.0, 160.0, 100.0)    # columns of Table D.1
PHOT_COLUMNS = ['source'] + [f'{k}_{int(b)}um' for b in BANDS_UM for k in ('flux_mjy', 'err_mjy')]
H_PLANCK, K_BOLTZ, C_LIGHT, L_SUN, MPC = 6.62607015e-34, 1.380649e-23, 2.99792458e8, 3.828e26, 3.0856775814913673e22
H0, OMEGA_M = 67.7, 0.31                        # flat LCDM, as in Reuter et al. 2020 (Planck 2015)
Z_COARSE = np.arange(0.4, 8.25, 0.05)           # the luminosity is fitted here and interpolated
T_GRID = np.arange(10.0, 160.25, 0.5)
OUT = RESULTS / 'redshift_joint_likelihood.json'
CHOICE = ROOT / 'results' / 'redshift_joint_likelihood.json'      # the development sample's file holds the choice


# ---------------------------------------------------------------- photometry and luminosity

def write_photometry(table):
    """Copy the rows of Table D.1 for the sources of this sample (those in truth.csv with a spectrum)."""
    rows = {}
    for line in Path(table).read_text().splitlines():
        if line.startswith('#') or not line.strip():
            continue
        p = line.split()
        rows[p[0]] = p[2:2 + 2 * len(BANDS_UM)]
    want = [s for s in sorted(read_truth()) if (DATA / 'spectra' / f'{s}.csv').exists()]
    out = []
    for s in want:
        if s not in rows:
            print(f'  {s}: not in the table'); continue
        if len(rows[s]) != 2 * len(BANDS_UM):
            raise SystemExit(f'{s}: expected {2 * len(BANDS_UM)} photometry columns, found {len(rows[s])}')
        out.append([s] + ['' if v.lower() == 'nan' else v for v in rows[s]])
    with open(DATA / 'photometry.csv', 'w', newline='') as fh:
        w = csv.writer(fh, lineterminator='\n')
        w.writerow(PHOT_COLUMNS); w.writerows(out)
    print(f'wrote {DATA / "photometry.csv"}: {len(out)} of {len(want)} sources')


def read_photometry():
    """{source: (flux, error)} in mJy per band of BANDS_UM; NaN where the table has no value."""
    path = DATA / 'photometry.csv'
    if not path.exists():
        return {}
    num = lambda v: float(v) if v else np.nan
    with open(path, newline='') as fh:
        return {r['source']: (np.array([num(r[f'flux_mjy_{int(b)}um']) for b in BANDS_UM]),
                              np.array([num(r[f'err_mjy_{int(b)}um']) for b in BANDS_UM]))
                for r in csv.DictReader(fh)}


def lum_dist_mpc(z):
    """Luminosity distance in Mpc (flat LCDM)."""
    zz = np.linspace(0.0, float(np.max(z)), 4001)
    inv_e = 1.0 / np.sqrt(OMEGA_M * (1 + zz) ** 3 + 1 - OMEGA_M)
    dc = np.concatenate([[0.0], np.cumsum(0.5 * (inv_e[1:] + inv_e[:-1]) * np.diff(zz))])
    return (1 + np.asarray(z, float)) * np.interp(z, zz, dc) * C_LIGHT / 1e3 / H0


def greybody(nu_rest_hz, temp):
    """Modified blackbody, emissivity index 2, optically thick below 100 micron (arbitrary units)."""
    nu0 = C_LIGHT / 100e-6
    return -np.expm1(-(nu_rest_hz / nu0) ** 2) * nu_rest_hz ** 3 / np.expm1(H_PLANCK * nu_rest_hz / (K_BOLTZ * temp))


_NU_INT = np.geomspace(C_LIGHT / 1000e-6, C_LIGHT / 8e-6, 600)                        # 8-1000 micron, rest
_SHAPE = greybody(_NU_INT[None, :], T_GRID[:, None])
_BOLO = (0.5 * (_SHAPE[:, 1:] + _SHAPE[:, :-1]) * np.diff(_NU_INT)).sum(1)            # integral for each temperature
_NU_FIR = np.geomspace(C_LIGHT / 500e-6, C_LIGHT / 42.5e-6, 400)                      # 42.5-500 micron, as Reuter+20
_SHAPE_FIR = greybody(_NU_FIR[None, :], T_GRID[:, None])
_BOLO_FIR = (0.5 * (_SHAPE_FIR[:, 1:] + _SHAPE_FIR[:, :-1]) * np.diff(_NU_FIR)).sum(1)


def infrared_luminosity(flux, err, z, fir=False):
    """Apparent infrared luminosity (solar) of the best modified blackbody at each redshift in z, and its
    temperature. 8-1000 micron rest, or 42.5-500 micron with fir=True (the range of Reuter et al. 2020)."""
    z = np.atleast_1d(np.asarray(z, float))
    nu_obs = C_LIGHT / (np.array(BANDS_UM) * 1e-6)
    f = np.where(np.isfinite(flux), flux, 0.0)                       # a non-detection: zero with its stated error
    e = np.sqrt(err ** 2 + (0.10 * f) ** 2)
    lum, temp = np.zeros(len(z)), np.zeros(len(z))
    dl = lum_dist_mpc(z) * MPC
    for i, zi in enumerate(z):
        use = np.isfinite(e) & (e > 0) & (np.array(BANDS_UM) / (1 + zi) >= 50.0)
        m = greybody(nu_obs[None, use] * (1 + zi), T_GRID[:, None])                   # temperature x band
        w = 1.0 / e[use] ** 2
        norm = (m * f[use] * w).sum(1) / (m * m * w).sum(1)
        chi2 = ((f[use] - norm[:, None] * m) ** 2 * w).sum(1)
        k = int(np.argmin(chi2))
        bolo = (_BOLO_FIR if fir else _BOLO)[k]
        # L = 4 pi D_L^2 x integral of S(nu_obs) d nu_obs; S in mJy = 1e-29 W m-2 Hz-1; d nu_obs = d nu_rest / (1+z)
        lum[i] = 4 * np.pi * dl[i] ** 2 * max(norm[k], 0.0) * 1e-29 * bolo / (1 + zi) / L_SUN
        temp[i] = T_GRID[k]
    return lum, temp


def line_peaks_mjy(lum, z, fwhm):
    """Peak flux density (mJy) of each line of TABLE7 for infrared luminosity lum at redshift z (arrays of one
    length), for a Gaussian line of the given FWHM. Returns redshift x line."""
    lprime = 10 ** ((np.log10(np.maximum(lum, 1.0))[:, None] - ICPT[None, :]) / SLOPE[None, :])   # K km/s pc2
    nu_obs = NU[None, :] / (1 + z[:, None])
    # Solomon et al. 1997: L' = 3.25e7 S dv nu_obs^-2 D_L^2 (1+z)^-3, with S dv in Jy km/s and D_L in Mpc
    sdv = lprime * nu_obs ** 2 * (1 + z[:, None]) ** 3 / (3.25e7 * lum_dist_mpc(z)[:, None] ** 2)
    return 1e3 * sdv / (1.0644670 * fwhm)


# ---------------------------------------------------------------- the likelihood

def profile(fwhm):
    s = fwhm / 2.3548 / PIX_KMS
    x = np.arange(-int(5 * s) - 1, int(5 * s) + 2)
    return x, np.exp(-0.5 * (x / s) ** 2)


class Joint:
    """The spectral likelihood of Bing et al. for one source: everything that does not depend on the noise draw."""

    def __init__(self, grid, lnz, flux, err):
        self.lnz = lnz
        z = np.exp(lnz) - 1
        P = np.rint((np.log(NU)[None, :] - lnz[:, None] - grid.lnmin) / DLN).astype(int)
        self.inr = (P >= 0) & (P < grid.npix)
        self.P = np.where(self.inr, P, 0)
        lum_c, _ = infrared_luminosity(flux, err, Z_COARSE)
        self.lum = np.exp(np.interp(z, Z_COARSE, np.log(np.maximum(lum_c, 1.0))))
        self.peaks = {w: line_peaks_mjy(self.lum, z, w) for w in WIDTHS}             # at the predicted strengths
        self.lag = int(round(np.log(NU[PAIR[1]] / NU[PAIR[0]]) / DLN))               # pixels between the close pair
        self.kern = {}
        for w in WIDTHS:
            x, g = profile(w)
            m = len(x) // 2 + self.lag
            xx = np.arange(-m, m + 1)
            s = w / 2.3548 / PIX_KMS
            self.kern[w] = (g, np.exp(-0.5 * (xx / s) ** 2) * np.exp(-0.5 * ((xx - self.lag) / s) ** 2))

    def loglike(self, d, sig, cov):
        """ln likelihood - ln likelihood of no lines, at every redshift, for each variant (variant x redshift).
        d: the spectrum in units of its noise; sig: noise per pixel (mJy); cov: observed pixels."""
        w1 = np.where(cov, d / np.where(cov, sig, 1.0), 0.0)         # flux / sigma^2
        w2 = np.where(cov, 1.0 / np.where(cov, sig, 1.0) ** 2, 0.0)
        out = np.empty((len(VARIANTS), len(self.lnz)))
        terms = {}
        for w in WIDTHS:
            g, g12 = self.kern[w]
            b0 = np.where(self.inr, fftconvolve(w1, g, mode='same')[self.P], 0.0)           # sum of model x data / var
            a0 = np.where(self.inr, np.clip(fftconvolve(w2, g * g, mode='same'), 0, None)[self.P], 0.0)   # model^2 / var
            # overlap of the two lines of the close pair: sum over pixels q of g(q - p_i) g(q - p_i - lag) / var
            x0 = np.clip(fftconvolve(w2, g12[::-1], mode='same'), 0, None)[self.P[:, PAIR[0]]]
            x0 = np.where(self.inr[:, PAIR[0]] & self.inr[:, PAIR[1]], x0, 0.0)
            terms[w] = (b0, a0, x0)
        for k, v in enumerate(VARIANTS):
            b0, a0, x0 = terms[v['fwhm']]
            amp = self.peaks[v['fwhm']] * 10 ** (v['offset'] * SCATTER[None, :])
            if v['strength'] == 'fixed':
                out[k] = (amp * b0 - 0.5 * amp ** 2 * a0).sum(1) - amp[:, PAIR[0]] * amp[:, PAIR[1]] * x0
            else:
                # each line's strength integrated over a log-normal of its relation's scatter (the overlap of the
                # close pair, under 1% of a line at these widths, is left out)
                out[k] = ln_marginal(amp * b0, amp ** 2 * a0, np.broadcast_to(SCATTER[None, :], amp.shape)).sum(1)
        return out

    def lines_at(self, z_true, d, sig, cov, fwhm=WIDTHS[0]):
        """For the lines inside the coverage at z_true: predicted and measured matched-filter S/N and their ratio."""
        i = int(np.argmin(np.abs(self.lnz - np.log1p(z_true))))
        g = self.kern[fwhm][0]
        w1 = np.where(cov, d / np.where(cov, sig, 1.0), 0.0)
        w2 = np.where(cov, 1.0 / np.where(cov, sig, 1.0) ** 2, 0.0)
        c1, c2 = fftconvolve(w1, g, mode='same'), np.clip(fftconvolve(w2, g * g, mode='same'), 0, None)
        out = []
        for j, name in enumerate(NAMES):
            p = self.P[i, j]
            if self.inr[i, j] and c2[p] > 0 and cov[p]:
                amp = self.peaks[fwhm][i, j]
                out.append(dict(line=name, predicted_snr=round(float(amp * np.sqrt(c2[p])), 2),
                                measured_snr=round(float(c1[p] / np.sqrt(c2[p])), 2),
                                measured_over_predicted=round(float(c1[p] / c2[p] / amp), 3)))
        return out


def ln_marginal(B, A, sdex):
    """ln of the integral over a line's strength factor s of p(s) exp(s B - s^2 A / 2), with log10 s normal of
    width sdex. B = sum of model x data / variance and A = sum of model^2 / variance at the predicted strength,
    so the exponent is the line's contribution to -chi2/2 relative to no line. Arrays of one shape; 0 where A = 0.

    As a function of s the exponent is a Gaussian centred on B/A with standard deviation 1/sqrt(A). For a line
    measured at high S/N that is far narrower than the spacing of any fixed set of nodes in log s, so the
    integral is then taken on nodes across that Gaussian; otherwise on nodes across the log-normal."""
    out = np.zeros(B.shape)
    act = A > 0
    snr = np.zeros(B.shape)
    snr[act] = B[act] / np.sqrt(A[act])
    hi = act & (snr > SNR_SWITCH)
    lo = act & ~hi
    if lo.any():
        s = 10 ** (T_NODES[None, :] * sdex[lo][:, None])
        out[lo] = logsumexp(s * B[lo][:, None] - 0.5 * s ** 2 * A[lo][:, None] + T_LOGW[None, :], axis=1)
    if hi.any():
        b, a, sd = B[hi], A[hi], sdex[hi]
        width = 1.0 / np.sqrt(a)
        s = (b / a)[:, None] + width[:, None] * U_NODES[None, :]
        pos = s > 0
        sp = np.where(pos, s, 1.0)
        lnp = -0.5 * (np.log10(sp) / sd[:, None]) ** 2 - np.log(sp * np.log(10.0) * sd[:, None] * np.sqrt(2 * np.pi))
        lnf = np.where(pos, lnp - 0.5 * U_NODES[None, :] ** 2, -np.inf)
        out[hi] = 0.5 * b * b / a + np.log(width * (U_NODES[1] - U_NODES[0])) + logsumexp(lnf, axis=1)
    return out


def dump(obj, path):
    """Write JSON with one-level indentation, but with every list of plain numbers on a single line."""
    import re
    text = json.dumps(obj, indent=1)
    flat = lambda m: '[' + re.sub(r'\s+', ' ', m.group(1)).strip() + ']'
    text = re.sub(r'\[\s*((?:-?\d+(?:\.\d+)?(?:e-?\d+)?,\s*)*-?\d+(?:\.\d+)?(?:e-?\d+)?)\s*\]', flat, text)
    path.write_text(text + '\n')


def correct(z, z_true):
    return abs(np.log1p(z) - np.log1p(z_true)) < TOL_LN


def variant_name(v):
    off = 'as predicted' if v['offset'] == 0 else f"{v['offset']:+g} sigma"
    if v['strength'] == 'scatter':
        return f"scatter in the likelihood, centred {off}, {v['fwhm']:.0f} km/s"
    return f"strengths {off}, {v['fwhm']:.0f} km/s"


# ---------------------------------------------------------------- as observed

def as_observed(phot):
    """Engine and every variant on every real spectrum that has a published redshift and photometry."""
    rows = []
    for src, t in sorted(read_truth().items()):
        path = DATA / 'spectra' / f'{src}.csv'
        pz = photoz_priors(src).get('photoz')
        if not path.exists() or src not in phot or pz is None:
            continue
        z_true = float(t['z_published'])
        g = to_grid_auto(read_spectrum(path))
        grid, d, sig, cov = g
        like = engine_likelihood(*g)
        e = like['e']
        J = Joint(grid, e.lnz, *phot[src])
        LJ = J.loglike(d, sig, cov)
        lines = 'one' if '(1 line)' in t['notes'] else 'none' if '(0 lines' in t['notes'] else 'several'
        i_true = int(np.argmin(np.abs(e.lnz - np.log1p(z_true))))
        row = dict(source=src, z_published=z_true, lines=lines, status=t.get('status', 'confirmed'),
                   log_lir_at_published=round(float(np.log10(J.lum[i_true])), 3),
                   lines_at_published=J.lines_at(z_true, d, sig, cov), engine={}, joint=[{} for _ in VARIANTS])
        for name, prior in (('uniform', None), ('spt', prior_spt), ('photoz', pz)):
            ze, ce = engine_posterior(like, prior)['sled']
            row['engine'][name] = dict(z=round(ze, 4), confidence=round(ce, 3), correct=bool(correct(ze, z_true)))
            for k in range(len(VARIANTS)):
                zj, cj, _ = e.posterior(LJ[k], prior=prior)
                row['joint'][k][name] = dict(z=round(float(zj), 4), confidence=round(float(cj), 3),
                                             correct=bool(correct(zj, z_true)))
        rows.append(row)
    return rows


def _count(answers, threshold=CONF):
    return dict(confident_correct=int(sum(a['correct'] and a['confidence'] >= threshold for a in answers)),
                confident_wrong=int(sum((not a['correct']) and a['confidence'] >= threshold for a in answers)),
                best_is_published=int(sum(a['correct'] for a in answers)))


def observed_summary(rows, chosen):
    """Counts over the real spectra, per group of sources and prior: confident and correct, confident and wrong,
    and best redshift is the published one, for the engine, the variant as published, and the chosen form (the
    survey prior uses the form chosen for no prior). Also the sources only one of the two settles."""
    out = {}
    for group in ('several', 'one', 'none', 'all'):
        rr = [r for r in rows if group == 'all' or r['lines'] == group]
        if not rr:
            continue
        g = dict(n=len(rr), sources=[r['source'] for r in rr], priors={})
        for p in ('uniform', 'spt', 'photoz'):
            ch = chosen['photoz' if p == 'photoz' else 'uniform']
            k, thr = ch['index'], ch['threshold']
            settled_e = [r['source'] for r in rr if r['engine'][p]['correct'] and r['engine'][p]['confidence'] >= CONF]
            settled_j = [r['source'] for r in rr if r['joint'][k][p]['correct'] and r['joint'][k][p]['confidence'] >= thr]
            g['priors'][p] = dict(
                engine=_count([r['engine'][p] for r in rr]),
                as_published=_count([r['joint'][AS_PUBLISHED][p] for r in rr]),
                chosen=_count([r['joint'][k][p] for r in rr], thr),
                engine_only=sorted(set(settled_e) - set(settled_j)), chosen_only=sorted(set(settled_j) - set(settled_e)),
                wrong_engine=[r['source'] for r in rr if not r['engine'][p]['correct'] and r['engine'][p]['confidence'] >= CONF],
                wrong_chosen=[r['source'] for r in rr if not r['joint'][k][p]['correct'] and r['joint'][k][p]['confidence'] >= thr],
                every_variant=[_count([r['joint'][j][p] for r in rr]) for j in range(len(VARIANTS))])
        out[group] = g
    return out


# ---------------------------------------------------------------- noise added

def qualifying_sources():
    """The multi-line sources the noise-added test can use: two published lines detected, the second at the
    highest target S/N or more (degrade.main; scripts/independent_noise_test.py), in the same order."""
    keep = []
    for src, z in D.multi_line_sources():
        r = run_engine(*to_grid_auto(read_spectrum(DATA / 'spectra' / f'{src}.csv')), prior=prior_spt)
        lines = D.published_lines(r['detections'], z)
        if len(lines) >= 2 and lines[1][2] >= max(D.TARGETS):
            keep.append((src, z, lines[1][1], lines[1][2]))
    return keep


def run_level(task):
    """All draws of one source at one noise level: the draws of degrade.run_level (same seed, same order)."""
    src, z_true, factor, n, seed, flux, err = task
    spec = read_spectrum(DATA / 'spectra' / f'{src}.csv')
    rng = np.random.default_rng(seed)
    priors = dict(uniform=None, photoz=photoz_priors(src)['photoz'])
    thr = np.array(THRESHOLDS)
    eng = {p: np.zeros((len(thr), 2), int) for p in PRIORS}                    # threshold x (correct, wrong)
    jnt = {p: np.zeros((len(VARIANTS), len(thr), 2), int) for p in PRIORS}
    top = {p: np.zeros(1 + len(VARIANTS), int) for p in PRIORS}                # best redshift is the published one
    J = None
    for _ in range(n):
        g = to_grid_auto(degrade(spec, factor, rng))
        grid, d, sig, cov = g
        like = engine_likelihood(*g)
        e = like['e']
        if J is None:
            J = Joint(grid, e.lnz, flux, err)
        LJ = J.loglike(d, sig, cov)
        for p, prior in priors.items():
            ze, ce = engine_posterior(like, prior)['sled']
            ok = correct(ze, z_true)
            eng[p][:, 0 if ok else 1] += ce >= thr
            top[p][0] += ok
            for k in range(len(VARIANTS)):
                zj, cj, _ = e.posterior(LJ[k], prior=prior)
                ok = correct(zj, z_true)
                jnt[p][k, :, 0 if ok else 1] += cj >= thr
                top[p][1 + k] += ok
    return src, {p: dict(engine=eng[p].tolist(), joint=jnt[p].tolist(), top=top[p].tolist()) for p in PRIORS}


def check_same_draws(sources, counts, n):
    """Stop unless the engine's counts on these draws equal those stored by degrade.py for the same sample."""
    path = RESULTS / 'redshift_degradation.json'
    if not path.exists():
        print('  (no stored noise-added test to check the draws against)'); return False
    o = json.loads(path.read_text())
    if o['n_draws'] != n or list(o['sources']) != list(sources):
        print(f"  (the stored test used {o['n_draws']} draws, this run {n}: draws not checked)"); return False
    for src in sources:
        for k, q in enumerate(o['sources'][src]['levels']):
            for p in PRIORS:
                got = counts[src][k][p]['engine'][KCONF]
                want = [round(q[p]['engine']['confident_correct'] * n), round(q[p]['engine']['confidently_wrong'] * n)]
                if got != want:
                    raise SystemExit(f'{src} level {k} prior {p}: engine counts {got} differ from the stored {want}')
    return True


def choose(a_e, a_j):
    """The variant to compare with, per prior, by the criterion in the module docstring (both at CONF).
    a_e[p]: source x level x threshold x 2; a_j[p]: source x level x variant x threshold x 2 (correct, wrong)."""
    out = {}
    for p in PRIORS:
        engine_wrong = int(a_e[p][:, :, KCONF, 1].sum())
        tot = a_j[p][:, :, :, KCONF, :].sum((0, 1))                 # variant x 2
        ok, wrong = tot[:, 0], tot[:, 1]
        eligible = wrong <= engine_wrong
        if eligible.any():
            how, k = 'no more wrong than the engine', int(np.argmax(np.where(eligible, ok, -1)))
        else:
            how, k = 'none qualified: fewest wrong', int(np.argmax(np.where(wrong == wrong.min(), ok, -1)))
        out[p] = dict(index=k, threshold=CONF, threshold_index=KCONF, variant=VARIANTS[k],
                      name=variant_name(VARIANTS[k]), how=how, engine_wrong_allowed=engine_wrong,
                      eligible_variants=int(eligible.sum()), variants=len(VARIANTS))
    return out


def sign_test(ahead, behind):
    """Two-sided exact sign test for `ahead` wins and `behind` losses (ties left out)."""
    from scipy.stats import binomtest
    n = ahead + behind
    return 1.0 if n == 0 else float(binomtest(ahead, n, 0.5).pvalue)


def noise_summary(counts, chosen, names, n):
    """Totals, rates per noise level (mean over sources) and per source, for the engine at CONF, the variant as
    published at CONF and the chosen form at its threshold; every variant at every threshold; and how many
    sources the engine is ahead of or behind the chosen form at three noise levels."""
    lv = lambda arr: [round(float(x), 4) for x in arr]
    out = {}
    for p in PRIORS:
        a_e, a_j, a_t = (np.array(counts[p][k]) for k in ('engine', 'joint', 'top'))
        ch = chosen[p]

        def block(c, best):                                         # c: source x level x 2
            return dict(total_correct=int(c[..., 0].sum()), total_wrong=int(c[..., 1].sum()),
                        best_is_published=int(best.sum()),
                        correct_by_level=lv((c[..., 0] / n).mean(0)), wrong_by_level=lv((c[..., 1] / n).mean(0)),
                        correct_by_source={s: lv(c[i, :, 0] / n) for i, s in enumerate(names)})
        e = block(a_e[:, :, KCONF, :], a_t[:, :, 0])
        c = block(a_j[:, :, ch['index'], ch['threshold_index'], :], a_t[:, :, 1 + ch['index']])
        ap = block(a_j[:, :, AS_PUBLISHED, KCONF, :], a_t[:, :, 1 + AS_PUBLISHED])
        per_source = []
        for target in (4.0, 3.5, 3.0):
            k = list(D.TARGETS).index(target)
            de = a_e[:, k, KCONF, 0] - a_j[:, k, ch['index'], ch['threshold_index'], 0]
            per_source.append(dict(target=target, engine_ahead=int((de > 0).sum()), engine_behind=int((de < 0).sum()),
                                   level=int((de == 0).sum()),
                                   sign_test_p=round(sign_test(int((de > 0).sum()), int((de < 0).sum())), 4)))
        tot = a_j.sum((0, 1))                                       # variant x threshold x 2
        operating = []
        for i, t in enumerate(THRESHOLDS):
            ec, ew = int(a_e[:, :, i, 0].sum()), int(a_e[:, :, i, 1].sum())
            elig = tot[..., 1] <= ew                                # with hindsight: the best variant and threshold
            best = None                                             # that is wrong no more often than the engine here
            if elig.any():
                k, j = np.unravel_index(int(np.argmax(np.where(elig, tot[..., 0], -1))), elig.shape)
                best = dict(name=variant_name(VARIANTS[int(k)]), index=int(k), threshold=THRESHOLDS[int(j)],
                            correct=int(tot[k, j, 0]), wrong=int(tot[k, j, 1]))
            operating.append(dict(threshold=t, engine_correct=ec, engine_wrong=ew,
                                  chosen_correct=int(tot[ch['index'], i, 0]), chosen_wrong=int(tot[ch['index'], i, 1]),
                                  best_with_hindsight_no_more_wrong=best))
        out[p] = dict(engine=e, chosen=dict(c, index=ch['index'], threshold=ch['threshold'], name=ch['name']),
                      as_published=ap, engine_against_chosen_by_source=per_source, operating_points=operating,
                      variants=[dict(index=k, name=variant_name(v), **v, best_is_published=int(a_t[:, :, 1 + k].sum()),
                                     by_threshold=[dict(threshold=t, correct=int(tot[k, i, 0]), wrong=int(tot[k, i, 1]))
                                                   for i, t in enumerate(THRESHOLDS)])
                                for k, v in enumerate(VARIANTS)])
    return out


def noise_added(phot, n):
    srcs = qualifying_sources()
    names = [s for s, *_ in srcs]
    tasks = [(s, z, snr2 / target, n, [D.SEED, i, k], *phot[s])
             for i, (s, z, _, snr2) in enumerate(srcs) for k, target in enumerate(D.TARGETS)]
    got = {s: [] for s in names}
    with mp.get_context('fork').Pool(2) as pool:
        for src, out in pool.imap(run_level, tasks):
            got[src].append(out)
    same = check_same_draws(names, got, n)
    counts = {p: {key: [[got[s][k][p][key] for k in range(len(D.TARGETS))] for s in names]
                  for key in ('engine', 'joint', 'top')} for p in PRIORS}
    if SAMPLE:
        chosen = json.loads(CHOICE.read_text())['noise_added']['chosen']
        how = f'fixed beforehand on the development sample (results/{CHOICE.name})'
    else:
        chosen = choose({p: np.array(counts[p]['engine']) for p in PRIORS},
                        {p: np.array(counts[p]['joint']) for p in PRIORS})
        how = 'chosen on these draws'
    return dict(draws_per_level=n, sources=names, targets=list(D.TARGETS), thresholds=list(THRESHOLDS),
                total_draws=n * len(names) * len(D.TARGETS), same_draws_as_degrade=same, chosen=chosen, choice=how,
                summary=noise_summary(counts, chosen, names, n),
                counts=dict(note='engine: source x level x threshold x (correct, wrong); joint: source x level x '
                                 'variant x threshold x (correct, wrong); top: source x level x (engine, then each '
                                 'variant), draws whose best redshift is the published one', **counts))


# ---------------------------------------------------------------- main

def run(n):
    t0 = time.time()
    phot = read_photometry()
    if not phot:
        raise SystemExit(f'{DATA / "photometry.csv"} not found: run the "photometry" step first')
    rows = as_observed(phot)
    na = noise_added(phot, n)
    ratios = [q['measured_over_predicted'] for r in rows for q in r['lines_at_published'] if q['measured_snr'] >= 5]
    out = dict(sample=SAMPLE or 'development', method='joint spectrum and photometry likelihood after Bing et al. '
               '2024 (arXiv:2403.00911), re-implemented; see scripts/compare_joint_likelihood.py',
               table7=[dict(line=t[0], slope=t[2], intercept=t[3], scatter_dex=t[4]) for t in TABLE7],
               variants=[dict(index=k, name=variant_name(v), **v) for k, v in enumerate(VARIANTS)],
               as_published=AS_PUBLISHED, confidence_threshold=CONF,
               as_observed=dict(
                   summary=observed_summary(rows, na['chosen']),
                   lines_measured_over_predicted=dict(
                       n=len(ratios), median_dex=round(float(np.median(np.log10(ratios))), 3),
                       scatter_dex=round(float(np.std(np.log10(ratios))), 3),
                       note='lines of the published redshift measured at 5 sigma or more, 500 km/s profile'),
                   sources=rows),
               noise_added=na)
    dump(out, OUT)
    report(out)
    print(f'[{time.time() - t0:.0f}s] wrote {OUT}')


def report(o):
    na, ch = o['noise_added'], o['noise_added']['chosen']
    q = o['as_observed']['lines_measured_over_predicted']
    print(f"\n{o['sample']} sample. Lines measured / predicted: median {q['median_dex']:+.2f} dex, scatter "
          f"{q['scatter_dex']:.2f} dex ({q['n']} lines)")
    print(f"form compared with ({na['choice']}; same draws as degrade: {na['same_draws_as_degrade']}):")
    for p in PRIORS:
        print(f"  {p}: {ch[p]['name']} [{ch[p]['how']}; {ch[p]['eligible_variants']} of {ch[p]['variants']} "
              f"variants qualify]")
    fmt = lambda c: f"{c['confident_correct']}/{c['confident_wrong']}/{c['best_is_published']}"
    print('As observed: confident and correct / confident and wrong / best is published')
    for group, g in o['as_observed']['summary'].items():
        print(f"  {group} ({g['n']} sources)")
        for p, c in g['priors'].items():
            print(f"    {p:8s} engine {fmt(c['engine'])}   as published {fmt(c['as_published'])}   chosen "
                  f"{fmt(c['chosen'])}   engine only {c['engine_only']}   chosen only {c['chosen_only']}   "
                  f"wrong: engine {c['wrong_engine']} chosen {c['wrong_chosen']}")
    print(f"Noise added: {na['total_draws']} draws, {len(na['sources'])} sources; confident and correct by second-line "
          f"S/N {na['targets']}; totals correct/wrong; draws whose best redshift is the published one")
    for p in PRIORS:
        s = na['summary'][p]
        for name in ('engine', 'chosen', 'as_published'):
            b = s[name]
            print(f"  {p:8s} {name:13s} {' '.join(f'{x:.2f}' for x in b['correct_by_level'])}  "
                  f"{b['total_correct']}/{b['total_wrong']}  best {b['best_is_published']}")
        print('           engine against chosen, per source: ' + '; '.join(
            f"{x['target']}: ahead {x['engine_ahead']}, behind {x['engine_behind']} (p={x['sign_test_p']})"
            for x in s['engine_against_chosen_by_source']))
        print('           by threshold, engine | chosen | best with hindsight and no more wrong than the engine:')
        for x in s['operating_points']:
            b = x['best_with_hindsight_no_more_wrong']
            print(f"             {x['threshold']}: {x['engine_correct']}/{x['engine_wrong']} | {x['chosen_correct']}/"
                  f"{x['chosen_wrong']} | " + ('none' if b is None else
                                               f"{b['correct']}/{b['wrong']} ({b['name']} at {b['threshold']})"))
        for v in s['variants']:
            print(f"           {v['name']:55s} " + ' '.join(f"{x['correct']}/{x['wrong']}" for x in v['by_threshold'])
                  + f"  best {v['best_is_published']}")


# ---------------------------------------------------------------- both samples together

SAMPLE_FILES = [('development', ROOT / 'results' / 'redshift_joint_likelihood.json'),
                ('independent', ROOT / 'results' / 'independent' / 'redshift_joint_likelihood.json')]
POOLED = ROOT / 'results' / 'redshift_joint_likelihood_pooled.json'


def pooled():
    """Both samples together, from their results files: the real spectra as observed, the noise-added totals,
    and for how many sources the engine is ahead of or behind the chosen variant at three noise levels.
    Writes results/redshift_joint_likelihood_pooled.json."""
    oo = [(name, json.loads(f.read_text())) for name, f in SAMPLE_FILES if f.exists()]
    out = dict(samples=[name for name, _ in oo], confidence_threshold=CONF, as_observed={}, noise_added={})
    for p in ('uniform', 'spt', 'photoz'):
        g = [o['as_observed']['summary']['all']['priors'][p] for _, o in oo]
        out['as_observed'][p] = dict(
            sources=sum(o['as_observed']['summary']['all']['n'] for _, o in oo),
            **{m: {k: sum(q[m][k] for q in g) for k in ('confident_correct', 'confident_wrong', 'best_is_published')}
               for m in ('engine', 'chosen', 'as_published')},
            engine_only=sum((q['engine_only'] for q in g), []), chosen_only=sum((q['chosen_only'] for q in g), []))
    for p in PRIORS:
        s = [o['noise_added']['summary'][p] for _, o in oo]
        levels = []
        for target in (4.0, 3.5, 3.0):
            ahead = behind = level = 0
            mean = dict(engine=[], chosen=[])
            for _, o in oo:
                na = o['noise_added']
                k = na['targets'].index(target)
                ch = na['chosen'][p]
                e = np.array(na['counts'][p]['engine'])[:, k, KCONF, 0]
                j = np.array(na['counts'][p]['joint'])[:, k, ch['index'], ch['threshold_index'], 0]
                ahead += int((e > j).sum()); behind += int((e < j).sum()); level += int((e == j).sum())
                mean['engine'] += list(e / na['draws_per_level']); mean['chosen'] += list(j / na['draws_per_level'])
            levels.append(dict(target=target, engine_ahead=ahead, engine_behind=behind, level=level,
                               sign_test_p=round(sign_test(ahead, behind), 4),
                               engine_mean=round(float(np.mean(mean['engine'])), 4),
                               chosen_mean=round(float(np.mean(mean['chosen'])), 4)))
        out['noise_added'][p] = dict(
            draws=sum(o['noise_added']['total_draws'] for _, o in oo),
            sources=sum(len(o['noise_added']['sources']) for _, o in oo),
            **{m: dict(correct=sum(q[m]['total_correct'] for q in s), wrong=sum(q[m]['total_wrong'] for q in s),
                       best_is_published=sum(q[m]['best_is_published'] for q in s))
               for m in ('engine', 'chosen', 'as_published')},
            engine_against_chosen_by_source=levels)
    dump(out, POOLED)
    print(f'wrote {POOLED}')
    for p, q in out['as_observed'].items():
        print(f"as observed, {p}: {q['sources']} sources; confident and correct / wrong: engine "
              f"{q['engine']['confident_correct']}/{q['engine']['confident_wrong']}, chosen "
              f"{q['chosen']['confident_correct']}/{q['chosen']['confident_wrong']}, as published "
              f"{q['as_published']['confident_correct']}/{q['as_published']['confident_wrong']}; engine only "
              f"{q['engine_only']}, chosen only {q['chosen_only']}")
    for p, q in out['noise_added'].items():
        print(f"noise added, {p}: {q['draws']} draws, {q['sources']} sources; correct/wrong: engine "
              f"{q['engine']['correct']}/{q['engine']['wrong']}, chosen {q['chosen']['correct']}/{q['chosen']['wrong']}, "
              f"as published {q['as_published']['correct']}/{q['as_published']['wrong']}")
        for x in q['engine_against_chosen_by_source']:
            print(f"   second line at {x['target']}: engine {x['engine_mean']:.3f}, chosen {x['chosen_mean']:.3f}; engine "
                  f"ahead for {x['engine_ahead']}, behind for {x['engine_behind']}, level {x['level']} "
                  f"(p = {x['sign_test_p']})")
    return out


# ---------------------------------------------------------------- figure

def plot():
    """Noise-added test, both samples: the engine, the chosen variant and the variant as published."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    ink, ink2, grid_c, surf = D.INK, D.INK2, D.GRID_C, D.SURF
    present = [(name.capitalize() + ' sample', json.loads(f.read_text())) for name, f in SAMPLE_FILES if f.exists()]
    if not present:
        raise SystemExit('no results to plot: run the "run" step first')
    pooled()
    titles = dict(uniform='no redshift prior', photoz="each source's photometric redshift as the prior")
    fig, axs = plt.subplots(len(present), 2, figsize=(12.4, 4.7 * len(present) + 0.5), facecolor=surf, sharey=True,
                            squeeze=False)
    for row, (name, o) in zip(axs, present):
        na = o['noise_added']
        x = na['targets']
        for ax, p in zip(row, PRIORS):
            s = na['summary'][p]
            scatter = VARIANTS[s['chosen']['index']]['strength'] == 'scatter'
            series = [('engine', D.COL[0], 'o', '-', s['engine']),
                      ('joint likelihood, ' + ('scatter of the relations included' if scatter else 'line strengths fixed')
                       + ' (the variant chosen)', D.COL[4], 'v', '-', s['chosen']),
                      ('joint likelihood, line strengths as published', D.COL[4], 'v', '--', s['as_published'])]
            for label, colr, mk, ls, b in series:
                dashed = ls == '--'
                ax.plot(x, b['correct_by_level'], color=colr, marker=mk, ms=7, lw=2, ls=ls,
                        mec=colr if dashed else surf, mfc=surf if dashed else colr, mew=1.4, label=label, clip_on=False)
            ax.set_title(f"{name}, {len(na['sources'])} sources: {titles[p]}", fontsize=10, color=ink, loc='left')
            ax.set_xticks(x); ax.set_ylim(0, 1.03)
            ax.grid(axis='y', color=grid_c, lw=0.8)
            ax.text(0.985, 0.66, f"a different redshift returned with confidence,\nof {na['total_draws']:,} draws: "
                    f"engine {s['engine']['total_wrong']}, chosen variant {s['chosen']['total_wrong']},\n"
                    f"as published {s['as_published']['total_wrong']:,}", transform=ax.transAxes, fontsize=8,
                    color=ink2, ha='right', va='top')
            D._style(ax)
        row[0].set_ylabel('fraction of noise draws in which the published\nredshift is returned with confidence',
                          fontsize=9, color=ink2)
    for ax in axs[-1]:
        ax.set_xlabel('matched-filter S/N of the second-brightest line after adding noise', fontsize=9, color=ink2)
    h, lab = axs[0, 0].get_legend_handles_labels()
    top = 1 - 0.38 / fig.get_figheight()
    fig.legend(h, lab, fontsize=8.5, frameon=False, loc='upper left', ncol=3, labelcolor=ink, handlelength=3.4,
               bbox_to_anchor=(0.006, top), columnspacing=1.6)
    fig.suptitle('The engine against a joint spectrum-and-photometry likelihood after Bing et al. (2024), re-implemented: '
                 'real spectra with noise added', fontsize=10.5, color=ink, x=0.012, y=0.995, ha='left')
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.62 / fig.get_figheight()))
    path = ROOT / 'figures' / 'redshift_joint_likelihood.png'
    fig.savefig(path, dpi=130, facecolor=surf)
    plt.close(fig)
    print(f'wrote {path}')


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else ''
    if cmd == 'photometry' and len(sys.argv) == 3:
        write_photometry(sys.argv[2])
    elif cmd == 'run':
        run(int(sys.argv[2]) if len(sys.argv) > 2 else 100)
    elif cmd == 'plot':
        plot()
    else:
        raise SystemExit(__doc__)
