"""Ladder redshifts as pitch detection: line catalogue, spectrum simulator and estimators.

The CO rotational ladder is a near-harmonic series: nu(J -> J-1) ~ J x 115.271 GHz.
After redshift its "fundamental" is 115.271/(1+z) GHz, and a spectral scan sees only a few
rungs. Finding z is therefore missing-fundamental pitch estimation, and its typical failure is
an octave-type error: (1+z_est)/(1+z_true) = J'/J.

Everything works on a log-frequency grid with constant velocity pixels (like SDSS / constant-Q),
so a redshift is a pure pixel shift and every line has the same width in pixels.

Estimators (all share one matched-filter front end):
  * naive_spacing  - detections >= 5 sigma, assume adjacent CO rungs; f0 = smallest spacing
                     ("pitch = partial spacing"). No estimate from a single line.
  * sudoku_rule    - enumerate (detection x transition) solutions, keep those explaining every
                     detection, reject those predicting a strong line where none is seen
                     (a simplified High-z Sudoku / gcd-style rule; Bakx & Dannerbauer 2022).
  * guided_rule    - sudoku_rule, then a look for a weaker line at the positions each surviving solution
                     predicts (what an observer does by eye); used by degrade.py.
  * xcorr_redshift - best redshift of a template cross-correlation, re-implemented after Yun et al. (2015);
                     used by degrade.py.
  * flat           - generalized-likelihood harmonic sum with EQUAL weights on the S/N of all covered
                     ladder lines (no line ratios, no noise weighting); a simplified stand-in for
                     line-sum methods such as Lupu et al. 2012, not their full procedure.
  * sled           - the same likelihood with SLED-shaped ("timbre") weights, averaged over a grid of
                     excitation templates, line widths and [CI] strength.
  * sled+photz     - sled with a photometric-redshift prior.

The likelihood statistic for template w is T = max(0, sum_i w_i y_i) / sqrt(w^T C w), where
y_i is the matched-filter S/N at the predicted position of line i and the sums run over lines
inside the frequency coverage. C_ij is the overlap of the matched filters at lines i and j
(C_ii = 1; nonzero off the diagonal only for lines closer than the filter width, such as
CO(7-6) and [CI](2-1), 998 km/s apart); without it one broad line seen from two predicted
positions would be counted twice. T^2/2 is the log generalized likelihood ratio for a line set
of known shape and unknown amplitude in white noise. Because the denominator counts every
covered line, a template that predicts a strong line where the data are empty is penalized
automatically: that is the "undertone check" of pitch detection, in likelihood form.

Each line may sit up to DV_MAX (150 km/s) from its predicted position: y_i is the largest filter
output within that range, i.e. the statistic is also maximised over one velocity offset per line
(except that CO(7-6) and [CI](2-1) keep their separation when both are covered).
Real lines of one galaxy differ by tens of km/s (kinematics, differential lensing), and the
ladder's anharmonicity separates octave-family readings of the same lines by 29-122 km/s; without
the tolerance, a bright spectrum would be "decided" by such offsets (docs/development-notes.md, item 21).

"""
import numpy as np
from scipy.ndimage import maximum_filter1d
from scipy.signal import fftconvolve
from scipy.special import logsumexp

C_KMS = 299792.458
PIX_KMS = 50.0                      # velocity pixel
DLN = PIX_KMS / C_KMS               # pixel in ln(nu)
TOL_LN = 0.003                      # "correct" if |ln(1+z_est) - ln(1+z_true)| < 0.003 (900 km/s)
CONF_WIN = 0.006                    # confidence = posterior mass within +/-0.006 (1800 km/s) of the best z, so that
                                    # blend partners (CO(7-6) vs [CI](2-1), 1000 km/s apart) are counted together
DV_MAX = 150.0                      # each line may sit this far (km/s) from its predicted position: wider than the
                                    # anharmonic mismatch between octave-family readings (29-122 km/s) and than
                                    # typical velocity offsets between the lines of one galaxy

# Rest frequencies (GHz). CO J=1-7 from ASPECS Table 1; J=8-13 standard catalogue values,
# cross-checked here against B0 = 57635.968 MHz, D0 = 0.1835 MHz to <0.01 MHz.
CO = np.array([115.271202, 230.538000, 345.795990, 461.040768, 576.267931, 691.473076,
               806.651806, 921.799700, 1036.912393, 1151.985452, 1267.014486, 1381.995105,
               1496.922909])
CI10, CI21, CII = 492.160651, 809.341970, 1900.536900
NCO = len(CO)


def line_list(include_cii=False):
    names = [f'CO({J}-{J-1})' for J in range(1, NCO + 1)] + ['[CI](1-0)', '[CI](2-1)']
    nus = list(CO) + [CI10, CI21]
    if include_cii:
        names.append('[CII]'); nus.append(CII)
    return names, np.array(nus)


def sled_flux(Jc, r_ci, r_ci21=0.45, f_cii=20.0, include_cii=False, scatter_dex=0.0, rng=None):
    """Velocity-integrated line fluxes relative to CO(1-0).

    CO brightness-temperature ratios r_J1 = exp(-(J-1)/Jc): Jc ~ 3.7 approximates the
    Bothwell+13 SMG average (0.84, 0.52, 0.41, 0.32, 0.21, 0.18 for J = 2..7); Jc ~ 6 matches
    the more excited Birkin+21 r41 = 0.63. Flux ratio = r x (nu/nu_1)^2.
    [CI](1-0): L' ratio r_ci to CO(1-0) (Walter+11: 0.29 +/- 0.12).
    [CI](2-1)/[CI](1-0) L' ratio r_ci21 = 0.45 (ILLUSTRATIVE; not verified).
    [CII] flux = f_cii x CO(4-3) flux (ILLUSTRATIVE; ~20 in SPT0452-50).
    """
    J = np.arange(1, NCO + 1)
    r = np.exp(-(J - 1) / Jc)
    if scatter_dex and rng is not None:
        r = r * 10 ** rng.normal(0, scatter_dex, NCO); r[0] = 1.0
    f = list(r * J ** 2)
    sc = (lambda: 10 ** rng.normal(0, scatter_dex)) if (scatter_dex and rng is not None) else (lambda: 1.0)
    f.append(r_ci * (CI10 / CO[0]) ** 2 * sc())
    f.append(r_ci * r_ci21 * (CI21 / CO[0]) ** 2 * sc())
    if include_cii:
        f.append(f_cii * f[3] * sc())
    return np.array(f)


def kernel(fwhm_kms):
    """Unit-norm Gaussian line profile on the pixel grid."""
    s = fwhm_kms / 2.3548 / PIX_KMS
    x = np.arange(-int(5 * s) - 1, int(5 * s) + 2)
    k = np.exp(-0.5 * (x / s) ** 2)
    return k / np.sqrt((k ** 2).sum())


def kernel_overlap(fwhm_kms, lag_px):
    """Correlation of two matched-filter outputs lag_px pixels apart (white noise): the kernel's autocorrelation."""
    k = kernel(fwhm_kms)
    lag = int(round(abs(lag_px)))
    if lag >= len(k):
        return 0.0
    return float(np.dot(k[lag:], k[:len(k) - lag]))


class Grid:
    """Log-frequency pixel grid spanning [nu_lo, nu_hi] GHz with a coverage mask."""

    def __init__(self, nu_lo, nu_hi):
        self.lnmin = np.log(nu_lo)
        self.npix = int(np.ceil((np.log(nu_hi) - self.lnmin) / DLN)) + 1
        self.lnnu = self.lnmin + DLN * np.arange(self.npix)

    def mask(self, windows):
        m = np.zeros(self.npix, bool)
        for w in windows:
            lo, hi = w[0], w[1]
            m |= (self.lnnu >= np.log(lo)) & (self.lnnu <= np.log(hi))
        return m

    def noise(self, windows):
        """Relative flux noise per pixel from windows given as (lo, hi, sigma)."""
        s = np.ones(self.npix)
        for w in windows:
            if len(w) > 2:
                s[(self.lnnu >= np.log(w[0])) & (self.lnnu <= np.log(w[1]))] = w[2]
        return s

    def pix(self, nu):
        return (np.log(nu) - self.lnmin) / DLN


def simulate(grid, z, line_amp, fwhm_kms, rng, include_cii=False, sigma=None):
    """Spectrum in units of the local per-pixel noise (an S/N spectrum).
    line_amp: line strengths (aligned with line_list) in flux units; a line's matched-filter
    S/N is line_amp / sigma at its position (sigma = relative noise per pixel, default 1, so
    by default line_amp IS the matched-filter S/N for a filter matched to the true width)."""
    names, nus = line_list(include_cii)
    s = fwhm_kms / 2.3548 / PIX_KMS
    knorm = 1.0 / np.sqrt(np.sqrt(np.pi) * s)        # unit-norm continuous Gaussian amplitude
    d = rng.normal(0, 1, grid.npix)
    x = np.arange(grid.npix)
    for nu, amp in zip(nus, line_amp):
        p = grid.pix(nu / (1 + z))
        if -5 * s < p < grid.npix + 5 * s and amp != 0:
            lo, hi = max(0, int(p - 5 * s)), min(grid.npix, int(p + 5 * s) + 2)
            sig = 1.0 if sigma is None else sigma[min(max(int(round(p)), 0), grid.npix - 1)]
            d[lo:hi] += (amp / sig) * knorm * np.exp(-0.5 * ((x[lo:hi] - p) / s) ** 2)
    return d


BANK = (250.0, 450.0, 800.0, 1200.0)  # matched-filter widths (km/s FWHM)


DETECT_FWHM = 450.0                   # profile used for line detection (detect, sudoku_rule, naive_spacing)


def detection_index(bank=BANK):
    """Row of a matched_filter output to use for detection: the profile closest to DETECT_FWHM."""
    return int(np.argmin(np.abs(np.asarray(bank, float) - DETECT_FWHM)))


def matched_filter(d, bank=BANK, cover=None):
    """Matched-filter S/N for each line profile in `bank` (one row per profile).

    cover: boolean mask of the observed pixels. When given (as every analysis in this package does), only
    observed pixels enter and the output keeps unit noise near edges and gaps:
        y(p) = sum_j k_j d(p+j) c(p+j) / sqrt(sum_j k_j^2 c(p+j)),
    so data outside the coverage cannot influence the result. A line centred at p then has its expected output
    reduced by filter_efficiency(cover, width)[p]; Engine applies that factor to the template weights.
    `cover` is required, so that data outside the coverage can never enter by accident; for a spectrum observed
    everywhere pass np.ones(len(d), bool)."""
    if cover is None:
        raise TypeError('matched_filter needs the coverage mask (cover=...); see Engine.filter')
    cover = np.asarray(cover, bool)
    dc = np.where(cover, d, 0.0)
    out = []
    for w in bank:
        k = kernel(w)
        num = fftconvolve(dc, k, mode='same')
        den = _efficiency(cover, k)
        out.append(np.where(cover & (den > 1e-6), num / np.where(den > 1e-6, den, 1.0), 0.0))
    return np.stack(out)


def _efficiency(cover, k):
    e2 = fftconvolve(np.asarray(cover, float), k * k, mode='same')
    return np.sqrt(np.clip(e2, 0.0, 1.0))


def filter_efficiency(cover, fwhm_kms):
    """Fraction of a unit line's matched-filter S/N that the coverage retains, per pixel: sqrt of the observed
    share of the profile's energy (1 inside a window, about 0.7 at a window edge, 0 outside)."""
    return _efficiency(cover, kernel(fwhm_kms))


def templates_grid(include_cii=False, robust=True):
    """Excitation x [CI] grid used by the 'sled' estimator (the timbre prior).
    Returns (templates, log prior weights). robust=True adds a weak-[CI] template (L' ratio 0.05,
    prior 0.2) so that one faint off-ladder line cannot make the wrong octave look certain;
    the first version (robust=False: [CI] ratios 0.15/0.45 only) was confidently wrong on
    HerBS-108-like data in up to a third of noise draws (see traps.py and docs/development-notes.md)."""
    rcis = ((0.05, 0.2), (0.15, 0.4), (0.45, 0.4)) if robust else ((0.15, 0.5), (0.45, 0.5))
    T, lw = [], []
    for Jc in (2.0, 2.8, 3.7, 5.0, 7.0):
        for rci, prc in rcis:
            fciis = (8.0, 40.0) if include_cii else (20.0,)
            for fcii in fciis:
                T.append(sled_flux(Jc, rci, f_cii=fcii, include_cii=include_cii))
                lw.append(np.log(prc / 5 / len(fciis)))
    return np.array(T), np.array(lw)


def _nearest_observed(cover, reach):
    """For every pixel, the nearest observed pixel within `reach` pixels (itself when observed), or -1."""
    n = len(cover)
    idx = np.arange(n)
    near = np.where(cover, idx, -1)
    for k in range(1, reach + 1):
        for cand in (idx - k, idx + k):
            ok = (near < 0) & (cand >= 0) & (cand < n)
            ok[ok] = cover[cand[ok]]
            near[ok] = cand[ok]
    return near


class Engine:
    """Harmonic-model likelihood over a redshift grid aligned to the pixel grid."""

    def __init__(self, grid, cover, templates, include_cii=False, zmin=0.5, zmax=8.0, sigma=None, logw=None,
                 dv_max=DV_MAX):
        """sigma: optional relative noise per pixel (flux units). Template weights become the
        expected S/N of each line, w_i = f_i / sigma(nu_i(z)), so lines in noisier windows count
        for less - required when windows come from different receivers or instruments.
        dv_max: velocity offset (km/s) allowed for each line around its predicted position (0: none)."""
        self.dv_max = float(dv_max)
        self.nslide = int(round(self.dv_max / PIX_KMS))                  # the same, in pixels
        names, nus = line_list(include_cii)
        self.names, self.nus = names, nus
        self.lnz = np.arange(np.log(1 + zmin), np.log(1 + zmax), DLN)
        self.cover = np.asarray(cover, bool)
        self.sigma = None if sigma is None else np.asarray(sigma, float)
        P = np.rint((np.log(nus)[None, :] - self.lnz[:, None] - grid.lnmin) / DLN).astype(int)
        self._inr = (P >= 0) & (P < grid.npix)
        self._P = np.where(self._inr, P, 0)                              # predicted pixel of each line at each z
        self.W = np.atleast_2d(templates)
        nt = self.W.shape[0]
        self.logw = np.full(nt, -np.log(nt)) if logw is None else np.asarray(logw, float)
        # line pairs close enough for their matched-filter outputs to overlap (e.g. CO(7-6) / [CI](2-1),
        # 998 km/s apart): the pixel lag between two lines is the same at every z on the log grid
        lnu = np.log(nus)
        self.pairs = [(i, j, (lnu[i] - lnu[j]) / DLN) for i in range(len(nus)) for j in range(i + 1, len(nus))
                      if abs(lnu[i] - lnu[j]) / DLN < 200]
        self._geo = {}
        g = self._geometry(self.nslide)
        self.P, self.valid, self.inv_sig, self.rigid = g['P'], g['valid'], g['inv_sig'], g['rigid']
        self.nlines_cov = self.valid.sum(1)

    def _geometry(self, nslide):
        """Which lines count as covered at each z, and where their noise and filter efficiency are read, when
        every line may move by nslide pixels.

        A line counts as covered when an observed pixel lies within nslide pixels of its predicted position (at
        least 1, so that the one-pixel gaps between adjacent spectral windows do not hide a line). Its noise
        and filter efficiency are taken at the nearest observed pixel. Requiring the predicted pixel itself to
        be observed would let a reading shift by a pixel or two, within the velocity tolerance, and park a
        line it predicts but the data lack inside such a gap, escaping the penalty."""
        if nslide not in self._geo:
            near = _nearest_observed(self.cover, max(nslide, 1))
            Pn = near[self._P]
            valid = self._inr & (Pn >= 0)
            Pn = np.where(valid, Pn, 0)
            inv_sig = valid.astype(float) if self.sigma is None else valid / self.sigma[Pn]
            # the two lines of a close pair keep their predicted separation wherever both are covered (no
            # separate velocity offsets): otherwise one broad line could pull both predicted positions onto itself
            rigid = np.zeros(valid.shape, bool)
            for i, j, _ in self.pairs:
                both = valid[:, i] & valid[:, j]
                rigid[:, i] |= both; rigid[:, j] |= both
            self._geo[nslide] = dict(P=Pn, Pc=np.where(valid, self._P, 0), valid=valid, inv_sig=inv_sig,
                                     rigid=rigid, scale={}, inv={})
        return self._geo[nslide]

    def _nslide(self, dv_max):
        return self.nslide if dv_max is None else int(round(float(dv_max) / PIX_KMS))

    def scale(self, spec, nslide=None):
        """Expected matched-filter S/N per unit template flux for each (z, line): 1/sigma at the line, times
        the filter efficiency there (below 1 when the line sits near a window edge or gap)."""
        g = self._geometry(self.nslide if nslide is None else nslide)
        if spec not in g['scale']:
            eff = filter_efficiency(self.cover, spec)
            g['scale'][spec] = g['inv_sig'] * np.where(g['valid'], eff[g['P']], 0.0)
        return g['scale'][spec]

    def _inv_norm(self, spec, nslide=None):
        """1 / sqrt(w^T C w) for one profile, with w the expected S/N of each predicted line (template flux x
        scale) and C the overlap of the filter outputs at the predicted lines. Without the overlap term, a broad
        line seen by two close predicted lines is counted twice. (C uses the untruncated profile; near an edge
        the true overlap differs slightly, which matters only for the CO(7-6)/[CI](2-1) pair.)"""
        g = self._geometry(self.nslide if nslide is None else nslide)
        if spec not in g['inv']:
            s = self.scale(spec, nslide)
            n2 = (s ** 2) @ (self.W ** 2).T
            for i, j, lag in self.pairs:
                c = kernel_overlap(spec, lag)
                if c > 1e-4:
                    n2 = n2 + 2 * c * np.outer(s[:, i] * s[:, j], self.W[:, i] * self.W[:, j])
            n = np.sqrt(np.maximum(n2, 0.0))
            g['inv'][spec] = np.where(n > 0, 1.0 / np.where(n > 0, n, 1.0), 0.0)
        return g['inv'][spec]

    def filter(self, d, bank=BANK):
        """Matched-filter an S/N spectrum with this engine's coverage (see matched_filter)."""
        return matched_filter(d, bank, self.cover)

    def line_outputs(self, y, nslide=None):
        """Filter output of one profile at every predicted line, nz x nlines: the largest value among the observed
        pixels within nslide pixels of the predicted position (the value at the predicted position, or at the
        observed pixel next to it, when nslide is 0). The two lines of a close pair (CO(7-6) and [CI](2-1))
        stay at their predicted positions when both are covered, so their separation is fixed."""
        nslide = self.nslide if nslide is None else nslide
        g = self._geometry(nslide)
        out = y[g['P']]
        if nslide > 0:
            ym = maximum_filter1d(np.where(self.cover, y, -np.inf), 2 * nslide + 1, mode='constant', cval=-np.inf)
            out = np.where(g['rigid'], out, ym[g['Pc']])
        return np.where(g['valid'], out, 0.0)

    def loglike(self, Ybank, bank=BANK, dv_max=None):
        """Ybank: matched-filter outputs from matched_filter(d, bank, cover) with this engine's coverage (or
        self.filter(d, bank)), one row per profile in `bank` (same order).
        dv_max: velocity offset allowed per line, km/s (default: the engine's, DV_MAX). With 0 every line sits
        exactly at its predicted position; posterior(refine=...) uses that to place z within a solution."""
        nt = self.W.shape[0]
        nslide = self._nslide(dv_max)
        assert len(bank) == len(Ybank), 'Ybank rows must match the profiles in bank'
        L = np.empty((len(self.lnz), len(Ybank) * nt))
        for b, y in enumerate(Ybank):
            G = self.line_outputs(y, nslide) * self.scale(bank[b], nslide)
            T = np.maximum(G @ self.W.T, 0.0) * self._inv_norm(bank[b], nslide)
            L[:, b * nt:(b + 1) * nt] = 0.5 * T * T + self.logw
        m = L.max(axis=1)
        return m + np.log(np.exp(L - m[:, None]).sum(axis=1)) - np.log(len(Ybank))

    def log_prior(self, z_phot=None, sig_ln=0.20, prior=None):
        """Log prior on the ln(1+z) grid. prior: None (uniform in z), or a function of z giving a density in z."""
        lp = self.lnz.copy()                                # uniform in z => density ~ (1+z) on ln(1+z) grid
        if prior is not None:
            lp = lp + np.log(np.maximum(prior(np.exp(self.lnz) - 1), 1e-300))
        if z_phot is not None:
            lp = lp - 0.5 * ((self.lnz - np.log(1 + z_phot)) / sig_ln) ** 2
        return lp

    def locate(self, i, exact):
        """Index of the best redshift within the solution that contains index i. With per-line velocity offsets
        the likelihood is flat over about +/-nslide pixels around a solution; `exact` (log-likelihood x prior
        with dv_max = 0) picks the redshift inside it at which the lines sit where predicted."""
        if exact is None or self.nslide == 0:
            return int(i)
        lo, hi = max(i - 2 * self.nslide, 0), min(i + 2 * self.nslide + 1, len(self.lnz))
        return int(lo + np.argmax(exact[lo:hi]))

    def posterior(self, loglike, z_phot=None, sig_ln=0.20, prior=None, refine=None):
        """prior: None (uniform in z), or a function of z giving a probability density in z.
        refine: optional loglike(..., dv_max=0) of the same data. The solution is chosen with `loglike`; the
        redshift returned is then the best one inside that solution according to `refine` (see locate). Without
        it the redshift can be off by up to dv_max, which is well inside the tolerance for a correct answer."""
        lp = self.log_prior(z_phot, sig_ln, prior)
        lpost = loglike + lp
        p = np.exp(lpost - lpost.max()); p /= p.sum()
        i = self.locate(int(np.argmax(p)), None if refine is None else refine + lp)
        win = np.abs(self.lnz - self.lnz[i]) < CONF_WIN
        return np.exp(self.lnz[i]) - 1, float(p[win].sum()), p


def prior_spt(z):
    """SPT-like DSFG redshift density (Reuter+20: 1.9 < z < 6.9, median 3.9)."""
    return np.exp(-0.5 * ((z - 3.9) / 1.1) ** 2)


def prior_herschel(z):
    """Bright-Herschel-like density (z-GAL: 0.8 < z < 6.55, median 2.56): lognormal in 1+z."""
    return np.exp(-0.5 * (np.log((1 + z) / 3.56) / 0.25) ** 2) / (1 + z)


def detect(y, cover, thresh=5.0, sep=12):
    """Peaks of a matched-filter S/N trace inside the coverage."""
    yy = np.where(cover, y, -np.inf)
    cand = np.where((yy > thresh) & (yy >= np.roll(yy, 1)) & (yy >= np.roll(yy, -1)))[0]
    cand = cand[np.argsort(-yy[cand])]
    keep = []
    for c in cand:
        if all(abs(c - k) > sep for k in keep):
            keep.append(c)
    out = []
    for c in sorted(keep):
        if 0 < c < len(y) - 1:          # parabolic refinement
            a, b, cc = y[c - 1], y[c], y[c + 1]
            den = a - 2 * b + cc
            off = 0.5 * (a - cc) / den if den != 0 else 0.0
        else:
            off = 0.0
        out.append((c + off, float(y[c])))
    return out


def naive_spacing(dets, grid):
    """'Pitch = partial spacing': assume all detections are adjacent CO rungs."""
    if len(dets) < 2:
        return None
    nu = np.sort([np.exp(grid.lnmin + p * DLN) for p, _ in dets])
    f0 = np.min(np.diff(nu))
    J = np.maximum(1, np.rint(nu / f0))
    return float(np.median(J * CO[0] / nu) - 1)


def _rule_solutions(dets, grid, cover, y, ztol_kms=300.0, zmin=0.5, zmax=8.0, ci_prunes=True):
    """ln(1+z) of every solution that explains all detections without predicting strong lines (expected
    S/N > 8 from the median SLED) where the data show < 3 sigma. Shared by sudoku_rule and guided_rule.

    ci_prunes=False keeps a solution whose only missing strong line is [CI]: the [CI] to CO ratio varies by a
    large factor from galaxy to galaxy, so an observer would not rule a redshift out on a missing [CI] line."""
    names, nus = line_list()
    tmpl = sled_flux(3.7, 0.29)
    lnd = np.array([grid.lnmin + p * DLN for p, _ in dets]); snr = np.array([s for _, s in dets])
    tol = ztol_kms / C_KMS
    cands = sorted({round(float(np.log(n) - l), 4) for l in lnd for n in nus})
    sols = []
    for lz in cands:
        if not (np.log(1 + zmin) <= lz <= np.log(1 + zmax)):
            continue
        lnpred = np.log(nus) - lz
        match = [np.argmin(np.abs(lnpred - l)) for l in lnd]
        if any(abs(lnpred[m] - l) > tol for m, l in zip(match, lnd)):
            continue
        scale = np.max(snr / tmpl[match])
        ok = True
        for i, lp in enumerate(lnpred):
            if i in match or (not ci_prunes and names[i].startswith('[CI]')):
                continue
            p = int(round((lp - grid.lnmin) / DLN))
            if 0 <= p < grid.npix and cover[p] and scale * tmpl[i] > 8 and y[p] < 3:
                ok = False; break
        if ok and not any(abs(lz - s) < tol for s in sols):
            sols.append(lz)
    return sols


def sudoku_rule(dets, grid, cover, y, ztol_kms=300.0, zmin=0.5, zmax=8.0):
    """Simplified High-z Sudoku: unique solution explaining all detections without predicting
    strong lines (expected S/N > 8 from the median SLED) where the data show < 3 sigma.
    Returns (z or None, number of surviving solutions)."""
    if not dets:
        return None, 0
    sols = _rule_solutions(dets, grid, cover, y, ztol_kms, zmin, zmax)
    if len(sols) == 1:
        return float(np.exp(sols[0]) - 1), 1
    return None, len(sols)


def guided_rule(dets, grid, cover, y, confirm=3.0, ztol_kms=300.0, zmin=0.5, zmax=8.0, prior=None, share=0.95,
                ci_prunes=True):
    """The detection rule followed by a guided look for weaker lines, as an observer would do by eye.

    If the 5 sigma detections leave one solution, that is the answer (as sudoku_rule). If they leave several,
    a solution counts as confirmed when the data reach `confirm` sigma within DV_MAX of a line it predicts
    inside the coverage and that no detection accounts for; the rule answers when exactly one solution is
    confirmed.

    prior: optional redshift density (a function of z), so that the rule can be given the same information as
    the engine. When the lines do not single out a solution, the rule then also answers if one of the remaining
    solutions (the confirmed ones, or all surviving ones when none is confirmed) holds at least `share` of the
    prior weight among them. The weight is the density per unit ln(1+z), as in Engine.log_prior.
    ci_prunes: see _rule_solutions.
    Returns (z or None, number of surviving solutions, number confirmed)."""
    if not dets:
        return None, 0, 0
    sols = _rule_solutions(dets, grid, cover, y, ztol_kms, zmin, zmax, ci_prunes)
    if len(sols) == 1:
        return float(np.exp(sols[0]) - 1), 1, 1
    nus = line_list()[1]
    lnd = np.array([grid.lnmin + p * DLN for p, _ in dets])
    tol, n = ztol_kms / C_KMS, int(round(DV_MAX / PIX_KMS))
    yc = np.where(cover, y, -np.inf)
    confirmed = []
    for lz in sols:
        for lp in np.log(nus) - lz:
            if np.min(np.abs(lp - lnd)) <= tol:
                continue                                   # this line is one of the detections
            p = int(round((lp - grid.lnmin) / DLN))
            if 0 <= p < grid.npix and np.max(yc[max(p - n, 0):p + n + 1]) >= confirm:
                confirmed.append(lz); break
    if len(confirmed) == 1:
        return float(np.exp(confirmed[0]) - 1), len(sols), 1
    pool = confirmed if confirmed else sols
    if prior is not None and pool:
        z = np.exp(np.array(pool)) - 1
        w = np.asarray(prior(z), float) * (1 + z)
        if w.sum() > 0 and w.max() >= share * w.sum():
            return float(z[int(np.argmax(w))]), len(sols), len(confirmed)
    return None, len(sols), len(confirmed)


# Line-ratio template for the cross-correlation: CO excitation and [CI] strength fitted to the composite spectrum
# of SPT galaxies in Spilker et al. (2014, Table 2), where L' in 1e9 K km/s pc2 is 296.6 for CO(1-0), 179.1 for
# CO(4-3), 136.3 for CO(6-5) and 47.5 for [CI](1-0): exp(-3/Jc) = 0.60 and exp(-5/Jc) = 0.46 give Jc = 6.0-6.4.
XCORR_TEMPLATE = (6.2, 0.16)


def _xcorr_setup(grid, weights, fwhm_kms, zmin, zmax):
    nus = line_list()[1]
    W = sled_flux(*XCORR_TEMPLATE) if weights is None else np.asarray(weights, float)
    lnz = np.arange(np.log(1 + zmin), np.log(1 + zmax), DLN)
    P = np.rint((np.log(nus)[None, :] - lnz[:, None] - grid.lnmin) / DLN).astype(int)
    inr = (P >= 0) & (P < grid.npix)
    return lnz, np.where(inr, P, 0), inr, W / W.max(), kernel(fwhm_kms)


def cross_correlation(d, grid, cover, weights=None, fwhm_kms=DETECT_FWHM, zmin=0.5, zmax=8.0):
    """S/N of a template cross-correlation at each trial redshift, after Yun et al. (2015, appendix A).

    zeta(z) = sum_i W_i s(p_i(z)): the spectrum `d` (in units of its noise), smoothed with the line profile,
    read at the pixel where redshift z puts line i of the list, and weighted by the relative strength W_i of
    that line. Lines outside the coverage contribute nothing. This is a re-implementation from the paper's
    description, not the authors' code, and it differs from it in three stated ways:
      - the paper takes the weights from a composite spectrum of SPT galaxies and gives no numbers; here they
        are this package's line-ratio template fitted to that composite (XCORR_TEMPLATE);
      - the template lines have the detection profile's width; the paper gives no width;
      - the spectrum is divided by its noise first, which helps on scans whose noise differs between windows.
    The paper takes the noise of zeta at each z from 10,000 random shuffles of the spectrum. This function
    returns the limit of that procedure in closed form, the mean and variance of a weighted sum under a random
    permutation of the observed pixels, so it needs no random numbers; shuffled_cross_correlation does the
    shuffles and the tests compare the two.
    Returns (ln(1+z) grid, S/N); the S/N is 0 where no line of the list falls in the coverage."""
    cover = np.asarray(cover, bool)
    lnz, P, inr, W, k = _xcorr_setup(grid, weights, fwhm_kms, zmin, zmax)
    x = d[cover]
    n, m, v = x.size, x.mean(), x.var()
    c = cover.astype(float)
    s, k1, k2 = (fftconvolve(a, b, mode='same') for a, b in ((np.where(cover, d, 0.0), k), (c, k), (c, k * k)))
    zeta = (np.where(inr, s[P], 0.0) * W).sum(1)
    a1 = (np.where(inr, k1[P], 0.0) * W).sum(1)
    a2 = (np.where(inr, k2[P], 0.0) * W ** 2).sum(1)
    var = n / (n - 1.0) * v * np.clip(a2 - a1 ** 2 / n, 0.0, None)
    ok = var > 1e-9
    return lnz, np.where(ok, (zeta - m * a1) / np.sqrt(np.where(ok, var, 1.0)), 0.0)


def shuffled_cross_correlation(d, grid, cover, n_shuffle, rng, weights=None, fwhm_kms=DETECT_FWHM, zmin=0.5,
                               zmax=8.0):
    """cross_correlation with the noise of zeta taken from `n_shuffle` random shuffles of the observed pixels,
    as the paper does. Slow; used to check the closed form."""
    cover = np.asarray(cover, bool)
    lnz, P, inr, W, k = _xcorr_setup(grid, weights, fwhm_kms, zmin, zmax)
    zeta = lambda arr: (np.where(inr, fftconvolve(arr, k, mode='same')[P], 0.0) * W).sum(1)
    z0, x = zeta(np.where(cover, d, 0.0)), d[cover]
    acc, acc2, buf = np.zeros(len(lnz)), np.zeros(len(lnz)), np.zeros(len(d))
    for _ in range(n_shuffle):
        buf[cover] = rng.permutation(x)
        zz = zeta(buf)
        acc += zz; acc2 += zz * zz
    mean = acc / n_shuffle
    std = np.sqrt(np.clip(acc2 / n_shuffle - mean ** 2, 0.0, None))
    ok = std > 1e-5
    return lnz, np.where(ok, (z0 - mean) / np.where(ok, std, 1.0), 0.0)


def xcorr_redshift(d, grid, cover, prior=None, **kw):
    """Best redshift of the template cross-correlation: the trial redshift with the highest S/N, or with the
    highest S/N x p(z) when a redshift density `prior` is given (the paper multiplies the S/N by the probability
    distribution of a photometric redshift). Returns (z, S/N of the cross-correlation there)."""
    lnz, snr = cross_correlation(d, grid, cover, **kw)
    score = snr if prior is None else np.clip(snr, 0.0, None) * np.asarray(prior(np.exp(lnz) - 1), float)
    i = int(np.argmax(score))
    return float(np.exp(lnz[i]) - 1), float(snr[i])


def classify(z_est, z_true, include_cii=False):
    """Label a redshift error by the line swap that explains it."""
    if z_est is None:
        return 'none'
    d = np.log(1 + z_est) - np.log(1 + z_true)
    if abs(d) < TOL_LN:
        return 'correct'
    names, nus = line_list(include_cii)
    best, lab = 1.0, 'other'
    for i in range(len(nus)):
        for j in range(len(nus)):
            if i == j:
                continue
            off = np.log(nus[j] / nus[i])      # true line i read as line j
            if abs(d - off) < TOL_LN and abs(d - off) < best:
                best = abs(d - off)
                if i < NCO and j < NCO:
                    from math import gcd
                    a, b = j + 1, i + 1; g = gcd(a, b); a //= g; b //= g
                    lab = f'{a}/{b}'
                else:
                    lab = 'CO<->[CI]'
    return lab
