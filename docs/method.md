# Method

## 1. The ladder as a harmonic series

A rigid rotor emits at ν_J = 2B·J for the transition J → J−1. For CO, ν_J = 2B₀J − 4D₀J³ with B₀ = 57635.968 MHz
and D₀ = 0.1835 MHz, so ν_J ≈ J × 115.271 GHz.

- **The CO lines.** They form an overtone series, slightly compressed by centrifugal stretching: CO(7–6) is
  0.53 cent flat of 7 × CO(1–0). The engine uses catalogue rest frequencies for J = 1–13.
- **Off-ladder lines.** The \[CI\](1–0) (492.161 GHz) and \[CI\](2–1) (809.342 GHz) lines are off-ladder, like
  inharmonic partials. [CII] (1900.537 GHz) is optional.
- **Redshift.** Redshift divides every frequency by 1+z, so on ln ν it is a rigid shift. The observed ladder is an
  overtone series on the fundamental F = 115.271/(1+z) GHz, and that fundamental is never observed.

## 2. The ambiguity is an undertone family

Two observed lines with ratio ν_b/ν_a = b/a (in lowest terms) are explained by the assignments (J_a, J_b) = (k·a, k·b)
for every k = 1, 2, 3, …:

    candidate fundamentals   F_k = ν_a / (k·a) = F₁ / k     (an undertone series)
    candidate redshifts      1 + z_k = k · (1 + z₁)         (an overtone series in 1+z)

This is the subharmonic matching of virtual-pitch theory (Terhardt), and the gcd structure behind High-z Sudoku.
For example, HerBS-108 has lines at 97.75 and 146.60 GHz, a ratio of 3/2. The candidates are z = 1.359, 3.717 and
6.073 (k = 1, 2, 3), and the true value is k = 2.

What breaks the tie:

1. **Missing rungs.** For k ≥ 2, rungs k·a+1, …, k·b−1 fall between the observed lines. If they are covered and
   absent, k is disfavoured.
2. **Off-ladder lines** ([CI], [CII]) at their own ratios.
3. **The spectral line energy distribution (SLED).** The relative strengths of the rungs act as a timbre.
4. **Anharmonicity** (section 6). It works in principle but needs very high S/N.

## 3. Likelihood

The spectrum is put on a log-frequency grid with constant 50 km/s pixels. Every line then has the same width in
pixels, and a redshift is a pure pixel shift.

**Observed pixels only.** Only pixels inside the frequency coverage enter the matched filter, and its output is
normalised so that it has unit noise everywhere, including near window edges and gaps:

    y(p) = Σⱼ kⱼ d(p+j) c(p+j) / √(Σⱼ kⱼ² c(p+j))

Here d is the S/N spectrum, k the unit-norm line profile and c the coverage mask (1 observed, 0 not). A line
centred at p then has its expected output multiplied by the *filter efficiency* e(p) = √(Σⱼ kⱼ² c(p+j)). This is 1
inside a window, about 0.7 at an edge, and 0 outside. Data outside the coverage cannot affect the result, and real
and simulated spectra go through the same observation model.

Let yᵢ(z) be this output at the predicted position of line i, for the lines whose centres are covered.

**The statistic.** Take a template with expected S/N wᵢ = eᵢ fᵢ/σᵢ (filter efficiency × expected flux over the
local noise) and an unknown amplitude A ≥ 0, in white noise. Maximising over A gives the generalized log-likelihood
ratio ln Λ(z) = ½ T², with

    T = max(0, Σᵢ wᵢ yᵢ) / √(wᵀ C w)

- **The overlap matrix C.** C_ij is the correlation between the filter outputs at lines i and j: the kernel
  autocorrelation at the pixel lag between them, with C_ii = 1. On the log grid that lag is the same at every z.
  C is nearly diagonal. The exception is CO(7–6) and \[CI\](2–1), 998 km/s apart, where 800 km/s filters correlate at
  0.115. Without the off-diagonal term, one broad line seen from both predicted positions is counted twice.
  Near a window edge the true overlap differs slightly from the untruncated value used here. This matters only
  for that one pair.
- **The missing-rung penalty.** Empty predicted rungs add to the denominator but not to the numerator, so the
  penalty is automatic. This is the "undertone check" of pitch detection, in likelihood form.
- **Velocity offsets between lines.** Each line may sit up to 150 km/s (3 pixels) from its predicted position:
  yᵢ is the largest filter output among the observed pixels within that range, so the statistic is also
  maximised over one velocity offset per line.
  - **Why.** The lines of one galaxy differ by tens of km/s (kinematics, differential lensing), and the ladder's
    anharmonicity separates octave-family readings of the same lines by 29–122 km/s (section 6). With exact
    positions, a bright spectrum was "decided" by such offsets: 25–76 km/s was enough to turn a correct answer at
    confidence 1.00 into a wrong one at 1.00 ([`development-notes.md`](development-notes.md), item 21).
  - **The redshift quoted.** The likelihood is now flat over about ±150 km/s around a solution. The solution is
    chosen with the tolerance; the redshift quoted for it is the one inside that range at which the lines sit
    where predicted (the same likelihood without offsets).
  - **The close pair.** CO(7–6) and \[CI\](2–1) keep their predicted separation wherever both are covered.
    Otherwise one broad line could pull both predicted positions onto itself.
  - **The cost.** Under pure noise the largest of seven neighbouring filter outputs is positive on average, so
    absent lines are penalised slightly less. In 600 pure-noise Band 3 scans (SPT-like prior) the confidence
    exceeded 0.5 in 4.8% of cases, against 2.2% with exact positions; it never reached 0.95 either way. For
    Band 3 + Band 4 scans the figures are 11.8% against 4.5%, with 3 and 2 of the 600 reaching 0.95
    (`results/redshift_noise_only.json`). Section 5 of [`results.md`](results.md) gives the cost on simulated
    sources.
- **Lines next to gaps.** A line counts as covered when an observed pixel lies within that tolerance of its
  predicted position (within one pixel when no tolerance is used), and its noise and filter efficiency are read
  at the nearest observed pixel. Adjacent spectral windows often leave one-pixel gaps. If the predicted pixel
  itself had to be observed, a reading could shift by a pixel or two and hide a line it predicts, but the data
  lack, inside such a gap.
- **Equal weights (the comparison baseline).** Setting wᵢ = 1 on every covered line gives an equal-weight sum of
  line S/N normalised by √N (with C). This is the analogue of Hermes' subharmonic summation, and it is close in
  spirit to the line-sum statistics of Lupu et al. (2012). Their Z-Spec procedure uses several statistics (a
  total-flux/noise sum and a pairwise S/N sum) and further significance tests, so this baseline is a simplified
  stand-in, not their method. It ignores both the expected line ratios and the noise differences between windows.
  The repository uses it as the "equal-weight sum" comparison.
- **Template cross-correlation (the closest earlier method).** Yun, Heyer & Aretxaga (2007) and Yun et al. (2015,
  appendix A) find redshifts from Redshift Search Receiver spectra by cross-correlating the spectrum with a
  template of many transitions, each weighted by its relative strength in a composite spectrum of SPT galaxies.
  They set the significance by shuffling the spectrum and multiply the result by a photometric-redshift
  probability. For one template, T above is the same weighted sum, normalised analytically by the noise and the
  coverage in place of the shuffles. What this engine adds is the bank of templates with a prior on excitation
  and \[CI\] strength, the conversion to a posterior with a confidence score and a list of competing redshifts,
  and the velocity tolerance. A re-implementation from the paper's description is compared with the engine
  on noise-added real spectra in [`results.md`](results.md), section 9; the authors' own code has not been run.

**Templates (the timbre prior).**

- **CO.** CO brightness-temperature ratios are r_J1 = exp(−(J−1)/J_c), with J_c ∈ {2.0, 2.8, 3.7, 5.0, 7.0}.
  J_c ≈ 3.7 approximates the average for submillimetre galaxies.
- **\[CI\](1–0).** The \[CI\](1–0)/CO(1–0) L′ ratio takes the values 0.05, 0.15 and 0.45, with prior weights 0.2, 0.4
  and 0.4. The weak-[CI] template stops one faint off-ladder line from making the wrong octave look certain.
- **\[CI\](2–1).** The \[CI\](2–1)/\[CI\](1–0) L′ ratio is fixed at 0.45.
- **Line widths.** The filter bank has 250, 450, 800 and 1200 km/s FWHM. The 1200 km/s profile was added after
  the pre-publication review. Without it, a single line broader than about 1000 km/s was better fit by a CO(7–6) +
  \[CI\](2–1) pair; see [`development-notes.md`](development-notes.md).

The likelihood ratio Λ is averaged over the 15 templates (with their prior weights) and the 4 widths.

**The redshift distribution and confidence.**

- **What is normalised.** P(z) ∝ π(z) Σ_θ π(θ) Λ_θ(z) on z = 0.5–8, normalised over z.
- **Redshift priors.** Uniform in z; SPT-like (normal, mean 3.9, s.d. 1.1, from Reuter et al. 2020); or
  Herschel-like (lognormal in 1+z, modelled on the z-GAL sample, median z ≈ 2.6). A source's own photometric
  redshift can be used in their place, as a normal distribution in z (section 10).
- **Not a full Bayesian posterior.** Λ is a *profile* likelihood (maximised over the line amplitude), not a
  marginal likelihood with an amplitude prior. The code calls P(z) the posterior for brevity.
- **Confidence.** The mass of P(z) within |Δ ln(1+z)| < 0.006 (1800 km/s) of its maximum. The window is wide
  enough that blend partners count together.
- **How to read it.** The confidence is a model-based score, not a verified probability. In simulated
  multi-tuning scans, answers above 0.95 are right in 99.4–99.9% of cases, and for a single Band 3 tuning in about
  98%, while the stated values average 0.998–1.000. On the real single-line sources the stated values are lower
  than the hit rate. Both directions are measured in [`results.md`](results.md), section 5. Between readings that
  predict the same line positions the score rests on line ratios alone and is not reliable for bright lines
  ([`development-notes.md`](development-notes.md), "Known open issues").
- **Confident and correct.** An answer is *confident* if its confidence is ≥ 0.95. It is *correct* if it lies within
  0.003 (900 km/s) of the published redshift.
- **Candidates.** The posterior is also summarised as a list of candidate redshifts with their mass.

## 4. Comparison estimators

All share the same matched-filter front end (`ladder.py`).

| Name | What it does |
|---|---|
| `sled` | The likelihood above (this work). |
| `flat` | The same likelihood with equal weights on the S/N of every covered line: no line ratios and no noise weighting. A simplified stand-in for line-sum methods such as Lupu et al. (2012), not their procedure. |
| `sudoku` | A simplified High-z Sudoku rule. Lines are detected at ≥ 5σ. A solution is kept if it explains every detection to within 300 km/s and predicts no line of expected S/N > 8 where the data show < 3σ. The rule answers only when a single solution survives. This is a re-implementation for comparison, not the authors' code. |
| `naive` | Assumes the detections are adjacent CO rungs and takes the smallest spacing as the fundamental. It fails whenever [CI] sits between rungs, and is kept as a baseline. |
| `guided` | The `sudoku` rule followed by a guided look for a weaker line. If the 5σ detections leave several solutions, a solution is confirmed when the data reach 3σ (or 4σ) within 150 km/s of a line it predicts and that no detection accounts for. The rule answers when exactly one solution is confirmed. Given a redshift prior, it also answers when one of the remaining solutions holds at least 0.95 of the prior weight among them. It stands for what an observer does by eye and is used in the degradation and photometric-prior tests (section 10). With `ci_prunes=False` a missing \[CI\] line does not count against a solution. |
| `joint` | A joint likelihood of spectrum and photometry, re-implemented after Bing et al. (2024, section 4.2) in `scripts/compare_joint_likelihood.py`, outside the package. At each trial redshift the flux of every CO and \[CI\] line is predicted from the infrared luminosity, and the model spectrum is compared with the whole spectrum by chi-square. In the form used for the comparison, each line's strength is integrated over the scatter of its luminosity relation. It is given the same prior as the engine and scored with the engine's confidence. Section 10 has the details. |
| `xcorr` | A template cross-correlation re-implemented after Yun et al. (2015, appendix A): ζ(z) is the sum over all lines of the smoothed spectrum at each line's predicted position, weighted by the line's relative strength, and its S/N is ζ less its mean, over its standard deviation, under random shuffles of the spectrum. The best guess is the redshift of highest S/N, or of highest S/N × p(z) when a prior is given. It has no confidence and never abstains. Used only in the degradation test (section 10). |

## 5. Ambiguity function

Plot the noiseless E[ln Λ(z′)] / E[ln Λ(z)] against (1+z′)/(1+z). It peaks at ratios J′/J. The tallest side lobes
are at 2 and 1/2, because every rung J ≤ 6 has an octave partner. The next are at 2/3, 3/2, 3 and 4/3 (panel B of
`figures/redshift_ambiguity.png`). This is the log-frequency autocorrelation of the set {ln J}, the same reason
octave errors dominate pitch tracking.

## 6. Precision and the anharmonic tie-breaker

**Cramér–Rao bound.** For a Gaussian line of width σ_v at matched-filter S/N s, the centroid bound is
σ(v) = √2·σ_v/s. For N lines, 1/σ²_{ln(1+z)} = Σᵢ sᵢ² / (2σ_v² c⁻²). For example, FWHM 450 km/s at S/N 30 gives
9.0 km/s per line and 12.7 km/s on a pair difference.

**Anharmonicity.** Suppose a pair is read at the wrong k. Because the ladder is not exactly harmonic, the two lines
then imply redshifts that differ by c·|ln(ν_{ka}/ν_{kb}) − ln(ν_a/ν_b)|:

| Pair → alternative | Mismatch (km/s) |
|---|---|
| (2,3) ↔ (4,6) | 28.6 |
| (1,3) ↔ (2,6) | 45.8 |
| (3,4) ↔ (6,8) | 40.1 |
| (4,5) ↔ (8,10) | 51.6 |
| (2,3) ↔ (6,9) | 76.4 |
| (1,3) ↔ (3,9) | 122.2 |

Octave confusions need S/N ≳ 30 in both lines to be broken this way. Factor-3 confusions become testable at
S/N ≳ 10–15, but systematic velocity offsets between transitions of tens of km/s (from kinematics or differential
lensing) would compete. The engine does not use this effect: the velocity tolerance of section 3 (150 km/s) is
wider than every mismatch in the table. Until that tolerance was added it did use it, at high S/N
([`development-notes.md`](development-notes.md), item 21).

## 7. Real spectra

- **Cube cut-outs.** `scripts/alma/extract_cube_spectra.py` cuts a box around the source in each spectral window,
  averages channels to about 50 km/s, and measures the per-channel noise in an annulus. It writes one `.npz` file per
  window.
- **Point-source spectrum.** `realspec.spectrum_point` takes the mean over pixels within 0.75″ of the published
  position, with noise from the annulus. It is used for cut-outs centred on a known position (2016.1.00672.S).
- **Aperture spectrum.** `realspec.aperture` keeps the connected region above 3σ in the channel-averaged image
  that touches the central 3″, grown by half a beam. It is used for SPT0027-50, a lensed source with two images.
- **Continuum.** Each window's continuum is removed with a sigma-clipped linear fit.
- **Noise scale.** The annulus lies 1.3–2.6 box radii from the source in a primary-beam-corrected cube, where
  the noise is higher than at the source. The annulus noise is therefore too large: by a factor of about 1.3
  for the 8″ boxes of 2016.1.00672.S (1.2 at 84 GHz, rising to 1.4 at 115 GHz) and about 1.7 for the 12″ boxes
  of SPT0027-50. `realspec.calibrate_noise` corrects the scale from the spectrum itself:
  - **Per window.** The scatter of flux/noise over the line-free channels (3σ-clipped) should be 1.
  - **Per source.** One window has about 120 channels, so its own value is uncertain by about 6%. The values of
    all windows are fitted with a straight line in frequency, and each window's noise is multiplied by the
    fitted value.
  - **Record.** The factors are in `data/noise_calibration.csv`. The channel-to-channel shape of the noise still
    comes from the annulus.
  - **Cause.** The primary-beam explanation fits the dependence on box size and on frequency, but it has not
    been checked against the cubes, which are not shipped. The correction does not depend on it.
- **Onto the grid.** `realspec.grid_windows`, called by `batch_real.to_grid_auto` and `realspec.to_grid`, puts
  each window's flux and noise on the log-frequency grid.
  - **Overlaps.** Where windows overlap (adjacent tunings do), the measurements are combined in flux units with
    inverse-variance weights, as independent measurements, so the result does not depend on the order of the
    windows.
  - **Gaps.** Gaps inside a window are not bridged.
  - **Noise weights.** The combined noise sets the template weights wᵢ = eᵢfᵢ/σᵢ, so a line predicted in a noisy
    window counts for less.
- **Channel width.** The grid assumes channels of about 50 km/s. Wider channels make neighbouring pixels correlated
  after interpolation. Runs of channels wider than 62.5 km/s therefore have their noise inflated by
  √(width/50 km/s), which keeps the matched-filter noise close to 1 (an approximation, with a warning). Channels
  much narrower than 50 km/s are sampled rather than averaged, so they should be binned first. All the shipped
  spectra have channels of about 50 km/s.
- **ADMIT pseudo-spectra.** Built from the archive's ADMIT summaries (`scripts/alma/admit_spectra.py`), these mark
  each window as *measured* (a spectrum on the source exists) or *zero-filled* (flux 0, noise = the cube's
  per-channel noise × 10). See `data/README.md`.

## 8. Follow-up planner

`followup.py` takes a real, ambiguous spectrum and, for every local-oscillator (LO) setting in ALMA Bands 3–7, works
out how likely one more tuning is to settle the redshift:

- **The tuning.** One tuning is two 3.75 GHz sidebands (four 1.875 GHz windows), centred in the receiver's
  intermediate-frequency range: 4.125–7.875 GHz either side of the LO in Bands 3, 4, 5 and 7 (IF 4–8 GHz) and
  5.375–9.125 GHz in Band 6 (IF 4.5–10 GHz). Every LO on the grid that keeps both sidebands inside the band is
  tried. Real set-ups can place the windows elsewhere in the IF, so each LO here stands for one common choice.
- **Simulation.** Each candidate redshift is taken in turn as the truth. Its predicted lines (low or typical
  excitation) are added, with noise, to the real data, and the engine is re-run. Line fluxes are scaled from
  the flux of the observed line: its matched-filter S/N times the noise where it sits.
- **Shared lines do not count.** A tuning can separate two candidates only if it covers a line that one of them
  predicts and the other does not (no line of the other within 300 km/s). A line both predict, such as CO(9–8)
  at z = 2.93 and CO(6–5) at z = 1.62, shows that one of them is right, not which. When a tuning leaves the true
  candidate unseparated from a rival, its success probability for that candidate is set to 0.
- **Score.** The result is the posterior-weighted probability of a confident and correct answer.
- **Approximations.** Each candidate and condition gets 10 noise draws. Candidates with less than 1% of the
  posterior are dropped and the rest renormalised. When even the noise-free follow-up leaves less than 0.5
  of the posterior on the true redshift, the draws are skipped and the success probability is set to 0. So the
  numbers are planning estimates, with Monte Carlo
  uncertainties of about ±0.1 per LO, not rigorous bounds.

## 9. Synthetic studies

`atlas.py` simulates 2,500 sources × 6 ALMA setups:

- **Redshift and brightness.** z ~ U(1, 7). S/N is log-uniform from 2 to 64, quoted as the S/N of the CO(4–3) line.
- **Excitation.** J_c is lognormal with median 3.7, and each line scatters by 0.1 dex.
- **[CI].** The \[CI\](1–0) L′ ratio is lognormal with median 0.29 and 0.18 dex scatter.
- **Width and photo-z.** FWHM is lognormal with median 450 km/s, and the photometric-redshift scatter is 0.20 in
  ln(1+z).

These parameters are illustrative. The atlas uses a random seed first used after the last change to the
likelihood (12 November 2026), so its sources played no part in shaping the method.

The simulated sources come from the same family of models that the engine fits: the same ladder shape with other
parameter values and 0.1 dex of scatter, Gaussian lines at their predicted positions, and independent Gaussian
noise of known size. The synthetic error rates therefore describe the engine when its model holds. They do not
test what happens when it fails (other species, line shapes, velocity offsets, baseline ripples); only real
spectra do that.

`traps.py` rebuilds three published octave-type misidentifications from their
frequency coverage. `tuning.py` chooses tuning positions greedily for ALMA Band 2, now and after the Wideband
Sensitivity Upgrade. It chooses them on one simulated population and reports their performance on an independent
one.

## 10. Two checks against other ways of working

**The High-z Sudoku code.** `scripts/compare_sudoku.py` runs the published code of Bakx & Dannerbauer (2022),
unmodified, from a copy the user supplies. Its function `RSGquality` classes a redshift by the lines that fall
inside a set of frequency windows; it reads no spectrum.

- **Real sources.** Each is classed at its published redshift with the windows of its own spectrum.
- **Atlas.** The six setups are classed on the redshifts of the atlas's sources with CO(4–3) at S/N ≥ 16.
- **A detail of the call.** `RSGquality` builds its redshift grid up to 1.5 times the largest redshift it is given
  and reads neighbouring solutions beyond the grid at its last point. Every call therefore includes a sentinel
  redshift of 12, whose own class is subtracted from the result. The grid has 40,000 points.

**Noise added to real spectra.** `degrade.py` takes the seven multi-line spectra of 2016.1.00672.S and adds
Gaussian noise, channel by channel, until the second-brightest line would have a chosen matched-filter S/N:
the noise is raised by the factor (observed S/N) / (chosen S/N), and the stated noise is raised with it. Every
method is then run on 100 draws per source and level. Each draw keeps the noise already in the real spectrum, so
the draws of one source are not independent.

- **Two kinds of outcome.** The engine (when confident), the detection rule and the guided rule give an answer or
  none; they are compared on the answers they give. The cross-correlation always has a best redshift; it is
  compared with the engine's best guess.
- **Two states of knowledge.** Everything is run with no redshift prior, where each method has only the spectrum,
  and again with the survey prior.

**A joint likelihood after Bing et al. (2024).** `scripts/compare_joint_likelihood.py` is written from section 4.2
of the paper; no code is published.

- **The model.** For a trial redshift z and infrared luminosity L_IR(z), each line of the paper's Table 7 has
  L′ = 10^((log L_IR − A) / N), which becomes a flux through the relation of Solomon et al. (1997) and a Gaussian
  of 500 km/s FWHM on the 50 km/s grid. The log-likelihood relative to no lines is Σ (m d − m²/2) / σ² over the
  observed pixels, with d the continuum-subtracted flux and m the model. This is the paper's exp(−χ²/2).
- **The luminosity.** At each redshift a modified blackbody (emissivity index 2, optically thick below 100 μm,
  temperature free) is fitted to the photometry of Reuter et al. (2020, Table D.1) at rest wavelengths of 50 μm
  or more, with 10% added to the errors, and integrated over 8–1000 μm. It is the apparent luminosity; the line
  fluxes are apparent too. The paper fits galaxy templates instead.
- **The scatter.** The paper fixes each line at its predicted strength. The variant used for the comparison
  integrates each line's strength factor s over a log-normal of the relation's scatter:
  ∫ p(s) exp(s B − s² A / 2) ds, with B = Σ m d / σ² and A = Σ m² / σ² for that line. For a line measured above
  4σ the integrand is far narrower than any fixed grid in log s, so the nodes are placed across the data's own
  Gaussian in s; otherwise across the log-normal. Against direct numerical integration the result is right to
  0.01 in the log, except where a line predicted above an S/N of 40 is absent, where the value is below −25.
- **The photometric term and the confidence.** The paper multiplies by the likelihood of its template fit and
  reports the maximum. Here the term is the same prior the engine gets, and the likelihood goes through
  `Engine.posterior`, so the confidence and its threshold are the engine's.
- **Sixteen variants,** and how one is chosen on the development sample and then fixed:
  [`results.md`](results.md), section 12.

**The cross-correlation, in detail.** `ladder.cross_correlation` is written from the description in Yun et al.
(2015), not from their code, and differs from it in three ways.

- **Weights.** The paper takes them from a composite spectrum of SPT galaxies and gives no numbers. Here they are
  this package's line-ratio template with J_c = 6.2 and a \[CI\](1–0) ratio of 0.16, fitted to the luminosities
  of CO(1–0), CO(4–3), CO(6–5) and \[CI\](1–0) in Table 2 of Spilker et al. (2014).
- **Line width.** The paper gives none. The template lines have the 450 km/s detection profile.
- **Noise.** The spectrum is divided by its noise before the cross-correlation, which helps on scans whose noise
  differs between windows. The paper estimates the noise of ζ(z) from 10,000 shuffles of the spectrum. The code
  uses the limit of that procedure in closed form: under a random permutation of N pixel values with mean m and
  variance v, a weighted sum Σ a_p x_p has mean m Σ a_p and variance N/(N−1) · v · (Σ a_p² − (Σ a_p)²/N). A test
  compares it with 1,500 real shuffles.
- **Prior.** The paper multiplies the S/N by the probability distribution of each source's photometric redshift,
  which it takes from a radio-to-millimetre spectral index. In the noise-added test the survey prior stands in
  for it; in the photometric-prior test below it is each source's dust-continuum photometric redshift.

**Each source's photometric redshift as the prior** (`photoz.py`; [`results.md`](results.md), section 10).

- **The prior.** A normal distribution in z, centred on the published photometric redshift, with the published
  error as its standard deviation. It replaces the survey prior; the two are not multiplied.
- **Three versions.** The values printed in Reuter et al. (2020, Table F.1); the collaboration's 2024 revision of
  that table; and the 2020 values recalibrated. The recalibration takes the pull, (z_spec − z_phot) / error,
  of the 60 SPT sources that are not in this sample: its mean is −0.479 and its standard deviation 1.101. The
  recalibrated centre is z_phot − 0.479 × error and the width 1.101 × error
  (`scripts/photoz_calibration.py`, `results/redshift_photoz_calibration.json`).
- **The rules with a photometric redshift.**
  - **Guided rule.** As above, with the same normal distribution: when the lines do not single out a solution,
    it answers if one remaining solution holds 0.95 of the prior weight.
  - **Window rule.** The guided rule with a top-hat prior of half-width 0.13 (1 + z_phot), which is the default
    photometric accuracy (`dzUncertainty`) of the High-z Sudoku code. It answers when exactly one solution lies
    in the window. It is not that code's procedure, which classes coverage at a known redshift.
  - **Two forms.** The detection rule drops a solution that predicts a line of expected S/N > 8 where the data
    show < 3σ, and by default that includes \[CI\]. With `ci_prunes=False` only CO lines count. The
    \[CI\]-to-CO ratio varies by a large factor between galaxies, so the second form is closer to what an
    observer does. Both are reported.
- **The cross-correlation** is multiplied by the same normal distribution, as in the paper.
- **In the noise-added test** the three versions are run as further priors for the engine, the guided rule (both
  forms) and the cross-correlation. The likelihood is computed once per draw and each prior applied to it
  (`realspec.engine_likelihood`, `engine_posterior`), so the results under the earlier priors are unchanged.

**A family of rule settings, and the engine's threshold** (`degrade.py`; [`results.md`](results.md), sections 9
and 10).

- **Why.** The guided rule has thresholds, and a comparison with one setting of them mostly measures that
  setting. The first comparison did exactly that and overstated the engine.
- **The family.** Every combination of four dials, 72 settings:
  - a line counts as detected from 5, 4.5 or 4σ;
  - a weaker line confirms a solution from 4, 3.5, 3 or 2.5σ;
  - a missing \[CI\] line counts against a redshift, or not;
  - a missing line counts against a redshift when its expected S/N is above 5, 8 or 12.

  With a prior, each is run with the share of prior weight at which one solution is taken set to 0.95, 0.90 or
  0.80: 216 settings.
- **Which setting is quoted.** Fixed before the family was run: over all draws, the setting with the most
  correct answers among those that return a wrong redshift no more often than the engine does with the same
  prior. A second choice allows up to 1% of draws wrong. The best setting at each noise level separately is
  given as an upper limit; it is chosen with knowledge of the outcome.
- **The engine at other thresholds.** Its answers are counted at confidence 0.5, 0.6, 0.7, 0.8, 0.9, 0.95 and
  0.99. For each, the rule setting with the most correct answers among those with no more wrong answers is
  found, so the two are compared at equal numbers of wrong answers.
- **A bias to keep in mind.** The quoted rule setting is picked by its score on the same draws it is then scored
  on. The engine's 0.95 was fixed beforehand. This favours the rule slightly.
- **Bound on the wrong-answer rate.** k wrong answers among n confident ones give an exact (Clopper-Pearson)
  one-sided 95% upper bound; for k = 0 it is 1 − 0.05^(1/n).
