# Pre-publication review

Reviewed 2 November 2026. This is an automated code and methods review, not external scientific peer review.

**Recommendation: address the four high-priority findings below before publishing the project with its current headline claims.** The existing analysis is reproducible, but reproducibility alone does not establish correctness. The window-order and comparator checks change real-data results without adding observations.

The implementation, input spectra, existing documentation, and canonical results were left unchanged. This report, its [evidence](pre-publication-evidence.json) and the [diagnostic script](../scripts/review_prepublication.py) are additions.

## Verification performed

| Check | Outcome |
|---|---|
| Requested dependency | Installed Astropy 8.0.2; FITS import and extractor CLI work. Astropy was already declared in `requirements.txt` and the `archive` extra. |
| Standard reproduction | `bash reproduce.sh` completed in a temporary copy, including all nine existing tests, real-data runs, stress tests, and the configured figure regeneration. The expensive `--full` simulations were not rerun. |
| Numerical agreement | Batch-result JSON/CSV files reproduced exactly. SPT0027's JSON differed in 12 floating-point values, with maximum absolute difference 5.72e-15. Eight regenerated figures differed under the local plotting version. |
| Package build | A normal isolated `pip wheel --no-deps .` build succeeded. An initial build without isolation failed because the host lacked `wheel`; that is not a project build failure. |
| FITS extraction | A generated descending-frequency cube produced the expected 12 × 17 × 17 cutout, with recovered noise 1.001 mJy/beam for injected 1 mJy/beam noise. |
| Input validation | All 52,360 cube-route CSV rows and 41,127 ADMIT-route rows had finite numeric fields, positive frequencies, and positive noise. |
| Release-content scan | No matches for the scanned private-key, GitHub-token, AWS-key, credential-assignment, or private-home-path patterns. This was a bounded pattern scan of current files, not a security guarantee or history scan. |

Environment: Python 3.11.1, NumPy 2.2.4, SciPy 1.17.1, Matplotlib 3.10.8, Astropy 8.0.2, pytest 8.3.5. Installing Astropy also installed pyerfa and astropy-iers-data and upgraded packaging from 24.2 to 26.3.

Reproduce the additional findings from the checkout root:

```bash
python scripts/review_prepublication.py
```

The command reports diagnostics; its successful exit is not a statement that the findings are fixed. The saved evidence was produced with `--output docs/pre-publication-evidence.json`.

## Findings

### P1 Overlapping windows make inference depend on input order

Locations: `ladder_redshifts/batch_real.py:85–95` and `ladder_redshifts/realspec.py:148–160`.

Each window assigns directly into `d[g]` and `sig[g]`, replacing measurements already present in overlapping channels. `read_spectrum` preserves first-appearance order of windows. Reordering the same observations therefore selects different data for inference.

| Source | Original window order: z, confidence | Reversed order: z, confidence |
|---|---|---|
| SPT2349-52 | 2.67512, 0.64877 | 3.90024, 0.90101 |
| SPT2152-40 | 3.85145, 0.99994 | 3.85145, 0.91514 |
| SPT2311-45 | 2.50685, 0.35099 | 2.50744, 0.45388 |

SPT2152-40 changes from confident to ambiguous. This means the seven-of-seven headline is not invariant to ordering of the input windows.

**Required change:** adopt an explicit, order-independent overlap policy. Independent measurements can be combined in flux units with propagated inverse-variance weights; duplicate or correlated products require identification rather than independent coaddition. Apply the same policy in both gridders. Add a real-spectrum permutation test, then regenerate the real results and follow-up inputs. The example does not establish that the reversed ordering is the right scientific answer.

### P1 The advertised equal-weight comparator is implemented with unequal weights

Locations: `ladder_redshifts/realspec.py:177`, `ladder_redshifts/ladder.py:188–190`, and `docs/method.md:52–54`.

The real-data comparator passes an all-ones flux template **and** `sigma=sig_e`. The engine converts that into weights proportional to `1/sigma_i` on already standardized line S/N. This is an equal-flux template, not the documented equal-S/N sum.

Changing only this comparator construction to unit weights on S/N, while retaining the existing filters, overlap correction, redshift prior, and confidence calculation, changes SPT0436-40 from a wrong z=2.6385 at confidence 0.905 to the published solution z=3.85145 at confidence 0.9999997. The confident-correct count becomes **four of seven**, rather than three. This diagnostic is not a complete implementation of Z-Spec.

The attribution also needs correction. [Lupu et al. 2012, section III.1](https://arxiv.org/html/1009.5983v3#S3.SS1) uses multiple test statistics, including a total-flux/noise statistic and a pairwise S/N statistic, with additional selection and significance tests. Neither the current comparator nor the unit-S/N diagnostic reproduces that complete procedure.

**Required change:** choose and name the intended baseline precisely, implement that definition consistently, and regenerate the comparison table and figure. Describe the existing simplified Sudoku rule as a local baseline, as the README already largely does. Reserve claims about the published methods for a faithful comparison.

### P1 Synthetic observations include information outside their coverage

Locations: `ladder_redshifts/atlas.py:72–77`, `ladder_redshifts/traps.py:95–99`, `ladder_redshifts/tuning.py:83–100`, and `ladder_redshifts/ambiguity.py:63–66`.

These paths filter a spectrum across the full grid before applying each observation's coverage. `Engine` restricts predicted line centres to covered pixels, but the convolution at those centres can still include signal and noise outside the observed windows. The tuning optimizer caches these unrestricted filter outputs for every candidate layout.

The diagnostic creates two spectra that are identical at every covered pixel, changes one pixel immediately outside coverage, and obtains a maximum log-likelihood change of **1063.17**. The large injected value is a deliberate invariance test, not an estimate of the bias in the shipped population. The effect on aggregate published simulation results remains unquantified.

**Required change:** make filtering and template normalization depend on the actual observed pixels. Merely zeroing excluded pixels is incomplete: the filter norm and correlations must also account for truncation at edges and gaps. Use one observation model for simulations and real spectra, add an invariance test for excluded pixels, and rerun affected synthetic studies.

### P1 The general calibrated-confidence claim contradicts the saved simulation

Locations: `README.md:8–9`, `CITATION.cff:7`, and `results/redshift_atlas.json` under `calibration / B3 single tuning / sled`.

For the 101 synthetic single-tuning Band 3 cases with confidence at least 0.99, the stored mean confidence is **0.9998825**, while only **94/101 = 0.9306931** are correct. All seven errors also lie outside the wider ±0.006 confidence window, so this discrepancy is not caused by the different correctness and confidence tolerances. These values were checked against the per-source raw archive, not just the summary JSON.

The small real sample and the seven confident successes cannot establish a general calibration claim. Moreover, the likelihood maximizes over amplitude before averaging templates and widths. Normalizing that profile-likelihood score over redshift does not by itself supply a marginalized Bayesian posterior or calibrated probability.

**Required change:** describe the current quantity as a model-based confidence score or normalized profile-likelihood mass, state the demonstrated calibration limitations, and avoid interpreting 0.95 as a verified 95% success probability. If a Bayesian posterior is intended, specify and integrate an amplitude prior. Validate or calibrate on held-out simulations and independent observations after the coverage and gridding fixes. Preserve both correctness and confidence-window metrics explicitly.

### P2 Extraction failures can be reported as a successful run

Locations: `scripts/alma/extract_cube_spectra.py:171–174`, `205–217`.

For a FITS file with a normal pipeline-style filename but no frequency axis, `process_cube` raises an error. The wrapper catches it, logs `FAILED`, prints `DONE`, and exits with status **0**, producing no spectra. Automation can therefore continue after an incomplete extraction, potentially using older saved spectra through the batch fallback.

**Required change:** count and summarize failures, return a nonzero status when required inputs fail, and distinguish a deliberately partial run from a complete one. Exercise both loose-FITS and archive-member paths.

### P2 A one-width filter bank crashes the public inference function

Location: `ladder_redshifts/realspec.py:182–184`.

`run_engine(..., bank=(450.0,))` passes the width bank to the likelihood successfully, but detection then indexes `Y[1]` unconditionally and raises `IndexError`. A reordered multi-width bank also silently changes the detector's width.

**Required change:** select a documented detection width by value, or validate the required bank shape explicitly. Keep the default behavior intact and add a one-width/reordered-width regression.

## Questions for scientific reviewers

- **Noise model and resampling.** The covariance correction treats filter overlap, but not interpolation-induced channel covariance. A 300-draw white-noise diagnostic finds matched-filter noise standard deviations of 0.96, 1.44, and 1.88 for native 50, 100, and 200 km/s channels resampled by the current gridder. The latter two are outside the approximately 50 km/s channel assumption; accepting arbitrary CSVs needs resolution checks or covariance propagation. Edge/gap truncation needs the same treatment even for nominally matched resolution.
- **Model completeness.** Evaluate velocity offsets between species, non-Gaussian or multiple line components, correlated noise, continuum residuals, multiple redshift components, atypical excitation, and redshifts outside the search range. The present tests mostly reproduce cases used during development.
- **Follow-up probabilities.** The planner uses ten draws per candidate/condition, discards candidates below 1% and renormalizes the remainder, and sets success probability to zero without Monte Carlo draws when a noiseless case has truth mass below 0.5. These are approximations, not rigorous zero-probability bounds. The tuning optimizer also selects and evaluates layouts on the same simulated population. Ask for uncertainty intervals, sensitivity analysis, and evaluation on an independent population before treating a chosen tuning as optimal.
- **Data extraction.** An ALMA specialist should assess point-source flux losses, beam handling, noise annuli, continuum removal, duplicate products, and the ADMIT zero-fill approximation. The synthetic FITS smoke test does not validate the original reductions.

## Publication and provenance

The directory contains no `.git` metadata, so there is no recoverable commit history, remote, or historical secret scan available here. Current files total about 11.7 MiB, with no large-file obstacle. The code/documentation license split and citation metadata exist. The general ALMA acknowledgement matches the current [ALMA instructions](https://almascience.eso.org/alma-data/publication-acknowledgement); those instructions additionally specify an NRAO/Green Bank acknowledgement for North American authors.

The [CDS catalogue description](https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/ApJ/902/78?format=html&tex=true) confirms the published table and its redshift/line-identification fields. Attempts to retrieve its machine-readable rows failed at CDS and the Harvard mirror. A complete independent row-by-row check of `truth.csv` and target positions remains outstanding; the existing tests use the shipped truth file. The Band 2 citation was checked against [Bakx 2026](https://arxiv.org/html/2604.01089v2).

Before the public release, record a reproducible environment, run the installed package in CI as well as tests from the checkout, and add tests for the failures above. The current CI tests Python 3.10 and 3.12, while the package advertises Python ≥3.9; either exercise the supported lower bound or revise that declaration based on testing. The normal reproduction command does not regenerate the trap figure, despite promising every figure; `--full` does.

After corrections, rerun the full simulations, regenerate all affected results and prose, verify the external truth table, and create an identifiable Git commit. External reviewers should receive the methods, data provenance, reproduction instructions, and any unresolved limitations together.
