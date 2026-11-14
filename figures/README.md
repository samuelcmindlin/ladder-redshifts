# Figures

| Figure | What it shows | Made by |
|---|---|---|
| `summary_real.png` | The headline result. For each of the 18 sources of 2016.1.00672.S with a published z, the confidence of this engine and of the equal-weight sum, whether each best z is the published one, and whether the detection rule finds a unique solution. | `summary` |
| `redshift_real_batch.png` | One row per real source (SPT0027-50 and the 21 sources of 2016.1.00672.S): the spectrum with the engine's lines at its best z, and the posterior mass per candidate (green = published z). | `batch_real` |
| `redshift_real_batch_admit.png` | The same for the ADMIT-route pseudo-spectra. Grey bands are windows covered where ADMIT gave no spectrum on the source (zero-filled). | `batch_real --admit` |
| `redshift_real_SPT0027.png` | SPT0027-50, Band 2 science verification, in four panels (below). | `realspec` |
| `redshift_followup_SPT0314-44.png` | Which single follow-up tuning settles SPT0314-44. It plots P(confident & correct) against LO for Bands 3–7, with the lines each hypothesis predicts. | `followup` |
| `redshift_ambiguity.png` | Redshift as pitch, in three panels (below). | `ambiguity` |
| `redshift_atlas.png` | Synthetic atlas of 2,500 sources × 6 setups, in six panels (below). | `atlas` |
| `redshift_traps.png` | HerBS-108, SPT0551-50 and SPT0452-50 rebuilt from their published coverage. It shows the posterior on each candidate before and after the decisive data. | `traps` |
| `redshift_tuning.png` | Tuning optimiser. Confident, correct redshifts against the number of tunings for Band 2 now and after the upgrade, at fixed depth and at fixed total time. The layouts are chosen on one simulated population and scored on an independent one. | `tuning` |
| `redshift_degradation.png` | The seven real multi-line spectra with noise added, means over the sources against the S/N left to the second line. Panel A, answers given: the engine with no prior and with the survey prior, the best of 72 rule settings (the one with the most correct answers among those wrong no more often than the engine), the single setting used in the first comparison, and the rule that needs every line at 5σ. Panel B, best guesses: the engine and the template cross-correlation, each with and without the survey prior. | `degrade` |
| `redshift_degradation_sources.png` | The same test per source, at 4.5σ, 4σ, 3.5σ and 3σ: the engine (no prior, survey prior), the best of 72 rule settings and the 5σ rule, with the count of sources for which the engine is ahead of and behind the best rule. | `degrade` |
| `redshift_operating_points.png` | Answers against wrong answers over all 5,600 noise draws. Each rule setting is a point; the engine is a curve along its confidence thresholds. Panel A: no prior (72 settings). Panel B: each source's 2020 photometric redshift as the prior (216 settings). Settings with more than 30 wrong answers are off the plot. | `degrade` (from `redshift_degradation.json` and `redshift_rule_family.csv`) |
| `redshift_photoz.png` | The 11 sources with one line in the 3 mm scan, each with its published photometric redshift as the prior. Panel A: the published redshift, the other redshifts the scan allows (engine with no prior, p ≥ 0.02), and the photometric redshift ± error from the 2020 paper and from the 2024 revision. Panel B: the engine's probability for the published redshift under the survey prior and under the three versions of the photometric prior, on a logit axis, with the 0.95 threshold. | `photoz` |
| `redshift_photoz_noise.png` | The noise-added test with each source's 2020 photometric redshift as the prior. Panel A: fraction of draws in which the published redshift is returned, for the engine with that prior and with none, for the best of 216 rule settings with that prior, and for the rule as first built. Panel B: fraction in which a different redshift is returned. | `photoz` (from `redshift_degradation.json`) |

## Panels

**`redshift_real_SPT0027.png`**

- **A.** The collapsed image with the aperture.
- **B.** The spectrum, with the engine's lines.
- **C.** The posterior when only part of the coverage is kept.
- **D.** The spectrum with extra noise added.

**`redshift_ambiguity.png`**

- **A.** HerBS-108's two lines and their undertone family of candidate fundamentals.
- **B.** The ladder's ambiguity function, which has side lobes at 2, ½, 3/2 and 2/3.
- **C.** The anharmonic tie-breaker.

**`redshift_atlas.png`**

- **A.** The setups.
- **B.** Confident-and-correct fractions for each estimator.
- **C.** Confidently-wrong fractions for each estimator.
- **D.** Errors fall on harmonic ratios.
- **E.** Calibration.
- **F.** Redshift deserts.

## `independent/`: the 23 sources of the independent test

| Figure | What it shows | Made by |
|---|---|---|
| `independent/summary.png` | One row per source, grouped by the number of lines published from the 3 mm scan: the engine's probability for the published redshift with the survey prior and with the source's 2020 photometric redshift, on a logit axis, and the answers of the two rule settings chosen beforehand. | `scripts/independent_report.py` |
| `independent/redshift_real_batch.png` | The 23 spectra with the engine's lines at its best redshift, as `redshift_real_batch.png`. | `batch_real` |
| `independent/redshift_photoz.png` | The nine one-line sources with their photometric redshifts as priors, as `redshift_photoz.png`. | `photoz` |
