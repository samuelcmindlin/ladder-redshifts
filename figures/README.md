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
| `redshift_degradation.png` | The seven real multi-line spectra with noise added, mean over the sources, against the S/N the second-brightest line is left with. Panel A, answers given: the engine without a prior and with the survey prior, the guided rule and the 5σ detection rule. Panel B, best guesses: the engine and the template cross-correlation, each without and with the prior. | `degrade` |
| `redshift_degradation_sources.png` | The same test per source, with the second line at 4.5σ, 4σ, 3.5σ and 3σ: the engine without a prior and with the survey prior, the guided rule and the 5σ rule, one row per source. | `degrade` |

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
