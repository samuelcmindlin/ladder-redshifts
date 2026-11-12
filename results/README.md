# Results files

Every file here is written by a module of `ladder_redshifts` or, for one, by a script. `reproduce.sh` regenerates
the fast ones, and `reproduce.sh --full` regenerates all of them. [`docs/results.md`](../docs/results.md) has the tables drawn from these
files.

| File | Contents | Written by |
|---|---|---|
| `redshift_real_batch.json`, `.csv` | Real spectra under `data/spectra/` (2016.1.00672.S cube cut-outs and SPT0027-50) scored against `data/truth.csv`. For each source: engine z and confidence; `p_published` (posterior mass on the published z); verdict; the equal-weight-sum, detection-rule and naive-spacing answers; detections; candidates. | `batch_real` |
| `redshift_real_batch_admit.json`, `.csv` | The same for the ADMIT-route spectra, plus the result from the measured windows alone. | `batch_real --admit` |
| `redshift_real_SPT0027.json` | SPT0027-50 (Band 2): detections, best z and confidence for each estimator and prior, candidates, and the added-noise and reduced-coverage stress tests. | `realspec` |
| `redshift_followup_SPT0314-44.json` | Follow-up planner for SPT0314-44: candidates, the lines each one predicts, and P(confident & correct) for every Band 3–7 LO under three conditions. `unseparated` marks a candidate that a tuning cannot tell from a rival, because it adds only lines both predict. | `followup` |
| `redshift_atlas.json` | Synthetic atlas: confident-correct and confidently-wrong fractions for each setup and estimator, by S/N and by redshift; error taxonomy; calibration tables; `reliability` (how often answers above 0.95 and 0.99 are right). | `atlas` |
| `redshift_atlas_raw.npz` | Per-source raw output behind the atlas (true z, S/N, estimates, confidences); `atlas --plot-only` re-plots from it. | `atlas` |
| `redshift_traps.json` | The three rebuilt traps, for each variant and step, with the adopted (robust) and the first (narrow) [CI] prior. | `traps` |
| `redshift_traps_examples.npz` | Mean posterior per step for the figure; `traps --plot-only` re-plots from it. | `traps` |
| `redshift_tuning.json` | Tuning optimiser: greedy tuning positions and outcomes for Band 2 now, Band 2 with the upgrade, and Band 3, under two time models and two source brightnesses. Each entry has the fraction on the selection population and, under `held_out`, on an independent population. | `tuning` |
| `redshift_ambiguity.json` | HerBS-108's undertone family of solutions, the anharmonic mismatches, and the sonification mapping. | `ambiguity` |
| `redshift_noise_only.json` | The engine on spectra that contain only noise: how often the confidence exceeds 0.5, 0.8 and 0.95, with the per-line velocity tolerance and with exact line positions. | `scripts/noise_only_check.py` |
