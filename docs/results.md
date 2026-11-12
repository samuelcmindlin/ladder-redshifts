# Results

Every number here is produced by the code in `ladder_redshifts/` and stored under `results/`. Definitions are in
[`method.md`](method.md). An answer is **confident** if its confidence score is at least 0.95: that share of the
normalised likelihood × prior (called the posterior below, for brevity) lies within ±0.006 in ln(1+z) of the best
redshift. It is **correct** if it lies within 0.003 of the published value. The score is not a calibrated
probability; section 5 measures how reliable it is.

## 1. Real spectra: ALMA project 2016.1.00672.S, from pipeline-cube cut-outs

The sample is 21 SPT sources, each observed in 5 tunings × 4 spectral windows covering 84.1–114.8 GHz.

- **Extraction.** Cut-outs of ±20″ around the published positions; a point-source spectrum at the centre, with
  its noise put on the scale of its own line-free channels ([`method.md`](method.md), section 7).
- **Truth.** Reuter et al. 2020, Table E1.
- **Prior.** SPT-like redshift prior.

Per-source output is in `results/redshift_real_batch.csv`. There, `p_published` is the posterior mass on the
published z. "Flat sum" is the equal-weight sum of line S/N on the same spectra. "Rule" is the simplified Sudoku-style detection
rule, not the authors' code. Run with `python -m ladder_redshifts.batch_real`.

| Source | Lines in the published 3 mm scan | Engine z (confidence) | Published z | P(published) | Flat sum | Rule | Verdict |
|---|---|---|---|---|---|---|---|
| SPT0403-58 | CO(4–3) CO(5–4) \[CI\](1–0) | 4.056 (1.00) | 4.0564 | 1.00 | 4.056 (1.00) | 4.055 | confident & correct |
| SPT0425-40 | CO(5–4) CO(6–5) | 5.135 (1.00) | 5.1353 | 1.00 | 5.135 (1.00) | 5.134 | confident & correct |
| SPT0436-40 | CO(4–3) \[CI\](1–0) | 3.851 (1.00) | 3.8519 | 1.00 | 3.851 (1.00) | 3.851 | confident & correct |
| SPT0520-53 | CO(4–3) \[CI\](1–0) | 3.778 (1.00) | 3.7785 | 1.00 | 4.972 (0.66) ✗ | 3.778 | confident & correct |
| SPT0652-55 | CO(4–3) \[CI\](1–0) | 3.346 (1.00) | 3.3466 | 1.00 | 2.260 (0.91) ✗ | 3.345 | confident & correct |
| SPT2152-40 | CO(4–3) \[CI\](1–0) | 3.851 (1.00) | 3.8507 | 1.00 | 2.639 (0.91) ✗ | 3.851 | confident & correct |
| SPT2203-41 | CO(5–4) CO(6–5) | 5.194 (1.00) | 5.1937 | 1.00 | 5.194 (1.00) | 5.194 | confident & correct |
| SPT2232-61 | CO(3–2) | 2.893 (0.90) | 2.8936 | 0.90 | 2.893 (0.90) | — | ambiguous, best = published |
| SPT0314-44 | CO(3–2) | 2.933 (0.90) | 2.9345 | 0.90 | 2.933 (0.90) | — | ambiguous, best = published |
| SPT0402-45 | CO(3–2) | 2.683 (0.87) | 2.6833 | 0.87 | 2.683 (0.91) | — | ambiguous, best = published |
| SPT2349-52 | CO(4–3) | 3.900 (0.85) | 3.9000 | 0.85 | 2.675 (0.91) ✗ | — | ambiguous, best = published |
| SPT0226-45 | CO(4–3) | 3.233 (0.78) | 3.233 | 0.78 | 3.233 (0.78) | — | ambiguous, best = published |
| SPT2129-57 | CO(4–3) | 3.260 (0.77) | 3.2601 | 0.77 | 3.260 (0.77) | — | ambiguous, best = published |
| SPT0528-53 | CO(5–4) | 4.736 (0.76) | 4.737 | 0.76 | 2.442 (0.49) ✗ | — | ambiguous, best = published |
| SPT0054-41 | CO(5–4) | 4.876 (0.69) | 4.8769 | 0.69 | 4.876 (0.69) | — | ambiguous, best = published |
| SPT2031-51 | CO(3–2) | 2.452 (0.45) | 2.4521 | 0.45 | 2.452 (0.50) | — | ambiguous, best = published |
| SPT2311-45 | CO(3–2) | 2.507 (0.42) | 2.5073 | 0.42 | 2.507 (0.56) | — | ambiguous, best = published |
| SPT2316-50 | CO(4–3) | 4.178 (0.11) | 3.1413 | 0.02 | 4.178 (0.14) ✗ | — | ambiguous, best ≠ published |

- **All 7 sources with two or three lines: confident and correct.** The detection rule also settles all 7 on the
  same spectra; the equal-weight sum settles 4. For the other three, the equal-weight sum prefers a wrong redshift
  (0.66–0.91 of its weight): it reads CO(4–3) as CO(3–2) or CO(5–4) and treats the \[CI\](1–0) line as noise.
- **All 11 single-line sources: ambiguous.** None is above 0.90, and the published z is on top in 10. All eleven
  redshifts are confirmed by follow-up in other bands: ten in Reuter et al., and SPT0314-44, which that paper
  lists as tentative, by Riechers (2025) (section 4).
- **No confidently wrong answer** in 18 sources, or 19 with SPT0027-50 (section 3).
- **The stated confidence of the ambiguous answers is cautious.** Across the 11 ambiguous answers the mean stated
  confidence is 0.68, while the top candidate is the published z in 10 of 11; the average weight on the published z
  is 0.67. Eleven sources are far too few to calibrate the score.
- **The faintest second lines.** After the matched filter, SPT2203-41's CO(6–5) is at 5.6σ, SPT0520-53's
  \[CI\](1–0) at 6.5σ, SPT2152-40's \[CI\](1–0) at 7.1σ and SPT0652-55's \[CI\](1–0) at 7.4σ (3–5σ per channel
  at the peak). All are above the detection rule's 5σ threshold, so these spectra do not test the engine below it.
- **A feature outside the line list.** SPT0402-45 has a second feature at 96.26 GHz (6.4σ). It matches HCN(4–3)
  at the published redshift, which the engine does not model.
- **Sources without a published redshift.** SPT0433-59, SPT0549-56 and SPT2334-50 show no line above 5σ. The
  engine's confidence is 0.52, 0.05 and 0.27. For SPT0433-59, three features of 2–3σ fall where z = 4.394 puts
  CO(4–3), \[CI\](1–0) and CO(5–4); pure noise does this in about 5% of simulated Band 3 scans
  (`results/redshift_noise_only.json`).
- **With a flat redshift prior** (`batch_real --prior uniform`): the 7 multi-line sources are still all confident and
  correct. The published z is on top for 8 of the 11 single-line sources, at a maximum confidence of 0.66.
- **Changes from earlier versions.**
  - **First release to the pre-publication review (2 November).** The review fixed how overlapping windows are
    combined, how the equal-weight baseline is weighted, and how window edges enter the filter
    ([`development-notes.md`](development-notes.md), items 13–16). SPT2349-52 then ranked its published
    redshift first (it was 2.675 before).
  - **Second review (12 November).** The noise of these spectra was about 1.3 times too high, and the likelihood
    now allows each line a velocity offset (items 20–22). The engine's verdicts did not change. The detection
    rule now settles all 7 multi-line sources; with the old noise it settled 4, because three second lines sat
    just below its 5σ threshold.

## 2. The same project from ADMIT summaries (first pass)

These are pseudo-spectra built from the archive's ADMIT summaries, without downloading cubes. Windows with no
on-source spectrum are zero-filled (noise × 10). Outputs: `results/redshift_real_batch_admit.*` and
`figures/redshift_real_batch_admit.png`; run with `python -m ladder_redshifts.batch_real --admit`.

| Source | Lines in the published 3 mm scan | ADMIT windows on the source | Engine z (confidence) | Published z | P(published) | Verdict |
|---|---|---|---|---|---|---|
| SPT0436-40 | CO(4–3), \[CI\](1–0) | 3 | 3.851 (1.00) | 3.8519 | 1.00 | confident & correct |
| SPT0314-44 | CO(3–2) | 2 | 2.934 (0.88) | 2.9345 | 0.88 | ambiguous, best = published |
| SPT2232-61 | CO(3–2) | 1 | 2.893 (0.76) | 2.8936 | 0.76 | ambiguous, best = published |
| SPT2129-57 | CO(4–3) | 1 | 3.260 (0.71) | 3.2601 | 0.71 | ambiguous, best = published |
| SPT0402-45 | CO(3–2) | 1 | 2.683 (0.64) | 2.6833 | 0.64 | ambiguous, best = published |
| SPT0054-41 | CO(5–4) | 3 | 4.875 (0.55) | 4.8769 | 0.55 | ambiguous, best = published |
| SPT0652-55 | CO(4–3), \[CI\](1–0) | 1 | 3.347 (0.51) | 3.3466 | 0.51 | ambiguous, best = published |
| SPT2152-40 | CO(4–3), \[CI\](1–0) | 1 | 3.851 (0.38) | 3.8507 | 0.38 | ambiguous, best = published |
| SPT2031-51 | CO(3–2) | 2 | 2.452 (0.26) | 2.4521 | 0.26 | ambiguous, best = published |
| SPT0403-58 | CO(4–3), CO(5–4), \[CI\](1–0) | 2 | 2.792 (0.62) | 4.0564 | 0.02 | ambiguous, best ≠ published |
| SPT0520-53 | CO(4–3), \[CI\](1–0) | 2 | 4.972 (0.52) | 3.7785 | 0.11 | ambiguous, best ≠ published |
| SPT2311-45 | CO(3–2) | 2 | 3.677 (0.30) | 2.5073 | 0.18 | ambiguous, best ≠ published |
| 9 others | — | 0 | — | — | — | no spectrum on the source |

- **No confidently wrong answer.**
- **The three misses come from ADMIT.** SPT0403-58's CO(5–4) and SPT0520-53's \[CI\](1–0) fall in windows where
  ADMIT's spectra lie elsewhere in the field. SPT2311-45 had a single line.
- **The zero-filled windows do most of the work.** With the measured windows alone, the published z is the top
  candidate only for SPT0436-40.
- **The cube cut-outs (section 1) turned 1 confident answer into 7.**
- **Noise.** The ADMIT-route spectra already had their noise scaled to the line-free scatter of each window, so the
  noise correction of section 1 does not apply to them.

## 3. SPT0027-50, ALMA Band 2 science verification (project 2011.0.00025.SV)

The data are 20 spectral windows covering 76.9–111.0 GHz, with a 0.6″ beam. The spectrum comes from an aperture of
264 pixels: the continuum above 3σ, grown by half a beam, covering both lensed images. Output is in
`results/redshift_real_SPT0027.json`; run with `python -m ladder_redshifts.realspec`.

**Detections** above 5σ (matched filter): 77.826 GHz (S/N 24.7), 103.760 GHz (35.1), 110.768 GHz (11.9).

**Estimators:**

| Estimator | z | Note |
|---|---|---|
| SLED likelihood | 3.4432 | confidence 1.000; published 3.444 |
| Equal-weight sum | 3.4432 | same as above |
| Detection rule | 3.4433 | unique solution |
| Naive spacing | 15.65 | wrong; fooled by [CI] |

**More noise** (30 draws per level):

| Noise | Confident & correct | Confidently wrong |
|---|---|---|
| ×1–4 | every time | — |
| ×6 | 93% | 0% |
| ×8 | 73% | 0% |

**Less coverage:**

| Coverage kept | Result |
|---|---|
| ≥ 84.2 GHz (Band 3-like) | correct, confidence 1.000 |
| CO(3–2) alone | ambiguous: 0.50 at the true z vs 0.34 at z = 5.32 |
| CO(4–3) alone | ambiguous, with the true z third at 0.25 |

The S/N values are about 1.7 times those of the first release, whose noise for this aperture spectrum was too high
by that factor ([`development-notes.md`](development-notes.md), item 20). The added-noise test starts from the
corrected noise, so each level is a fainter source than before.

## 4. Follow-up planner on a real single-line source: SPT0314-44

`python -m ladder_redshifts.followup` uses the cube cut-out spectrum. Here the line straddles the edge between two
tunings.

- **Detection.** One line at 87.922 GHz, matched-filter S/N 52.
- **Engine.** With the SPT-like prior: z = 2.9332 (CO(3–2)) at 0.897 and z = 1.6222 (CO(2–1)) at 0.103. That is
  ambiguous, with the best candidate equal to the published identification, z = 2.9345, which was tentative when
  published and has been confirmed since (end of this section). With a flat prior the split is 0.60 / 0.40.
- **Planner settings.** One tuning of 2 × 3.75 GHz in Bands 3–7, with the sidebands where each receiver's IF range
  puts them ([`method.md`](method.md), section 8), LO in 1 GHz steps, 10 noise draws. Candidates with at least 1%
  posterior are included, weighted by it.
- **Shared lines do not count.** The two candidates share every third harmonic of the detected line: 175.8, 263.7
  and 351.4 GHz (for example CO(9–8) at z = 2.933 and CO(6–5) at z = 1.622). A tuning that adds only such a
  line cannot settle the redshift and scores 0 ([`method.md`](method.md), section 8).

LO settings that reach P(confident & correct) ≥ 0.9:

| Condition | Band 4 (23 LOs) | Band 5 (38) | Band 6 (45) | Band 7 (83) | Total |
|---|---|---|---|---|---|
| Same noise per channel as the scan, low excitation | 11 | 12 | 7 | 19 | 49 |
| 3× noisier (≈ 1/9 of the time), low excitation | 11 | 12 | 7 | 10 | 40 |
| Same noise, typical excitation | 11 | 12 | 7 | 19 | 49 |

Band 3 is already covered, apart from 1 GHz at its top edge that holds no line of either candidate. Of the 190 LO
settings in all, 52 cover a line that only one candidate predicts. At the scan's noise level 49 of them settle the
redshift whichever candidate is true. The other three catch such a line only at the very edge of a window; they
settle it if z = 2.933 is right and not otherwise (P = 0.90). The remaining 138 score 0.

The best choice is Band 4 with the LO at 139 GHz. CO(3–2) at 131.9 GHz (if z = 1.622) falls in the lower sideband
and CO(5–4) at 146.5 GHz (if z = 2.933) in the upper one, so either outcome is a detection. LOs of 136–138 GHz
cover only the first of these and 140–142 GHz only the second; the planner still scores them at 1, because a line
that bright cannot be missed if it is there. All eleven Band 4 settings survive noise three times higher.

**What settled it in the end.** Riechers (2025, ApJ 980, 59) observed SPT0314-44 with the ACA in Band 6, with the
sidebands centred near 249 and 263 GHz, and found CO(9–8) together with H₂O emission and OH⁺ and NH absorption,
giving z = 2.9338 ± 0.0002. The engine's best candidate from the single Band 3 line, 2.9332, is 46 km/s from
that. The planner would not have picked this tuning: of those lines it knows only CO(9–8), which is the harmonic
the two candidates share, so it scores the LO settings near 256 GHz at 0. The redshift was decided by species
outside the engine's line list. For bright lensed sources the planner is therefore conservative: tunings that
work through H₂O or OH⁺ are invisible to it, and so is their evidence to the engine itself.

## 5. Synthetic atlas (2,500 sources × 6 setups)

`python -m ladder_redshifts.atlas`. S/N is the CO(4–3)-equivalent. The population is described in
[`method.md`](method.md), section 9. Its parameters are illustrative, and its random seed was first used after the
last change to the likelihood. The simulated lines have no velocity offsets between them. The table gives the
fraction with a confident, correct redshift at S/N 8–16, with the fraction confidently wrong in brackets.
"Equal-weight sum" and "Sudoku-style rule" are the simplified baselines of `method.md`, section 4.

| Setup | Naive spacing | Sudoku-style rule | Equal-weight sum | SLED likelihood | SLED + photo-z |
|---|---|---|---|---|---|
| Band 3 single tuning | 0.01 (3.6%) | 0.04 (0%) | 0.03 (0%) | 0.05 (0%) | 0.05 (0%) |
| Band 3 scan (84–115 GHz) | 0.16 (33%) | 0.48 (0%) | 0.46 (0%) | 0.51 (0%) | 0.59 (0.7%) |
| Band 4 scan (125–163 GHz) | 0.31 (30%) | 0.59 (0%) | 0.59 (0.2%) | 0.65 (0.2%) | 0.69 (0.2%) |
| Band 2 full, WSU (67–116 GHz) | 0.18 (47%) | 0.63 (0%) | 0.70 (2.2%) | 0.72 (0.2%) | 0.80 (0.4%) |
| Band 3 + Band 4 scans | 0.30 (55%) | 0.79 (0.4%) | 0.87 (5.8%) | 0.90 (0%) | 0.94 (0%) |
| Band 2 + Band 4 | 0.23 (64%) | 0.81 (0%) | 0.93 (2.4%) | 0.94 (0.2%) | 0.96 (0.2%) |

- **Low S/N.** At S/N 2–4 the SLED likelihood gives 3.3–5.8× more confident, correct redshifts than the detection
  rule in the five multi-tuning setups. At S/N 4–8 it gives 1.3–1.5× more. In every S/N bin it is confidently wrong
  in at most 0.8% of the scan sources (`results/redshift_atlas.json`, `by_rho`).
- **The equal-weight sum in wide setups.** It is confidently wrong in up to 2.7% (Band 2) and 6.4% (Band 3 + Band 4)
  of the sources in an S/N bin.
- **Naive spacing.** It fails because \[CI\] lines sit between CO rungs.

**How reliable is a confident answer?** This counts all S/N together (`reliability` in the JSON):

| Setup | SLED: answers ≥ 0.95 | of which correct | mean stated confidence | Equal-weight sum: of ≥ 0.95, correct |
|---|---|---|---|---|
| Band 3 single tuning | 93 | 97.8% (2 wrong) | 0.998 | 97.1% (70) |
| Band 3 scan | 1236 | 99.8% (3 wrong) | 0.999 | 99.4% (1136) |
| Band 4 scan | 1530 | 99.4% (9 wrong) | 0.999 | 99.4% (1436) |
| Band 2 full, WSU | 1750 | 99.7% (6 wrong) | 0.999 | 97.0% (1738) |
| Band 3 + Band 4 scans | 2122 | 99.6% (8 wrong) | 0.999 | 95.3% (2152) |
| Band 2 + Band 4 | 2208 | 99.9% (3 wrong) | 1.000 | 98.0% (2226) |

- **Scans.** A confident answer from the SLED likelihood is right in 99.4–99.9% of cases. That is below the stated
  mean of about 0.999, so the score is overconfident in the tail.
- **A single Band 3 tuning.** 2 of 93 confident answers are wrong: one CO line read as another rung, and one
  CO/\[CI\] confusion.
- **The price of the velocity tolerance.** The simulated lines have no velocity offsets, which is the case where
  exact line positions help most. On the same 2,500 sources, exact positions give 20 wrong answers among 8,931
  confident ones (all six setups together); the tolerance gives 31 among 8,939. Several of the wrong ones are
  off by a factor of 2 or 3 in the rung: readings that predict the same covered lines are separated only by line
  ratios and the prior ([`development-notes.md`](development-notes.md), "Known open issues").
- **Reading the score.** Treat it as a ranking with a conservative threshold in multi-tuning scans, not as a
  probability.
- **Earlier populations.** With seed 7, before the fixes of items 15–16, there were 7 wrong answers among 101
  single-tuning answers ≥ 0.99. With seed 20261102 and exact line positions (the release of 2 November), confident
  answers were right in 99.7–100% of cases in scans and in 101 of 104 for a single tuning.

## 6. Known octave traps, rebuilt from their published coverage

`python -m ladder_redshifts.traps`. Each source is rebuilt with 150 noise draws per variant (excitation or line
strength), using the adopted (robust) [CI] prior. Line strengths are illustrative where the papers do not give them.

| Trap | Before the decisive data | After |
|---|---|---|
| HerBS-108 (3 mm + 2 mm scans; z = 3.717 vs 1.358) | Weak-[CI] variants: flagged ambiguous in 49–95% of draws and never confidently wrong. With typical [CI] it is already settled. | + extra 3 mm setup (CO(3–2) at 73.3 GHz): confident & correct in 99–100% |
| SPT0551-50 (one Band 3 line; z = 3.164 vs 2.123) | Ambiguous in 97–100%, confidently wrong in at most 3.3%. Still ambiguous after adding the 1 mm window (96–100%; confidently wrong in at most 3.3%). | + [CII] window: confident & correct in 100% |
| SPT0452-50 (ATCA + Band 3 with CO(5–4) flagged; z = 5.016 vs 2.011) | Confident & correct in 0–43%, otherwise ambiguous; never confidently wrong | + empty [CI] windows: 83–91%; + re-reduced CO(5–4): 99–100%; + [CII]: 100% |

With the first, narrower [CI] prior (no weak-[CI] template), HerBS-108-like data are confidently wrong in up to a
quarter of draws; both versions are in `results/redshift_traps.json`.

The single faint line of SPT0551-50 (S/N 6) is where the velocity tolerance costs most: with exact line positions
the confidently wrong fraction was at most 0.7% before the 1 mm window and 2% after it.

## 7. Tuning optimiser for ALMA Band 2

`python -m ladder_redshifts.tuning`. The population has SPT-like redshifts.

- **How layouts are chosen.** Tunings are placed greedily on a 1 GHz grid, using one simulated population of 200
  sources.
- **How they are scored.** The chosen layouts are then scored on an independent population of 200 (a new seed), and
  those held-out fractions are the ones below. Values are uncertain by about ±0.035, and the fractions measured on
  the selection population are up to 0.08 higher.

The table gives the fraction of confident, correct redshifts:

| Strategy | Bright, same depth per tuning | Bright, fixed total time | Faint, same depth | Faint, fixed total time |
|---|---|---|---|---|
| Band 2 now, 1 tuning (best block starts at 82 GHz) | 0.06 | 0.06 | 0.04 | 0.04 |
| Band 2 now, 2 tunings | 0.21 | 0.21 | 0.18 | 0.14 |
| Band 2 now, 4 tunings | 0.68 | 0.56 | 0.51 | 0.18 |
| Band 2 now, 5 tunings (greedy) | 0.76 | — | 0.63 | — |
| Band 2 now, 5 contiguous (science-verification layout) | 0.77 | 0.65 | 0.64 | 0.16 |
| Band 2 now, full band (7 tunings) | 0.86 | 0.74 | 0.77 | 0.13 |
| Band 2 + WSU, 1 tuning | 0.66 | 0.66 | 0.61 | 0.61 |
| Band 2 + WSU, 2 tunings | 0.85 | 0.82 | 0.71 | 0.59 |
| Band 3 standard 5-tuning scan | 0.68 | 0.56 | 0.54 | 0.11 |

"Bright" means CO(4–3) at S/N 12 per tuning, and "faint" means S/N 5. Under fixed total time, per-line S/N falls as
1/√N. The receiver models follow C. De Breuck, "ALMA Band 2 as a redshift machine" (ESO, February 2026).

Greedy Band 3 layouts of one to four tunings, with the sidebands 4.125–7.875 GHz from the LO, are in
`results/redshift_tuning.json` and in the figure. Four tunings reach 0.49 for bright sources and 0.43 for faint
ones, below the five contiguous tunings of the standard scan.

Bakx (2026) asks the same question for Band 3 and for Band 2 with 8 and 14 GHz of bandwidth, using the High-z
Sudoku tool, which counts the lines a layout covers at each redshift. This optimiser scores simulated spectra at a
set S/N, and its receiver models differ. The two sets of answers have not been compared.
