# Development notes: what went wrong and how it was fixed

This log records every error found during development, in the order found, with its effect on the results. All the
results in this repository were produced after these fixes. The regression tests in `tests/` guard the ones that real
data exposed.

## Found on real data

1. **Close lines were counted twice.** The likelihood normalisation √(Σwᵢ²) assumed independent matched-filter
   outputs. CO(7–6) and \[CI\](2–1) are only 998 km/s apart, so their filter outputs overlap.
   - **What happened.** One broad line could pass as a confident blend pair. SPT0054-41 has a single flat-topped
     CO(5–4) line of about 1100 km/s, and its ADMIT spectrum came out as z = 7.23 at confidence 1.00 instead of
     4.877.
   - **Fix.** T = Σwy/√(wᵀCw), with C the kernel overlap ([`method.md`](method.md), section 3).
   - **Effect on the synthetic atlas.** At S/N 32–64, confident errors fell from 4.5% to 1.2% for a single Band 3
     tuning, and from 1.4% to 0.2% for the Band 4 scan. Confident-correct fractions moved by 0.01 or less.
   - **Guarded by** `test_broad_line_not_split_into_blend_pair`.
2. **Channels in descending frequency.** The gridding passed channels to `np.interp` in file order. For
   lower-sideband windows, where frequency decreases, that returns nonsense without raising an error.
   - **What happened.** It surfaced only with the multi-tuning cubes, as six confident errors in the first
     point-source run.
   - **Fix.** Channels are now sorted on load and again before gridding. Earlier spectra were already ascending, so
     no earlier result changed.
3. **Aperture spectra on continuum-subtracted cubes with a 4″ beam.** The detection-grown aperture summed several
   beams of noise. Cut-outs centred on a known position now use a point-source spectrum at the centre
   (`realspec.spectrum_point`).
4. **ADMIT spectra are not spectra at the source.** A first reading took the ADMIT CubeSpectrum tables as on-source
   spectra in every window.
   - **What happened.** In line-free windows those spectra sit on noise spikes anywhere in the field (in a
     primary-beam-corrected cube, usually near the edge).
   - **Fix.** Such windows are now zero-filled, and the engine is also run without them to show what they
     contribute.
5. **The ADMIT zero-fill was too confident.** The first noise scale (×2.5) assumed ADMIT would catch any line above
   about 5σ per channel.
   - **Why that is wrong.** ADMIT needs at least 4 channels above 5σ, so lines with an integrated S/N of up to ~20
     slip through.
   - **Effect.** At ×2.5 the true redshift was under-weighted in three sources: SPT0403-58, SPT0652-55 and
     SPT2152-40. The scale is now ×10.
6. **ADMIT's detection flag is not reliable.** A brief rule kept ADMIT spectra only where ADMIT itself reported
   emission.
   - **What happened.** SPT0403-58's published CO(4–3) line makes the cube's brightest voxel (5.2σ) right on the
     source, yet ADMIT reports neither a line nor a source.
   - **Fix.** The rule was reverted. A spectrum within half a beam of the target is now used whatever ADMIT's flag
     says. `--require-detection` restores the strict rule.
7. **Archive products and downloads.**
   - **Duplicate windows.** ARI-L re-processed cubes duplicate a pipeline tuning, so they are skipped.
   - **Continuum images.** These were briefly processed as cubes; `.cont.` files are now excluded.
   - **Downloads.** The archive's own download script stalled with curl error 35 (a TLS failure) on one mirror.
     `scripts/alma/get_admit_files.sh` and `get_cutouts.sh` try all three ALMA mirrors in turn and skip complete
     files.

## Found on synthetic data, before any reported run

8. **S/N treated as flux.** The first trap reconstructions treated S/N ratios as flux ratios. That is wrong when
   windows come from different instruments (ATCA vs ALMA). The engine now takes per-window noise, and the weights are
   expected S/N.
9. **The [CI] prior was too narrow.** The first prior allowed only [CI]/CO L′ ratios of 0.15 and 0.45.
   - **What happened.** On HerBS-108-like data (two CO lines, no [CI]) the estimator was confidently wrong in 8–33%
     of noise draws.
   - **Fix.** A weak-[CI] template (ratio 0.05, prior weight 0.2) cut this to 0–3%.
   - **Both versions** are in `results/redshift_traps.json`. With the current likelihood (item 21) the figures
     are 1–25% and 0%.
10. **The confidence window split blend partners.** The first window (±0.003 in ln(1+z)) split the posterior between
    the CO(7–6) and \[CI\](2–1) readings, which undercounted confident errors. It is now ±0.006.
11. **An unrealistic redshift prior in the traps.** A flat prior up to z = 12 gave high-z \[CI\](2–1) readings of a
    single line too much weight. The traps now use SPT-like or Herschel-like priors, and also report the flat-prior
    result.
12. **HerBS-108's alternative redshift.** Cox et al. (2023) quote z = 1.25 for the CO(2–1)/(3–2) reading and print
    the upper line as 145.60 GHz. The line frequencies imply 1.358, and CO(6–5) at z = 3.7166 is at 146.60 GHz.
    The repository uses the recomputed values and says so where it cites the paper.

## Found in the pre-publication review (2 November 2026)

An automated code and methods review ([`pre-publication-review.md`](pre-publication-review.md)) found the following.
All are fixed. The real-data and synthetic results in this repository were regenerated afterwards. The point-by-point
response is in [`review-response.md`](review-response.md).

13. **Overlapping windows overwrote each other.** Adjacent tunings overlap by a few channels. The gridder wrote
    each window over the previous one, so the result depended on the order of the windows in the file.
    - **What happened.** Reversing the order changed SPT2349-52's best redshift from 2.675 to 3.900. It also cost
      SPT2152-40 its confident status.
    - **Fix.** Overlapping channels are now combined by inverse-variance weighting in flux units, and gaps inside
      a window are no longer bridged. Both gridders share the function `realspec.grid_windows`.
    - **Effect.** SPT2349-52 now ranks its published redshift first (3.900, 0.66), so 10 of the 11 single-line
      sources do (9 before). The detection rule now also settles SPT2152-40 (4 of 7 instead of 3).
    - **Guarded by** `test_gridding_is_independent_of_window_order` and
      `test_overlapping_windows_are_combined_by_inverse_variance`.
14. **The equal-weight baseline was noise-weighted.** It was documented as an equal-weight sum of line S/N, but it
    inherited the per-window noise weighting of the main likelihood.
    - **Fix.** It now weights every covered line's S/N equally.
    - **Effect.** It settles 4 of the 7 multi-line sources instead of 3, the extra one being SPT0436-40.
    - **Wording.** Its description no longer equates it with Lupu et al. (2012), whose procedure is more elaborate.
    - **Guarded by** `test_equal_weight_baseline_weights_line_snr_equally`.
15. **Simulated data outside the coverage reached the matched filter.** The synthetic studies filtered the whole
    simulated spectrum before applying each setup's coverage. Signal and noise just outside a window could
    therefore enter the filter output at the window's edge.
    - **Fix.** The filter now uses only observed pixels and is normalised to unit noise near edges and gaps. The
      expected S/N of a line near an edge is reduced by the matching filter efficiency (`docs/method.md`,
      section 3). Real spectra go through the same model; there, zeros outside the coverage used to make the
      filter output slightly too quiet at edges.
    - **Guarded by** `test_data_outside_the_coverage_have_no_effect` and
      `test_filter_output_has_unit_noise_at_window_edges`.
16. **Broad single lines were mistaken for a blend pair.** The widest matched-filter profile was 800 km/s.
    - **What happened.** A single line broader than about 1000 km/s was then fit better by CO(7–6) + \[CI\](2–1),
      which are 998 km/s apart. In the development population this caused most of the confident errors for a
      single Band 3 tuning (7 in 100 confident answers).
    - **Fix.** A 1200 km/s profile was added. On the same population the errors fell to 3 in 96.
17. **"Calibrated" was overstated.** The confidence was called calibrated, but the development atlas had 7 wrong
    answers among the 101 single-tuning cases with confidence ≥ 0.99.
    - **Fix: what it is called.** It is now described as a confidence score. Its reliability is measured and
      reported for every simulated setup (`results/redshift_atlas.json`, `reliability`) and for the real data.
    - **Fix: fresh simulations.** After the fixes above, the atlas was re-run with a new random seed, so the
      reported synthetic numbers come from sources that played no part in development.
    - **Fix: the tuning optimiser.** It now reports its chosen layouts on an independent population.
18. **The detection profile was chosen by position.** Detection used the second row of the filter bank, so a
    one-width bank crashed and a reordered bank changed the detector. It is now chosen by width (450 km/s).
19. **Extraction failures were reported as success.** `extract_cube_spectra.py` logged failed cubes but exited with
    status 0. It now prints a failure summary and exits with status 1, unless `--allow-failures` is given.

## Found in the second review (12 November 2026)

A second review, of the release as a whole, found the following. All are fixed. Every real-data and synthetic result
in this repository was regenerated afterwards, and the atlas was re-run with a new random seed. Items 25–27 were
closed in a second pass on the same day; the follow-up planner, the tuning optimiser and the real-data batch were
regenerated again after them.

20. **The noise of the cube-route spectra was too high.** In line-free channels, flux/noise should scatter by 1.
    It scattered by 0.76 on average for 2016.1.00672.S (0.80 at 84–90 GHz, falling to 0.72 at 108–115 GHz) and by
    0.60 for SPT0027-50. Adjacent channels were uncorrelated, so only the scale was wrong.
    - **Likely cause.** The extractor measures noise in an annulus 1.3–2.6 box radii from the source, in
      primary-beam-corrected cubes, where the noise rises away from the pointing centre. The error grows with box
      size (0.58 for SPT0549-56 with a 12″ box, 1.0 for SPT2334-50 with a 6″ box) and with frequency, as that
      explanation predicts. It has not been checked against the cubes.
    - **What happened.** Every S/N was about 1.3 times too low (1.7 for SPT0027-50). The second lines of
      SPT0520-53, SPT0652-55 and SPT2203-41 fell below the detection rule's 5σ threshold; they are at 6.5σ, 7.4σ
      and 5.6σ. The headline had the detection rule settling 4 of the 7 multi-line sources. It settles all 7.
    - **Fix.** `realspec.calibrate_noise` measures the scatter of flux/noise in each window's line-free channels,
      fits a straight line in frequency through the windows of a source, and multiplies the noise by the fitted
      value. `spectrum_point` and `spectrum` now do this. The shipped spectra, extracted earlier, were corrected
      once with `scripts/recalibrate_noise.py`; only their `noise_mjy` column changed. The factors are in
      `data/noise_calibration.csv`.
    - **Effect.** The engine's verdicts are unchanged: 7 of 7 multi-line sources confident and correct, no
      single-line source confident, the published z first in 10 of 11, nothing confidently wrong. The ADMIT-route
      spectra were not affected: their noise was already scaled to the line-free scatter.
    - **Guarded by** `test_shipped_cube_spectra_have_unit_noise_scatter`,
      `test_calibrate_noise_recovers_a_wrong_noise_scale` and `test_point_source_spectrum_has_calibrated_noise`.
21. **Velocity offsets of tens of km/s decided between octave-family readings.** The likelihood placed every line
    exactly at its predicted position. Readings in the same octave family predict the same lines 29–122 km/s apart
    (the ladder's anharmonicity), and at S/N of tens that difference dominated the likelihood.
    - **What happened.** Take SPT0314-44's real spectrum and add one simulated Band 6 tuning that covers only the
      line both candidates share (CO(9–8) at z = 2.933, CO(6–5) at z = 1.622, 76 km/s apart). Displacing that line
      by 25–76 km/s relative to the Band 3 line turned a correct answer at confidence 1.00 into a wrong answer at
      confidence 1.00. Offsets of that size between the lines of one galaxy are ordinary. `method.md` said the
      engine did not rely on the effect.
    - **Fix.** Each line may sit up to 150 km/s from its predicted position (`ladder.DV_MAX`): its filter output
      is the largest within that range. CO(7–6) and \[CI\](2–1) keep their separation when both are covered. The
      redshift quoted for a solution is the one at which its lines sit where predicted
      ([`method.md`](method.md), section 3).
    - **Effect on real data.** No verdict changed. In the test above the answer no longer depends on the offset.
    - **Effect on the synthetic atlas.** The simulated lines have no velocity offsets, so there the tolerance
      can only cost. On the same 2,500 sources and six setups, exact positions give 20 wrong answers among 8,931
      confident ones; the tolerance gives 31 among 8,939. The number of confident, correct answers is the same
      to within 3.
    - **The cost under noise.** On pure-noise Band 3 scans the confidence exceeds 0.5 in 4.8% of cases, against
      2.2% with exact positions; it reaches 0.95 in neither (`scripts/noise_only_check.py`).
    - **Guarded by** `test_line_offsets_within_the_tolerance_leave_the_likelihood_unchanged` and
      `test_shared_harmonic_answer_does_not_depend_on_velocity_offsets`.
22. **A predicted line could hide in a one-pixel gap.** A line counted as covered only if its predicted pixel was
    observed. Adjacent spectral windows often leave one-pixel gaps.
    - **What happened.** With the tolerance of item 21, a reading could shift by a pixel or two and park a line it
      predicts, but the data lack, inside such a gap. SPT0314-44 briefly gained a third candidate at z = 5.56 with
      9% of the weight: CO(5–4) for the detected line, with the missing CO(6–5) in a gap.
    - **Fix.** A line counts as covered when an observed pixel lies within the tolerance of its predicted position
      (within one pixel when no tolerance is used). Its noise and filter efficiency are read at the nearest
      observed pixel.
    - **Guarded by** `test_line_predicted_in_a_one_pixel_gap_still_counts`.
23. **The follow-up planner trusted lines that cannot decide.** Two errors.
    - **S/N treated as flux, again.** Predicted lines were scaled from the S/N of the observed line as if its
      noise were the median noise of the scan. SPT0314-44's line lies where two tunings overlap, where the noise
      is 0.72 of the median, so its predicted lines were 1.4 times too bright. They are now scaled from its flux.
    - **Shared lines counted as decisive.** Tunings that add only a line both candidates predict were scored at
      P = 1, through the effect of item 21. The planner now requires a line that one candidate predicts and the
      other does not; otherwise that candidate's success probability is 0.
    - **Effect.** For SPT0314-44, 48 of the 200 LO settings reached P ≥ 0.9 at the scan's noise level (59 before),
      and the tunings near 176, 264 and 352 GHz no longer count. The recommended Band 4 tuning was unchanged.
      These counts were for the sideband positions corrected in item 25; the current ones are given there.
    - **Guarded by** `test_planner_does_not_count_tunings_that_only_add_shared_lines`.
24. **CI did not do what the review response said.** The response to the first review stated that CI tested
    Python 3.9 and ran the installed package outside the checkout. The workflow still tested 3.10 and 3.12 from the
    checkout only, and did not install astropy, so the extraction test returned without testing. The workflow now
    does what was stated.
25. **Follow-up tunings that no ALMA receiver can make.** `followup.py` and the Band 3 model of `tuning.py` put
    the two sidebands 6–9.75 GHz either side of the LO in every band. ALMA's IF range is 4–8 GHz in Bands 3, 4, 5
    and 7 and 4.5–10 GHz in Band 6, and the Band 3 data in this repository sit 4.1–7.9 GHz from their LOs. Band 5
    was also taken to start at 163 GHz; it starts at 158 GHz.
    - **Fix.** Each tuning is two 3.75 GHz sidebands centred in its band's IF range: 4.125–7.875 GHz from the LO,
      or 5.375–9.125 GHz in Band 6. Every LO that keeps both sidebands inside the band is tried.
    - **Effect on the planner.** For SPT0314-44, 49 of 190 LO settings reach P ≥ 0.9 at the scan's noise level
      (48 of 200 with the old geometry). The recommendation is still Band 4, now at an LO of 139 GHz
      (136–142 GHz all work; it was 137–141 GHz).
    - **Effect on the tuning optimiser.** Only the greedy Band 3 layouts change. On the held-out population,
      four tunings now reach 0.49 for bright sources and 0.43 for faint ones (0.62 and 0.48 with the old
      geometry). The five-tuning Band 3 scan and every Band 2 number are unchanged.
    - **Guarded by** `test_follow_up_tunings_fit_the_alma_receivers`.
26. **SPT0314-44 was listed as unconfirmed.** `truth.csv` and the text called its published redshift tentative, as
    Reuter et al. (2020) do. Riechers (2025, ApJ 980, 59) confirmed it: CO(9–8), H₂O emission and OH⁺ and NH
    absorption in an ACA Band 6 observation give z = 2.9338 ± 0.0002.
    - **Fix.** Its status is now `confirmed`, with the confirming paper in the reference and notes.
      `z_published` keeps Reuter et al.'s 2.9345, so that `scripts/verify_truth_table.py` still checks every
      SPT row against one table; the two values differ by 53 km/s.
    - **Effect.** No number changes. The engine's verdict was and is "ambiguous, best = published". What changes
      is the statement about the truth: all 11 single-line redshifts, not 10, are confirmed by data the engine
      never saw.
    - **What it shows about the planner.** The confirming tuning is one the planner scores at 0
      ([`results.md`](results.md), section 4).
27. **The closest earlier method was not cited.** The template cross-correlation used with the Redshift Search
    Receiver on the Large Millimeter Telescope (Yun, Heyer & Aretxaga 2007; Yun et al. 2015, appendix A) weights
    many transitions by their relative strengths, as this engine does. It is now described in the README and in
    [`method.md`](method.md), section 3, with what this engine adds and the statement that the two have not been
    compared. Bakx (2026), which studies scan strategies for Band 3 and Band 2 with the High-z Sudoku tool, is
    now cited next to `tuning.py`. "Riechers 2024" is now cited as published (Riechers 2025, ApJ 980, 59).

## Comparisons added on 12 November 2026

These are new tests, not fixes. The second review asked for the first. The second was added because the real
spectra, as observed, cannot show the gain that the simulations claim.

- **The High-z Sudoku code on the same sources.** `scripts/compare_sudoku.py` runs the authors' `RSGquality`
  function, unmodified, at the published redshift of each real source and on the atlas setups
  ([`results.md`](results.md), section 8).
  - **What it showed.** No contradiction. Every source the engine settles is classed as robust by that code once
    \[CI\] is counted. The stand-in rule tracks the class it imitates to within 0.02 in four of six setups and is
    higher in the two Band 2 setups.
  - **What it cannot show.** That code reads no spectrum and has no notion of S/N, so error rates cannot be
    compared.
  - **A detail found on the way.** That code places ALMA's sidebands 4–7.75 GHz from the LO, which agrees with
    the geometry of item 25.
- **Noise added to the real multi-line spectra.** `degrade.py` ([`results.md`](results.md), section 9). As observed,
  the second lines of all seven sources are above 5σ, so the real data could not show the gain the simulations
  claim. Noise is added until the second line sits at a chosen strength.
  - **The guided rule was added to make the test fair.** Beating a rule that needs every line at 5σ, below 5σ,
    shows nothing. `ladder.guided_rule` accepts a weaker line where a surviving solution predicts one.
  - **The first write-up led with the wrong comparison.** It quoted the engine with the survey prior (0.92 at 4σ)
    against a rule that had no prior (0.80). Like for like, with no prior, the engine stands at 0.87 against
    0.80, and per source it is ahead for three of seven and behind for four. The README and
    [`results.md`](results.md) now lead with the no-prior figures and show the per-source picture.
  - **The rule was then given the prior too.** It answers when one remaining solution holds 0.95 of the prior
    weight. With a broad survey prior that almost never happens, so its figures do not move.
  - **A cross-correlation after Yun et al. (2015) was added** (`ladder.cross_correlation`), since it is the closest
    earlier method. Compared on best guesses with no prior, the engine is at 1.00 and the cross-correlation at
    0.96 with the second line at 4σ; 0.98 and 0.92 at 3.5σ. Multiplied by the survey prior, the
    cross-correlation fails for the two sources at z > 5, but the paper uses a per-source photometric redshift
    there, which this repository does not have.
  - **What the test supports.** Somewhat more answers than an observer-style rule below 5σ, unevenly across
    sources, at the same rate of wrong answers; and a best guess that is right slightly more often than a
    fixed-template cross-correlation's.
  - **What it does not.** Anything about unseen spectra, or about a practised observer. Seven sources, all from
    the development sample; draws that share the real spectrum's noise; comparison methods of our own making.
  - **Guarded by** `test_guided_rule_confirms_a_weak_line_and_otherwise_equals_the_detection_rule`,
    `test_guided_rule_uses_a_prior_only_when_one_solution_dominates`,
    `test_cross_correlation_noise_matches_shuffles_and_two_lines_give_the_redshift`,
    `test_degradation_reaches_the_requested_second_line_snr` and
    `test_sudoku_comparison_file_matches_the_current_engine_results`.

## Photometric redshifts as priors, added on 13 November 2026

A new test, not a fix ([`results.md`](results.md), section 10). The open issue it closes was: the Sudoku classes
call five single-line sources robust given a photometric redshift, and the engine had never been given one.

- **Where the values came from.** Table F.1 of Reuter et al. (2020), read twice: from the arXiv copy of the
  paper, and from the machine-readable table the SPT-SMG collaboration posts. All 81 rows agree.
  - **One row needs care.** The posted file replaces the SPT0452-50 row with values recomputed for its corrected
    redshift and keeps the printed row as a comment. `scripts/photoz_calibration.py` takes the printed
    photometric redshift (4.4 ± 0.7) from the comment.
  - **A second posted file** revises the photometric redshifts (2024). Both versions are run.
  - **The tables are not copied here** (no licence file). `data/photoz.csv` holds three columns for the 21
    sources this repository uses.
- **The printed errors were checked before being trusted.** For the 60 SPT sources outside this sample the
  2020 values run high by 0.48 errors and scatter by 1.10 errors. A third version of the prior applies that
  shift and widening. It lowers the engine's count from seven to four.
- **What it showed.** With the photometric prior the engine is confident and correct for four to seven of the
  eleven one-line sources, and confident and wrong for none. Rules given the same values settle about as many
  and made six wrong calls over ten combinations of rule and prior.
- **What it showed against the engine.** With the 2024 value for SPT2349-52, which is 2.8 errors off, the engine
  gives a wrong redshift 0.93. Nothing in the engine protects it from a prior that is far off; the 0.95
  threshold held by a small margin.
- **The comparison rule was weaker than it should have been.** The detection rule drops a redshift when a
  predicted \[CI\] line is missing. Given a photometric prior, the guided rule then returned a wrong redshift
  in 249 of 5,600 noise draws, nearly all for two sources whose \[CI\] line had sunk below 3σ. A form that
  counts CO lines only (`ci_prunes=False`) gives 6 wrong answers and more correct ones, and is the form the
  engine is now compared with when a photometric prior is used.
  - **Effect on the earlier comparison.** With no prior the CO-only form stands at 0.82 and 0.62 with the
    second line at 4σ and 3.5σ, against 0.80 and 0.60 for the form used until now. The engine's 0.87 and 0.73
    are unchanged; the README quotes both forms.
- **The cross-correlation with a per-source p(z)** is now run as the paper describes (S/N × p(z)). With these
  narrow dust-continuum priors it fails for four sources that have two detected lines. The paper's p(z) is
  broader, so this says how the product behaves, not how the paper's results would change.
- **`run_engine` was split** into `engine_likelihood` and `engine_posterior` so that several priors share one
  likelihood. Every number in the earlier `results/redshift_degradation.json` is reproduced exactly by the new
  code; only new columns were added.
- **Guarded by** `test_photometric_redshift_table_matches_its_calibration_record`,
  `test_photometric_redshift_priors_settle_single_line_sources_without_a_wrong_call` and the photometric part
  of `test_degradation_reaches_the_requested_second_line_snr`.

## The comparison rule re-examined, 13 November 2026: a claim withdrawn

Until this date the README and section 9 of the results said that, with no prior, the engine gave more answers
than "a rule that works as an observer would" once the second line fell below 5σ: 87% against 80% at 4σ. That
claim is withdrawn.

- **How it was found.** Giving the rule a photometric prior exposed a fault in it: it ruled a redshift out when
  a \[CI\] line was missing. Removing that assumption raised the rule from 80% to 82%. That raised the
  question of what its other thresholds were worth, so 72 settings were run.
- **What the family showed.** A setting that counts lines from 4.5σ and confirms from 2.5σ returns the published
  redshift in 89% of draws at 4σ, with the same 4 wrong answers in 5,600 as the engine. Per source the engine is
  behind it for five of seven sources at 4σ.
- **What was wrong with the first comparison.** One rule, with thresholds chosen by us before any test (5σ and
  3σ), was treated as "the" rule. The engine was never tuned against it, but the rule was never given a chance
  to be good either.
- **What survives.**
  - With a photometric prior the engine is ahead of all 216 settings at every number of wrong answers. The same
    scrutiny was applied before saying so: the share of prior weight at which a rule answers was added to the
    family for this reason.
  - With no prior, the engine at a lower threshold returns more answers than any setting once 6 or more wrong
    answers in 5,600 are tolerated.
  - The engine gives a number where a rule gives yes or no, and it needs no thresholds chosen per survey.
- **What does not.** Any statement that the engine extracts more redshifts than a detection rule from the same
  spectrum alone.
- **Guarded by** `test_the_rule_the_engine_is_compared_with_is_the_best_of_its_family_by_the_stated_criterion`.

## A test on unseen spectra, planned on 13 November 2026

- **The plan** is [`independent-test-plan.md`](independent-test-plan.md): 23 sources of ALMA project
  2015.1.00504.S, with the sample, the data handling, the outcomes and their reading fixed in tag `v0.2.0`.
- **Nothing of that sample has been downloaded.** The archive was queried only for the list of targets and the
  kinds of file it holds. The project has re-imaged (ARI-L) cubes, not pipeline cubes, so the manifest builder
  and the extractor learned their file names.
- **`LADDER_REDSHIFTS_SAMPLE`** redirects the data, results and figures folders to a sub-folder, so that the new
  sample is scored by the same code with no edits.
- **Bounds.** `photoz.upper_bound` gives the exact 95% upper bound on the wrong-answer rate: 31% from 8 confident
  answers, 18% from 15.

## The independent test, run on 14 November 2026

[`results.md`](results.md), section 11, has the outcome. For the record of how it was done:

- **Order of events.** Tag `v0.2.0` and the plan were pushed at 23:50 UTC on 13 November. The list of cubes was
  fetched at 23:51 and the 460 cut-outs between 23:51 and 23:58. Extraction and the run followed.
- **The archive cut these cubes.** Each cut-out is 80 × 80 pixels and 4 to 16 MB; 3.0 GB in all.
- **No change to the extraction.** The one thing prepared in advance, reading ARI-L file names, was in the tag.
- **One command of the plan could not run.** `degrade` assumes every multi-line source shows two lines at 5σ.
  That was true of the development sample and of 5 of the 13 here. A wrapper outside the package restricts the
  sources; the module is untouched. The assumption should have been caught when the plan was written.
- **What the test taught about the claims.**
  - The 0.95 threshold held: no confident wrong answer in 11 and 17 more.
  - "All multi-line sources are settled" was a property of a deep sample. Here it is 9 of 13.
  - The rule setting picked on the development draws lost ground on new sources, as expected of a setting
    picked by its own score. With no prior, engine and rule remain level when both samples are pooled.
  - The lead with a photometric prior carried over: 17 sources against 10 as observed, and all five sources at
    every noise level.
- **Guarded by** `test_the_independent_sample_results_are_those_of_the_current_code`.

## Known open issues

- **The comparison rule.** It is a simplified Sudoku-style rule, not the High-z Sudoku code. That code has now been
  run on the same coverage and redshifts (see "Comparisons added on 12 November" above). What remains undone:
  - **The reading by eye.** The code's plotting tool (`RSGplot`) shows the redshifts allowed by a list of detected
    lines, the coverage and a photometric redshift. Its authors' reading of the 18 spectra is not reproduced; the
    stand-in rule and the guided rule are this repository's substitutes.
  - **Photometric redshifts for the real sources.** Done on 13 November (see above): with the published
    photometric redshifts the engine settles four of the five sources that the Sudoku classes call robust given
    a photometric redshift. What remains is a prior that allows for a photometric redshift being far off, for
    example a mixture with a broad component; none is implemented.
  - **Licence.** The code carries no licence file, so it is not bundled; `scripts/compare_sudoku.py` takes the
    path of the user's own copy.
- **Point-source spectra.** These lose flux for extended or multiply-imaged sources.
- **SPT0549-56.** It has no published position or redshift and no detected continuum, so its cut-out is centred on the
  pointing centre.
- **Published redshifts.** These were transcribed from Reuter et al. 2020, Table E1, using `pdftotext -layout`. The
  same table is also in CDS catalogue J/ApJ/902/78. All 18 redshifts, their 3 mm line lists and the Table A1
  positions were re-checked row by row against the published PDF on 2 November 2026: 0 mismatches. Anyone with
  the PDF can repeat the check with `python scripts/verify_truth_table.py PDF`.
- **The confidence is a score, not a probability.** The likelihood is maximised over line amplitude rather than
  marginalised with an amplitude prior. Simulated setups are overconfident: 97.8% of answers above 0.95 are right
  for a single tuning and 99.4–99.9% in scans, while the stated values average 0.998–1.000. The real single-line
  sources are underconfident. Calibrating it needs a larger independent sample.
- **Line ratios between readings with the same line positions.** When two readings predict lines at the same
  frequencies (an octave pair whose distinguishing rungs all lie outside the coverage), only the prior and the
  line ratios separate them. The excitation prior is a grid of five J_c values, and line fluxes are taken as known
  to the noise level. At S/N of tens, the reading whose grid happens to match the measured ratio can then win by a
  large factor, so the confidence overstates what the ratio can tell.
  - **Example.** In the test of item 21, with z = 1.622 as the truth, the shared line moves the wrong reading
    (z = 2.933) from 0.90 to about 0.95.
  - **What is needed.** An uncertainty on each line flux (calibration between tunings and bands, scatter of
    real line ratios about the templates) inside the likelihood. A simple cap on each line's S/N was tried and
    rejected: it let one bright line, split into CO(7–6) + \[CI\](2–1), beat the single-line reading.
- **Window placement in the planner.** Each LO stands for one placement of the four spectral windows, centred in
  the IF range. Real set-ups can move them within it, and the planner does not search over that.
- **Lines outside the line list.** Only CO, \[CI\] and optionally \[CII\] are modelled. SPT0402-45 shows HCN(4–3)
  at 6.4σ; the engine can only treat such a line as noise or assign it to the wrong species. Reuter et al. (2020)
  report 29 such weaker features in the SPT spectra. SPT0314-44's redshift was settled by H₂O and OH⁺ lines that
  neither the engine nor its planner knows about (item 26).
- **Correlated channels.** The noise model treats grid pixels as independent. For channels near 50 km/s this
  holds closely: the measured filter noise is 0.96 for 50 km/s channels. Wider channels get an approximate
  √(width) correction; full covariance propagation is not implemented.
- **Model completeness.** Velocity offsets between lines beyond 150 km/s, non-Gaussian or multi-component lines,
  continuum residuals, several sources along one line of sight, and redshifts outside 0.5–8 are not modelled.
- **The cross-correlation is a re-implementation.** The method of Yun et al. (2015) (item 27) is compared only
  through `ladder.cross_correlation`, written from the paper's description, and without the per-source
  photometric redshifts the paper uses. The authors' code has not been run.
- **No external review yet.** All checks so far are internal or automated, including both reviews above.
