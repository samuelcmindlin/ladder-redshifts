# Response to the pre-publication review

The automated review of 2 November 2026 is in [`pre-publication-review.md`](pre-publication-review.md), with its
evidence in [`pre-publication-evidence.json`](pre-publication-evidence.json). Every finding was reproduced before it
was fixed, using the reviewer's script, [`scripts/review_prepublication.py`](../scripts/review_prepublication.py). Its output on the fixed code is
[`post-fix-evidence.json`](post-fix-evidence.json).

**Note, 12 November 2026.** The numbers in this response are those of 2 November. A second review then corrected
the noise scale of the spectra and changed the likelihood (items 20–27 of
[`development-notes.md`](development-notes.md)), and every result was regenerated again. Current numbers are in
[`results.md`](results.md). Two things below changed as a result: the detection rule and this engine now both
settle all 7 multi-line sources, and [`post-fix-evidence.json`](post-fix-evidence.json) was regenerated with the
current code. The CI row describes the workflow as it is now; on 2 November the workflow had not been updated to
match it (item 24).

After the fixes, all real-data and synthetic results were regenerated:

- **The synthetic atlas** was re-run with a random seed not used during development.
- **The tuning optimiser** now reports its chosen layouts on an independent population.
- **New regression tests.** Eight were added (17 in all), one or more per finding.

Items 13–19 of [`development-notes.md`](development-notes.md) record each fix with its effect.

## P1 findings

| Finding | Status | What changed | Effect on the results |
|---|---|---|---|
| Overlapping windows made inference depend on input order | **Fixed** | Overlapping channels are combined in flux units with inverse-variance weights, as independent measurements (adjacent tunings and spectral windows have independent noise). Both gridders share `realspec.grid_windows`. Gaps inside a window are no longer bridged. | The real results no longer depend on window order (tested with reversed and shuffled windows). The 7 of 7 multi-line result holds. SPT2349-52 now ranks its published z first (3.900, 0.66), so single-line sources rank the published z first in 10 of 11 (9 before). The detection rule now settles 4 of 7 (3 before). |
| The equal-weight comparator was implemented with unequal weights | **Fixed** | The baseline now puts equal weights on the S/N of every covered line, with no noise weighting. It is described as a simplified stand-in for line-sum methods, not as Lupu et al.'s procedure, and the same definition is used in the real-data runs and the traps. | It now settles 4 of the 7 multi-line sources (3 before); the headline table and figure are updated. It is still never confident on a single-line source. |
| Synthetic observations used information outside their coverage | **Fixed** | The matched filter now uses only observed pixels. Its output is normalised to unit noise near edges and gaps, and the expected S/N of each line is reduced by the filter efficiency at its position. The reviewer noted that zeroing alone is incomplete; this covers the norm and the weights too. Real and simulated spectra use the same observation model, and the atlas, traps, tuning optimiser and ambiguity function were all re-run. | The reviewer's invariance test now gives a log-likelihood change of exactly 0. A test changes every unobserved pixel and checks that nothing moves. A Monte Carlo test checks unit noise at a window edge. |
| The general calibrated-confidence claim contradicted the simulations | **Fixed in the wording; measured** | "Calibrated confidence" was removed from the README and citation metadata. The quantity is now described as a confidence score, built from a profile likelihood (maximised over amplitude) and normalised over redshift. Its reliability is measured for every simulated setup on a fresh seed (`reliability` in `results/redshift_atlas.json`, and `docs/results.md` section 5) and on the real data. A cause of most single-tuning errors was found and fixed: broad single lines were fit as a CO(7–6)/\[CI\](2–1) pair, so a 1200 km/s profile was added. | On the fresh simulations, confident answers are right in 99.7–100% of cases in scans. For a single Band 3 tuning the figure is 97.1% (3 of 104 wrong), against a stated mean of 0.999, so the score remains overconfident there. The README states this. On the real single-line sources the score is cautious instead. |

## P2 findings

| Finding | Status | What changed |
|---|---|---|
| Extraction failures could be reported as success | **Fixed** | `extract_cube_spectra.py` counts failures, logs `INCOMPLETE`, and exits with status 1. A deliberately partial run can pass `--allow-failures`. Both the loose-FITS path and the `--tag freq` path are tested. |
| A one-width filter bank crashed `run_engine` | **Fixed** | Detection uses the profile closest to 450 km/s, chosen by value (`ladder.detection_index`). A one-width bank and a reordered bank are tested. |

## Other points raised

| Point | Response |
|---|---|
| **Channel covariance from resampling.** | `grid_windows` estimates each window's channel width. For runs of channels wider than 62.5 km/s it inflates the noise by √(width/50 km/s), with a warning. For channels much narrower than 50 km/s it warns that they should be binned. With this correction the reviewer's diagnostic gives matched-filter noise of 0.96, 0.95 and 0.97 for native 50, 100 and 200 km/s channels; before the fix the figures were 0.96, 1.44 and 1.88. All shipped spectra are 50–53 km/s. Full covariance propagation is listed as an open issue. |
| **Profile likelihood vs marginal posterior.** | Stated explicitly in `docs/method.md`. An amplitude prior has not been introduced. That is a modelling choice to make together with an independent calibration sample, not one to tune on the present sources. |
| **Follow-up planner approximations.** | The planner's approximations are listed in `docs/method.md` section 8: 10 draws, candidates below 1% dropped, success set to 0 when the noise-free case fails. Its outputs are described as planning estimates. |
| **Tuning layouts selected and scored on the same population.** | Each layout is now also scored on an independent population (new seed). Those held-out fractions are the ones reported and plotted. |
| **Truth table independently checked?** | All 18 redshifts and 3 mm line lists (Table E1), and all 18 positions (Table A1), were re-checked row by row against the published PDF of Reuter et al. 2020: 0 mismatches. `scripts/verify_truth_table.py` repeats the check for anyone with the PDF. CDS remained unreachable. |
| **CI and Python versions.** | CI tests Python 3.9 (the lowest supported version) and 3.12, with astropy installed so that the extraction test runs. It also installs the package and runs it outside the checkout. The test suite and `reproduce.sh` were also run locally on Python 3.9 and 3.13. |
| **`reproduce.sh` did not regenerate the trap figure.** | The traps now save the data behind their figure (`results/redshift_traps_examples.npz`), and `reproduce.sh` re-plots it. |
| **ALMA acknowledgement for North American authors.** | The NRAO / Green Bank sentence that ALMA asks North American authors to add is now included in the README and `data/README.md`. |

## Still open

- **External review.** No external scientific review has taken place yet.
- **Independent sample.** The calibration test needs about 100 new sources, run with a frozen version.
- **Comparison with the published methods.** Done since: the High-z Sudoku code, run on the same coverage and
  redshifts (it reads no spectrum; [`results.md`](results.md), section 8), and re-implementations of an
  observer-style rule and of the template cross-correlation of Yun et al. (2015) on noise-added real spectra
  (section 9). Still to do: the Sudoku authors' own reading of the spectra, a faithful Z-Spec procedure, and the
  cross-correlation authors' own code with per-source photometric redshifts.
- **Modelling gaps.** No amplitude prior, no full channel covariance, and no velocity offsets between species.
