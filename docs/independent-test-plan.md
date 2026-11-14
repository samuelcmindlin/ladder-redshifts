# Plan for a test on spectra the engine was not developed on

Written on 13 November 2026 and fixed in tag `v0.2.0`, before any spectrum of the sample was downloaded. The
outcome will be reported whatever it is, in `docs/results.md`, with every departure from this plan listed.

**Status:** run on 14 November 2026. The outcome and the departures are in
[`results.md`](results.md), section 11. The rest of this file is as it was in the tag.

## Why

Every result in this repository so far is on the 18 SPT sources of ALMA project 2016.1.00672.S and on SPT0027-50.
The engine's likelihood, noise calibration and thresholds were developed on those spectra. A result on them
cannot show how the engine does on spectra it has not seen.

## The sample

- **23 sources:** every science target of the Band 3 scans of ALMA project 2015.1.00504.S (Cycle 3) that has a
  redshift in Reuter et al. (2020, Table F.1) and is not in `data/truth.csv`. They are listed in
  `data/independent/truth.csv`, which `scripts/independent_sample.py` writes from the collaboration's posted tables.
- **By the lines Reuter et al. list from the 3 mm scan** (their Table E.1, read from the arXiv copy):
  - 13 sources with two or three lines;
  - 9 sources with one line;
  - 1 source with none in the scan used here: SPT0112-55, whose published lines come from a deeper Cycle 7 scan.
- **One redshift is tentative.** SPT0150-59 had one line and no second observation in 2020; its published
  redshift rested on the photometric redshift, as SPT0314-44's did.
- **Left out in advance:**
  - SPT0027-50, observed in the same project, because its Band 2 spectrum is in the development sample;
  - SPT0114-59, SPT2008-58 and SPT2052-56, observed in the same project, because Table F.1 has no redshift for
    them.
- **What was known when this was written.** The published redshifts, line lists, positions and photometric
  redshifts of these sources, all from the paper. No spectrum of any of them had been downloaded or looked at.
  The engine does not read the published redshifts; they are used only to score its answers.

## The data

- **Cubes.** The project predates per-window pipeline cubes. The archive holds re-imaged (ARI-L) cubes, four
  windows for each of five tunings per source: 460 cubes for the 23 sources.
  `bash scripts/alma/get_independent.sh` lists them, builds the manifest and downloads cut-outs around the
  published positions, or whole cubes if the archive does not cut these.
- **Extraction.** The script and settings used for the development sample:

  ```bash
  python scripts/alma/extract_cube_spectra.py alma_cutouts_independent --out data/independent --box 8 --tag freq \
         --targets data/cutouts/2015.1.00504.S_targets.csv
  ```

- **What may be adapted, and what may not.**
  - **May:** reading the files. Header keywords, axis order, units, channel width, and whether the cubes still
    contain continuum (the extraction already removes a linear continuum per window).
  - **May not:** anything in `ladder_redshifts/`. That covers the likelihood, the line list, the templates, the
    priors, the 0.95 threshold, the confidence window, the noise calibration and the comparison rules.
- **Exclusions.** A source is dropped only if its cubes cannot be obtained or do not cover its position. A source
  with a missing tuning is kept with the coverage it has.
- **Differences known in advance.** The Cycle 3 scans are shallower than the Cycle 4 ones (about 6 minutes on
  source against 12 to 15; 0.5 to 0.8 mJy/beam per channel against 0.4 to 0.6, by Reuter et al.'s account), and
  ARI-L imaging is not the Cycle 4 pipeline's.

## What is run, once

```bash
export LADDER_REDSHIFTS_SAMPLE=independent
python -m ladder_redshifts.batch_real     # spectra, noise calibration, engine with the survey prior
python -m ladder_redshifts.photoz         # each source's photometric redshift as the prior
python -m ladder_redshifts.degrade        # noise added to the multi-line spectra
```

`LADDER_REDSHIFTS_SAMPLE` only redirects the folders that are read and written (`data/independent/`,
`results/independent/`, `figures/independent/`). The code is that of tag `v0.2.0`.

## Outcomes, fixed in advance

A redshift is correct when it is within 0.003 in ln(1 + z) of the published one. The engine is confident at 0.95
or more.

**Primary**

1. **Confident answers that are wrong,** with the survey prior and with the 2020 photometric prior. Reported as
   k of N with the exact one-sided 95% upper bound on the rate, for this sample and pooled with the development
   sample (0 of 8 and 0 of 15 there).
2. **Multi-line sources (13):** number for which the engine is confident and correct with the survey prior.
   Development sample: 8 of 8.
3. **One-line sources (9):**
   - confident and correct with the 2020 photometric prior (development sample: 7 of 11);
   - best redshift equal to the published one with the survey prior (development sample: 10 of 11).

**Secondary**

4. **The comparison methods on the same spectra,** with the same priors: the two rule settings chosen on the
   development draws (below), the 5σ rule and the cross-correlation.
5. **Noise added to the multi-line spectra** (100 draws per source and level, as in section 9 of the results):
   - the engine at 0.95 against the rule setting chosen on the development draws, with no prior and with the 2020
     photometric prior, at 4σ, 3.5σ and 3σ;
   - sources for which the engine is ahead and behind, and a sign test pooled with the seven development sources.
6. **The 2024 revision of the photometric redshifts.** The recalibrated version is not interpreted for this
   sample: its two constants were fitted on sources that include it.

**The rule settings, chosen on the development draws and not re-chosen**

| Prior | Lines count from | A weaker line confirms from | Missing \[CI\] counts against | Missing line counts if expected above | Prior share to answer |
|---|---|---|---|---|---|
| none | 4.5σ | 2.5σ | no | 5σ | — |
| photometric | 4.5σ | 3.0σ | no | 8σ | 0.95 |

`degrade` evaluates all settings of the family on any sample; these two are the ones that will be quoted.

## How the outcomes will be read

- **The 0.95 threshold.** At a wrong-answer rate of 5% about one wrong answer is expected among 15 to 20
  confident ones. None or one is consistent with the threshold. Two or more with either prior would say the
  confidence is overstated on unseen spectra.
- **SPT0150-59 and SPT0112-55** are reported separately and are not counted in outcomes 2 and 3. A confident
  answer for SPT0150-59 that differs from its published redshift is a disagreement with a tentative redshift,
  and is reported as that.
- **No tuning afterwards.** If the engine is changed in the light of these spectra, that is a new version, and
  this sample is no longer independent for it.

## What this test cannot show

- **Small numbers.** 23 sources cannot rank methods that differ by a few sources. Bounding the wrong-answer rate
  below 5% takes about 60 confident answers with none wrong.
- **Same survey.** The sources are SPT-selected, observed with the same telescope and tuning scheme, and
  identified by the same team. A sample from another survey (z-GAL, for example) would be a harder test.
- **The author of the plan knew the published answers.** The code does not use them, but the test is not blind
  in the sense that a third party running it would be.
