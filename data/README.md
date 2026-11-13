# Data

Everything the analysis reads is in this folder. The large intermediate products (cube cut-outs, ADMIT archives) are
not shipped, but the scripts in `scripts/alma/` rebuild them from the public ALMA archive.

| Path | Contents |
|---|---|
| `truth.csv` | Published redshifts. Columns: `source`, `z_published`, `prior` (`spt` or `herschel`), `reference`, `notes` (which lines the published scan itself showed, and which came from follow-up), `status` (`confirmed`, or `tentative` for a published but unconfirmed identification). |
| `photoz.csv` | Photometric redshifts of the 21 SPT sources in `truth.csv`, from Reuter et al. (2020, ApJ 902, 78), Table F.1. Columns: `source`; `z_phot`, `z_phot_err` (as printed in the 2020 paper); `z_phot_2024`, `z_phot_2024_err` (the collaboration's 2024 revision). Written and checked by `scripts/photoz_calibration.py` from the tables the collaboration posts at <https://github.com/spt-smg/publicdata>, which are not copied here. For SPT0452-50 the 2020 value is the one printed in the journal (4.4 ± 0.7), which the posted file keeps as a comment. |
| `independent/truth.csv`, `independent/photoz.csv` | The sample of the planned independent test ([`docs/independent-test-plan.md`](../docs/independent-test-plan.md)): 23 sources of ALMA project 2015.1.00504.S, with the same columns as `truth.csv` and `photoz.csv`. Written by `scripts/independent_sample.py` from the SPT collaboration's posted tables. No spectra yet. |
| `cutouts/2015.1.00504.S_targets.csv` | Positions of those 23 sources (Reuter et al. 2020, Table A.1) for the cut-outs and the extraction. The list of cubes and the manifest are written by `scripts/alma/get_independent.sh` when it is run. |
| `spectra/<SOURCE>.csv` | One spectrum per source from the cube route. The 21 sources are from project 2016.1.00672.S, as point-source spectra at the published position; SPT0027-50 is an aperture spectrum. These are what `batch_real` scores. |
| `admit/spectra/<SOURCE>.csv` | ADMIT-route pseudo-spectra of the same 21 sources (`batch_real --admit`). |
| `admit/<SOURCE>_admit.json` | Per-window summary of what ADMIT gave for each source: target position, which windows were measured, zero-filled or dropped, and why. |
| `admit/2016.1.00672.S_admit_native_files.txt` | The 436 ADMIT archive files of the project, the input list for `get_admit_files.sh`. |
| `cutouts/2016.1.00672.S_targets.csv` | Cut-out centre, priority, cut-out radius and extraction box for each source. |
| `cutouts/2016.1.00672.S_cutout_manifest.tsv` | The 435 archive files to cut (420 cubes plus 15 continuum images), the input for `get_cutouts.sh`. |
| `cutouts/extract_log.txt` | Log of the extraction that produced `spectra/` (channels, frequency range, box and noise per window). |
| `noise_calibration.csv` | The factor applied to the extractor's noise in each window of every cube-route spectrum (see "Noise" below). Columns: `spectrum`, `spw`, `freq_ghz`, `channels`, `scatter_measured` (scatter of line-free flux/noise in that window before the correction), `scale_applied`. |
| `SPT0027-50_band2_spectrum.csv` | SPT0027-50, ALMA Band 2 science verification: aperture spectrum (`realspec`). |
| `SPT0027-50_band2_spectrum_circle3arcsec.csv` | The same source from a fixed 3″ circle, for comparison. |
| `SPT0027-50_band2_collapsed.npz` | Collapsed image and aperture mask of SPT0027-50, used for the figure. |

## Spectrum format

The spectra are CSV files with one row per channel:

```
spw,freq_ghz,flux_mjy,noise_mjy,continuum_mjy[,admit]
```

- **Rows.** Channels of about 50 km/s, grouped by spectral window (`spw`), in any frequency order. Windows may
  overlap; overlapping channels are combined by inverse-variance weighting, independent of the order of the
  windows. Channels much wider or narrower than 50 km/s should be rebinned (see `docs/method.md`, section 7).
- **Flux.** `flux_mjy` is continuum-subtracted, and `continuum_mjy` is the continuum that was removed.
- **Noise.** `noise_mjy` is the per-channel noise at the source. In the cube-route spectra it is the extractor's
  annulus noise times the factor in `noise_calibration.csv`. The annulus sits off the source in
  primary-beam-corrected cubes, so its noise is about 1.3 times too high (about 1.7 for SPT0027-50); the factor
  puts it on the scale of the spectrum's own line-free channels, where flux/noise must scatter by 1
  (`docs/method.md`, section 7).
- **ADMIT column.** The optional `admit` column marks ADMIT-route windows as `measured` or `zero-filled`.

Any spectrum in this format can be scored. Put it in `spectra/`, add a row to `truth.csv` if the redshift is known, and
run `python -m ladder_redshifts.batch_real`. The noise column has to be right: `batch_real` prints a note when the
line-free channels of flux/noise do not scatter by about 1, and `realspec.calibrate_noise` rescales it.

## Sources of truth and positions

- **Published redshifts for 2016.1.00672.S.** These come from Reuter et al. 2020 (ApJ 902, 78), Table E1, which
  covers 18 of the 21 sources; SPT0433-59, SPT0549-56 and SPT2334-50 are not in that paper.
  - **How they were taken.** Extracted from the published PDF with `pdftotext -layout`. The table is also
    machine-readable at CDS, catalogue J/ApJ/902/78.
  - **Lines per source.** In 11 of the 18 the 3 mm scan itself showed a single line; the other 7 showed two or three.
    The `notes` column records which, and which follow-up lines in other bands confirmed the redshift. The tests
    and the summary figure take the single-line and multi-line groups from this column.
- **Positions.** These come from Reuter et al. 2020, Table A1. They agree with ADMIT's own source-finder positions to
  within 0.7″. Some sources sit well off the pointing centre: SPT0226-45 by 18″, SPT0425-40 by 9″ and SPT2316-50
  by 8″.
  - **Sources not in that paper.** SPT0433-59 and SPT2334-50 are centred on their continuum peak. SPT0549-56 has no
    detected continuum, so it is centred on the pointing centre.
- **Other entries in `truth.csv`.** SPT0027-50 is from Bakx 2026 (arXiv:2604.01089). The traps use
  SPT0452-50 (Riechers 2025), SPT0551-50 (Strandet et al. 2016) and HerBS-108 (Cox et al. 2023).

## Rebuilding `spectra/`: the cut-out route

Every pipeline cube in the ALMA archive is listed by DataLink with two services: the whole file
(`https://almascience.eso.org/dataPortal/<cube>.fits`, 62.4 MB per window for this project), and a SODA cut-out
service (`https://almascience.eso.org/soda/sync?ID=<cube>.fits&CIRCLE=<ra> <dec> <radius_deg>`). The cut-out
service makes the job about seven times smaller.

```bash
bash scripts/alma/get_cutouts.sh --test   # one cut-out; reports whether the archive honoured the cut (98 x 98 x 238 px)
bash scripts/alma/get_cutouts.sh          # all 435 files, priority 1 first; ~3.9 GB into alma_cutouts/; safe to rerun
python scripts/alma/extract_cube_spectra.py alma_cutouts --out data --box 8 --tag freq \
       --targets data/cutouts/2016.1.00672.S_targets.csv
python -m ladder_redshifts.batch_real
```

- **Downloading.** `get_cutouts.sh` tries the three ALMA mirrors in turn (ESO, NAOJ, NRAO), two files at a time, and
  checks each file is complete FITS. `FULL=1` falls back to whole cubes.
- **Extracting.** The extractor centres each box on the source, averages channels to about 50 km/s, measures the
  per-channel noise in an annulus, and writes `data/<SOURCE>/<SOURCE>_fNNN.NNN.npz`, one file per window (about
  266 MB in all; ignored by git). With `--tag freq` the files are named by centre frequency, because the project's
  five tunings reuse the same window numbers.
- **Scoring.** `batch_real` then takes a point-source spectrum at the centre of each cut-out, rescales its noise,
  and rewrites `spectra/<SOURCE>.csv` and that spectrum's rows of `noise_calibration.csv`.
- **Spectra extracted before the noise calibration existed.** The shipped spectra were first written with the
  annulus noise. `python scripts/recalibrate_noise.py` applied the calibration to them once (only `noise_mjy`
  changed). Spectra rebuilt from cut-outs do not need it.

The manifest was built with:

```bash
python scripts/alma/make_cutout_manifest.py --files data/admit/2016.1.00672.S_admit_native_files.txt \
       --targets data/cutouts/2016.1.00672.S_targets.csv \
       --cont SPT0433-59,SPT0549-56,SPT2334-50 --out data/cutouts/2016.1.00672.S_cutout_manifest.tsv
```

**Two positions changed after the download.** The shipped manifest records what was actually downloaded. At that time
SPT0433-59 and SPT2334-50 had no known position, so their 30″ cut-outs were centred on the pointing centre. Their
continuum peaks were then found in the downloaded continuum images and written into the targets file, and the
extractor centres its boxes there. Rebuilding the manifest from the current targets file gives the same 435 files,
radii and priorities, centred on those peaks instead.

## The ADMIT route: spectra without downloading cubes

Each pipeline cube has an ADMIT summary next to it (`*.admit.native.tgz`, about 3 MB per window). For this project the
436 files come to about 1.2 GB.

```bash
bash scripts/alma/get_admit_files.sh data/admit/2016.1.00672.S_admit_native_files.txt alma_admit
python scripts/alma/admit_spectra.py alma_admit --out data/admit/spectra     # repeat until it prints ALL DONE
python -m ladder_redshifts.batch_real --admit
```

**What ADMIT gives:**

- The cube's per-channel noise.
- The source finder's detections, which give the source position.
- A spectrum at the source, but **only when** the cube's brightest voxel or moment-map peak falls on the source.

**What it cannot give.** In line-free windows its two spectra sit on noise spikes anywhere in the field. So the
output is a pseudo-spectrum:

- **Measured windows.** A spectrum within half a beam of the source is used, with its baseline refitted and its noise
  rescaled.
- **Zero-filled windows.** Every other window gets flux 0 and noise equal to the cube's per-channel noise × 10.

**Yield for 2016.1.00672.S.** On-source spectra for 12 of 21 sources. Of the 7 sources with two or three lines, ADMIT
caught both lines for only one, so the cube route is needed when the decisive line is faint.

## SPT0027-50 (ALMA Band 2 science verification, 2011.0.00025.SV)

- **Download.** Get `SPT0027_X57_Product.tgz` (16 GB) from the science-verification page:
  https://almascience.eso.org/alma-data/science-verification/spt0027_band2
- **Extract.** Run `python scripts/alma/extract_cube_spectra.py SPT0027_X57_Product.tgz --out data`. It streams the
  archive and keeps a ±12″ cut-out of each window; disk use peaks at about 1.3 GB.
- **Analyse.** `python -m ladder_redshifts.realspec` then builds the aperture spectrum and the figure. Without the
  cut-outs it runs from the saved spectrum in this folder.

## Other public data that would extend the test

| Project | Content | Notes |
|---|---|---|
| 2015.1.00504.S | SPT Band 3 scans (Cycle 3) | pipeline products expected |
| 2018.1.01254.S, 2019.1.00486.S | SPT deeper 3 mm scans and 2 mm follow-up | pipeline products expected |
| 2011.0.00957.S | SPT Cycle 0 Band 3 scans (Weiß et al. 2013; includes SPT0452-50 and SPT0551-50) | raw data only; needs re-processing |
| 2016.2.00133.S, 2018.1.00804.S, 2019.1.01477.S | BEARS (Herschel sources), Band 3 and Bands 3 + 4 | from Urquhart et al. 2022 |
| IRAM M18AB, D20AB (NOEMA) | z-GAL, including HerBS-108 | through the z-GAL collaboration |

## Acknowledgement

This work makes use of the following ALMA data: ADS/JAO.ALMA#2011.0.00025.SV and ADS/JAO.ALMA#2016.1.00672.S. ALMA
is a partnership of ESO (representing its member states), NSF (USA) and NINS (Japan), together with NRC (Canada),
NSTC and ASIAA (Taiwan), and KASI (Republic of Korea), in cooperation with the Republic of Chile. The Joint ALMA
Observatory is operated by ESO, AUI/NRAO and NAOJ. The National Radio Astronomy Observatory and Green Bank
Observatory are facilities of the U.S. National Science Foundation operated under cooperative agreement by
Associated Universities, Inc.
