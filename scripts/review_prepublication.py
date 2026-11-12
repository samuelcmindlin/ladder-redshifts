#!/usr/bin/env python3
"""Reproduce the diagnostic checks in docs/pre-publication-review.md.

Updated after the fixes (docs/review-response.md) for the new API only: matched_filter now requires the coverage
mask, so the filter calls pass it (through Engine.filter or cover=...), and since 12 November 2026 the comparator
check asks posterior() for the exact-position redshift (refine=...), as run_engine does. The checks themselves are
unchanged; run on the fixed code they show the findings resolved.

Run from a checkout: python scripts/review_prepublication.py
Optionally save evidence with --output PATH. This is a diagnostic report, not a
passing regression suite. It does not modify the analysis or canonical results.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np

from ladder_redshifts.atlas import load_raw
from ladder_redshifts.batch_real import read_spectrum, read_truth, to_grid_auto
from ladder_redshifts.ladder import (
    Engine, Grid, TOL_LN, line_list, matched_filter, prior_spt, templates_grid,
)
from ladder_redshifts.realspec import run_engine


def window_order():
    rows = []
    for source in ("SPT2349-52", "SPT2152-40", "SPT2311-45"):
        spec = read_spectrum(ROOT / "data/spectra" / f"{source}.csv")
        forward = run_engine(*to_grid_auto(spec), prior=prior_spt)["sled"]
        reverse = run_engine(*to_grid_auto(spec[::-1]), prior=prior_spt)["sled"]
        rows.append(dict(source=source, forward=forward, reverse=reverse))
    return rows


def outside_coverage():
    grid = Grid(90, 110)
    cover = grid.mask([(95, 100)])
    templates, logw = templates_grid()
    engine = Engine(grid, cover, templates, logw=logw)
    original = np.zeros(grid.npix)
    changed = original.copy()
    changed[np.flatnonzero(cover)[-1] + 1] = 100
    before = engine.loglike(engine.filter(original))
    after = engine.loglike(engine.filter(changed))
    return dict(
        covered_data_identical=bool(np.array_equal(original[cover], changed[cover])),
        max_loglike_change=float(np.max(np.abs(after - before))),
        note="Only an unobserved pixel changes; the log likelihood should be invariant.",
    )


def equal_snr_comparator():
    truth = read_truth()
    saved = {r["source"]: r for r in json.loads(
        (ROOT / "results/redshift_real_batch.json").read_text())}
    rows = []
    for source, entry in truth.items():
        if not any(label in entry["notes"] for label in ("(2 lines)", "(3 lines)")):
            continue
        spec = read_spectrum(ROOT / "data/spectra" / f"{source}.csv")
        grid, data, sigma, cover = to_grid_auto(spec)
        engine = Engine(grid, cover, np.ones(len(line_list()[0])))
        filtered = engine.filter(data)
        z, confidence, _ = engine.posterior(
            engine.loglike(filtered), prior=prior_spt, refine=engine.loglike(filtered, dv_max=0))
        correct = bool(abs(np.log1p(z) - np.log1p(float(entry["z_published"]))) < TOL_LN)
        rows.append(dict(source=source, saved_flat_z=saved[source]["z_flat"],
                         saved_flat_confidence=saved[source]["conf_flat"],
                         equal_snr_z=float(z), equal_snr_confidence=float(confidence),
                         equal_snr_correct=correct))
    return dict(rows=rows, confident_correct=sum(
        r["equal_snr_correct"] and r["equal_snr_confidence"] >= 0.95 for r in rows))


def saved_calibration():
    raw = load_raw(ROOT / "results/redshift_atlas_raw.npz")
    result = raw["out"]["B3 single tuning"]["sled"]
    mask = result["conf"] >= 0.99
    error = abs(np.log1p(result["z"][mask]) - np.log1p(raw["z"][mask]))
    return dict(setup="B3 single tuning", n=int(mask.sum()),
                mean_confidence=float(result["conf"][mask].mean()),
                accuracy=float((error < TOL_LN).mean()),
                errors_outside_confidence_window=int((error >= 0.006).sum()))


def channel_spacing():
    rng = np.random.default_rng(42)
    rows = []
    for spacing in (50, 100, 200):
        freq = 90 * np.exp(np.arange(600) * spacing / 299792.458)
        window = dict(spw="noise", freq=freq, S=np.zeros(len(freq)),
                      sig=np.ones(len(freq)), cont=np.zeros(len(freq)))
        _, _, _, cover = to_grid_auto([window])
        indices = np.flatnonzero(cover)
        pixel = indices[len(indices) // 2]
        values = []
        for _ in range(300):
            window["S"] = rng.normal(size=len(freq))
            _, data, _, _ = to_grid_auto([window])
            values.append(matched_filter(data, cover=cover)[1, pixel])
        rows.append(dict(native_kms=spacing, measured_filter_noise_std=float(
            np.std(values, ddof=1))))
    return rows


def extraction_failure():
    from astropy.io import fits
    with tempfile.TemporaryDirectory(prefix="ladder-extraction-review-") as directory:
        directory = Path(directory)
        cube = directory / "test.TEST_sci.spw0.cube.I.pbcor.fits"
        fits.PrimaryHDU().writeto(cube)
        proc = subprocess.run([
            sys.executable, str(ROOT / "scripts/alma/extract_cube_spectra.py"), str(cube),
            "--out", str(directory / "out"), "--scratch", str(directory / "scratch"),
        ], capture_output=True, text=True, check=False)
        return dict(exit_code=proc.returncode, logged_failure="FAILED" in proc.stdout,
                    produced_spectra=len(list((directory / "out").rglob("*.npz"))))


def single_filter():
    spec = read_spectrum(ROOT / "data/spectra/SPT0314-44.csv")
    try:
        run_engine(*to_grid_auto(spec), prior=prior_spt, bank=(450.0,))
        return dict(error=None)
    except IndexError as error:
        return dict(error=str(error))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if os.environ.get("LADDER_REDSHIFTS_ROOT") and Path(
            os.environ["LADDER_REDSHIFTS_ROOT"]).resolve() != ROOT:
        parser.error("Unset LADDER_REDSHIFTS_ROOT so all checks use this checkout.")
    report = dict(window_order=window_order(), outside_coverage=outside_coverage(),
                  equal_snr_comparator=equal_snr_comparator(),
                  saved_calibration=saved_calibration(), channel_spacing=channel_spacing(),
                  extraction_failure=extraction_failure(), single_filter=single_filter())
    content = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(content)
    print(content, end="")


if __name__ == "__main__":
    main()
