"""The noise-added test on a sample in which some multi-line sources show fewer than two lines at 5 sigma.

degrade.main() takes every source whose published scan has two or three lines, and for each it needs two
published lines detected in our spectrum, the second at 5.5 sigma or more, because noise can only be added.
All seven development sources meet that. In the independent sample 5 of 13 do: the Cycle 3 scans are shallower.
This wrapper restricts degrade.multi_line_sources() to the sources that qualify and then calls degrade.main()
unchanged. Nothing in ladder_redshifts/ is modified.

It is a departure from docs/independent-test-plan.md, which named `python -m ladder_redshifts.degrade`; that
command stops with an error on this sample.

Usage: LADDER_REDSHIFTS_SAMPLE=independent python scripts/independent_noise_test.py [N_draws]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ladder_redshifts import degrade as D  # noqa: E402
from ladder_redshifts.paths import DATA  # noqa: E402
from ladder_redshifts.batch_real import read_spectrum, to_grid_auto  # noqa: E402
from ladder_redshifts.realspec import run_engine  # noqa: E402
from ladder_redshifts.ladder import prior_spt  # noqa: E402


def qualifies(src, z):
    r = run_engine(*to_grid_auto(read_spectrum(DATA / 'spectra' / f'{src}.csv')), prior=prior_spt)
    lines = D.published_lines(r['detections'], z)
    return len(lines) >= 2 and lines[1][2] >= max(D.TARGETS)


if __name__ == '__main__':
    every = D.multi_line_sources()
    keep = [(s, z) for s, z in every if qualifies(s, z)]
    print(f'{len(keep)} of {len(every)} multi-line sources have two published lines at 5 sigma, the second at '
          f'{max(D.TARGETS)} or more: ' + ', '.join(s for s, _ in keep))
    D.multi_line_sources = lambda: keep
    D.main()
