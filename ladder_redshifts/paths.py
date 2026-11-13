"""Repository paths. Scripts read data/ and write figures to figures/, numbers to results/ and sound to audio/.

The repository root is, in order: $LADDER_REDSHIFTS_ROOT if set; the folder above this package when it holds
data/truth.csv (a checkout or an editable install); otherwise the current directory.

$LADDER_REDSHIFTS_SAMPLE names a second sample kept apart from the first: with it set to NAME, every module reads
data/NAME/ (truth.csv, photoz.csv, spectra/) and writes to results/NAME/, figures/NAME/ and audio/NAME/. No
other code changes, so a sample the engine was not developed on is scored by exactly the same code
(docs/independent-test-plan.md).
"""
import os
from pathlib import Path


def _root():
    env = os.environ.get('LADDER_REDSHIFTS_ROOT')
    if env:
        return Path(env).resolve()
    here = Path(__file__).resolve().parents[1]
    if (here / 'data' / 'truth.csv').exists():
        return here
    return Path.cwd()


ROOT = _root()
SAMPLE = os.environ.get('LADDER_REDSHIFTS_SAMPLE', '')
DATA = ROOT / 'data' / SAMPLE
RESULTS = ROOT / 'results' / SAMPLE
FIG = ROOT / 'figures' / SAMPLE
AUDIO = ROOT / 'audio' / SAMPLE
if DATA.exists():
    for _d in (RESULTS, FIG, AUDIO):
        _d.mkdir(parents=True, exist_ok=True)
