"""Repository paths. Scripts read data/ and write figures to figures/, numbers to results/ and sound to audio/.

The repository root is, in order: $LADDER_REDSHIFTS_ROOT if set; the folder above this package when it holds
data/truth.csv (a checkout or an editable install); otherwise the current directory.
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
DATA = ROOT / 'data'
RESULTS = ROOT / 'results'
FIG = ROOT / 'figures'
AUDIO = ROOT / 'audio'
if DATA.exists():
    for _d in (RESULTS, FIG, AUDIO):
        _d.mkdir(parents=True, exist_ok=True)
