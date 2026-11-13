#!/bin/bash
# Regenerate every result and figure from the data in this repository.
#   bash reproduce.sh          real-data results, the sonification and all figures (1-2 min)
#   bash reproduce.sh --full   also re-run the synthetic studies from scratch (about 1.7 hours on 2 cores)
# Logs go to logs/. Random seeds are fixed: with the same library versions a clean run leaves git status unchanged.
set -euo pipefail
cd "$(dirname "$0")"
PY="${PYTHON:-python3}"
mkdir -p logs
run() { echo "== $*"; "$PY" -m "$@" > "logs/$(echo "$*" | tr ' /' '__' | sed 's/ladder_redshifts\.//; s/--//g').txt" 2>&1; }

echo "== tests"
"$PY" -m pytest -q | tee logs/tests.txt

run ladder_redshifts.batch_real                 # real spectra: 2016.1.00672.S cube cut-outs + SPT0027-50
run ladder_redshifts.batch_real --admit         # the ADMIT-route pseudo-spectra
run ladder_redshifts.realspec                   # SPT0027-50, Band 2 science verification
run ladder_redshifts.ambiguity                  # undertone family, ambiguity function, sonification
run ladder_redshifts.summary                    # headline figure
if [ "${1:-}" = "--full" ]; then
  run ladder_redshifts.atlas 2500               # ~7 min
  run ladder_redshifts.traps 150                # ~7 min
  run ladder_redshifts.degrade 100              # ~10 min on 2 cores
  run ladder_redshifts.followup                 # ~15 min
  run ladder_redshifts.tuning 200               # ~55 min on 2 cores
  echo "== scripts/noise_only_check.py"; "$PY" scripts/noise_only_check.py > logs/noise_only.txt 2>&1   # ~5 min
else
  run ladder_redshifts.atlas --plot-only        # re-plot from results/redshift_atlas_raw.npz
  run ladder_redshifts.traps --plot-only        # re-plot from results/redshift_traps*.json/npz
  run ladder_redshifts.degrade --plot-only      # re-plot from results/redshift_degradation.json
  run ladder_redshifts.followup --plot-only     # re-plot from results/redshift_followup_SPT0314-44.json
  run ladder_redshifts.tuning --plot-only       # re-plot from results/redshift_tuning.json
  echo "(synthetic studies re-plotted from saved results; use --full to re-run them)"
fi
echo
grep -h "^summary:" logs/batch_real.txt
echo "done; logs in logs/"
