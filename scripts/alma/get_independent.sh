#!/bin/bash
# Fetch the 3 mm scans of the independent test sample: the 23 sources of ALMA project 2015.1.00504.S named in
# docs/independent-test-plan.md. Run from your own terminal; the ALMA archive must be reachable.
#
#   bash scripts/alma/get_independent.sh          list the cubes, build the manifest, fetch one cut-out as a test,
#                                                 then fetch all 460 cut-outs if the archive honours cut-outs
#   FULL=1 bash scripts/alma/get_independent.sh   fetch whole cubes instead (about 44 MB each, about 20 GB in all)
#
# Files go to alma_cutouts_independent/ in the repository root (ignored by git). Safe to rerun.
set -u
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT" || exit 1
P=2015.1.00504.S
LIST="data/cutouts/${P}_cube_files.txt"
export MANIFEST="$ROOT/data/cutouts/${P}_cutout_manifest.tsv"
export DEST="${DEST:-$ROOT/alma_cutouts_independent}"
[ "${FULL:-0}" = "0" ] && export FULL=0 || export FULL=2      # FULL=1 here means whole cubes only

if [ ! -s "$LIST" ]; then
  echo "Listing the cubes of $P in the archive (about a minute) ..."
  bash scripts/alma/list_project_cubes.sh "$P" 3 > "$LIST.part" && [ -s "$LIST.part" ] && mv "$LIST.part" "$LIST" \
    || { rm -f "$LIST.part"; echo "Could not list the cubes. Nothing was downloaded."; exit 1; }
fi
echo "$(grep -c . "$LIST") cubes listed in $LIST"
python3 scripts/alma/make_cutout_manifest.py --ari-l --files "$LIST" --targets "data/cutouts/${P}_targets.csv" \
        --out "$MANIFEST" || exit 1
ROWS=$(($(grep -c . "$MANIFEST") - 1))
echo "$ROWS cubes in the manifest (460 expected: 23 sources x 5 tunings x 4 windows)"

if [ "$FULL" = "0" ]; then
  bash scripts/alma/get_cutouts.sh --test | tee "$DEST.test.log"
  if ! grep -q "The cut-out service works" "$DEST.test.log"; then
    echo
    echo "Cut-outs are not available for these cubes. Whole cubes are about 44 MB each, about 20 GB in all."
    echo "To fetch them:  FULL=1 bash scripts/alma/get_independent.sh"
    exit 2
  fi
fi
bash scripts/alma/get_cutouts.sh
