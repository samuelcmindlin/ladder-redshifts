#!/bin/bash
# Download ALMA ADMIT summary files (*.admit.native.tgz) listed in a text file or in an archive download script,
# trying the three ALMA mirrors in turn (Europe/ESO, East Asia/NAOJ, North America/NRAO), two files at a time.
# Safe to rerun: complete files (gzip integrity check passes) are skipped; interrupted ones are resumed.
#
# Usage:  bash get_admit_files.sh <list.txt | downloadscript.sh> <destination folder>
# e.g.    bash scripts/alma/get_admit_files.sh data/admit/2016.1.00672.S_admit_native_files.txt alma_admit
# Dry run (only counts what is missing):  DRY=1 bash get_admit_files.sh ...
set -u
[ $# -eq 2 ] || { echo "usage: bash get_admit_files.sh <list.txt|downloadscript.sh> <destination folder>"; exit 1; }
LIST_SRC="$1"
export DEST="$2"
export DRY="${DRY:-0}"
mkdir -p "$DEST"
NAMES=$(grep -o 'member\.uid___[^"[:space:]]*\.admit\.native\.tgz' "$LIST_SRC" | sort -u)
TOTAL=$(echo "$NAMES" | grep -c .)
echo "$TOTAL files in the list"

fetch() {
  f="$1"
  cd "$DEST" || exit 1
  if [ -s "$f" ] && gzip -t "$f" 2>/dev/null; then return 0; fi
  if [ "$DRY" = "1" ]; then echo "missing $f"; return 0; fi
  for host in almascience.eso.org almascience.nao.ac.jp almascience.nrao.edu; do
    for try in 1 2; do
      if curl -fsS -C - --connect-timeout 30 --max-time 900 -o "$f" "https://$host/dataPortal/$f" 2>/dev/null; then
        if gzip -t "$f" 2>/dev/null; then echo "got     $f  ($host)"; return 0; fi
        rm -f "$f"                      # transfer finished but the file is corrupt: start it again
      fi
      sleep 15
    done
  done
  echo "FAILED  $f  (run the script again later)"
}
export -f fetch
echo "$NAMES" | xargs -P 2 -n 1 bash -c 'fetch "$0"'

ok=0; for f in $NAMES; do [ -s "$DEST/$f" ] && gzip -t "$DEST/$f" 2>/dev/null && ok=$((ok+1)); done
echo "complete: $ok of $TOTAL"
[ "$ok" -eq "$TOTAL" ] && echo "ALL FILES DOWNLOADED" || echo "Some files are still missing: run the script again."
