#!/bin/bash
# List the per-window cubes of a project's science targets in the ALMA archive: file names only, one per line,
# pipeline cubes and ARI-L cubes alike. The names are the input of make_cutout_manifest.py.
#
# Usage:  bash scripts/alma/list_project_cubes.sh <project code> [band] > list.txt
# e.g.    bash scripts/alma/list_project_cubes.sh 2015.1.00504.S 3 > data/cutouts/2015.1.00504.S_cube_files.txt
#
# It asks the archive's TAP service for the project's observing units in that band, and its DataLink service
# for the files of each. Tries the three ALMA mirrors in turn. Progress goes to stderr.
set -u
[ $# -ge 1 ] || { echo "usage: bash list_project_cubes.sh <project code> [band]" >&2; exit 1; }
PROJECT="$1"; BAND="${2:-3}"
QUERY="select distinct member_ous_uid from ivoa.obscore where proposal_id='$PROJECT' and band_list='$BAND'"
UIDS=""
for host in almascience.eso.org almascience.nao.ac.jp almascience.nrao.edu; do
  UIDS=$(curl -fsSL --connect-timeout 30 --max-time 180 -G "https://$host/tap/sync" \
           --data-urlencode "LANG=ADQL" --data-urlencode "FORMAT=csv" --data-urlencode "QUERY=$QUERY" 2>/dev/null \
         | tr -d '\r' | grep '^uid://' | sort -u)
  [ -n "$UIDS" ] && break
done
[ -n "$UIDS" ] || { echo "the archive did not answer the query for $PROJECT (band $BAND)" >&2; exit 1; }
echo "$(echo "$UIDS" | wc -l | tr -d ' ') observing units of $PROJECT in band $BAND (from $host)" >&2
missing=0
for u in $UIDS; do
  page=""
  for try in 1 2 3; do
    page=$(curl -fsSL --connect-timeout 30 --max-time 300 -G "https://$host/datalink/sync" \
             --data-urlencode "ID=$u" 2>/dev/null) && [ -n "$page" ] && break
    sleep 5
  done
  names=$(echo "$page" | grep -o 'member\.uid___[A-Za-z0-9_]*\.[^"<>/ ]*_sci\.spw[^"<>/ ]*\.cube\.I\.pbcor\.fits' | sort -u)
  n=$(echo "$names" | grep -c .)
  echo "  $u: $n cubes" >&2
  [ "$n" -gt 0 ] && echo "$names" || missing=$((missing + 1))
done | sort -u
[ "$missing" -eq 0 ] || echo "$missing observing units listed no cubes: run again, or check them in the archive" >&2
