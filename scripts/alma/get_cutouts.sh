#!/bin/bash
# Download small cut-outs of ALMA pipeline cubes through the archive's SODA cut-out service:
# one cut-out per spectral window, centred on each source's published position (rows of a manifest built by
# make_cutout_manifest.py). Instead of a 62 MB cube, each cut-out is a few MB.
# Safe to rerun: complete files are skipped. Tries the three ALMA mirrors in turn, two files at a time.
#
# Usage, from the repository root:
#   bash scripts/alma/get_cutouts.sh --test   fetch ONE cut-out and check it (do this first; takes under a minute)
#   bash scripts/alma/get_cutouts.sh          fetch everything (priority 1 first, then 2, then 3)
#   bash scripts/alma/get_cutouts.sh 1        fetch priority 1 only
# Settings (optional, as environment variables):
#   MANIFEST=path/to/manifest.tsv     default: data/cutouts/2016.1.00672.S_cutout_manifest.tsv
#   DEST=path/to/folder               default: alma_cutouts/ in the repository root (about 3.9 GB; ignored by git)
#   FULL=1                            if a cut-out fails on every mirror, download the whole cube instead (62 MB)
set -u
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
export MANIFEST="${MANIFEST:-$ROOT/data/cutouts/2016.1.00672.S_cutout_manifest.tsv}"
export DEST="${DEST:-$ROOT/alma_cutouts}"
export FULL="${FULL:-0}"
[ -f "$MANIFEST" ] || { echo "manifest not found: $MANIFEST"; exit 1; }
mkdir -p "$DEST"

# a complete FITS file: starts with SIMPLE, size a multiple of 2880 bytes, more than 20 kB
valid() {
  f="$1"; [ -s "$f" ] || return 1
  [ "$(head -c 6 "$f")" = "SIMPLE" ] || return 1
  n=$(wc -c < "$f" | tr -d ' ')
  [ "$n" -gt 20000 ] && [ $((n % 2880)) -eq 0 ]
}
# image size from the FITS header, e.g. "98 x 98 x 238"
dims() {
  head -c 14400 "$1" | LC_ALL=C tr -c '[:print:]' ' ' | fold -w 80 | LC_ALL=C grep -a -E '^NAXIS[123] ' \
    | LC_ALL=C sed -E 's/^NAXIS[123] *= *([0-9]+).*/\1/' | paste -sd'x' - | sed 's/x/ x /g'
}
export -f valid dims

fetch() {
  line=$(sed -n "${1}p" "$MANIFEST")
  IFS=$'\t' read -r file kind target ra dec radius prio <<< "$line"
  out="$DEST/$file"
  valid "$out" && return 0
  rdeg=$(awk -v r="$radius" 'BEGIN { printf "%.6f", r / 3600 }')
  for host in almascience.eso.org almascience.nao.ac.jp almascience.nrao.edu; do
    if [ "$kind" = "cutout" ]; then
      url="https://$host/soda/sync?ID=$file&CIRCLE=$ra%20$dec%20$rdeg"
    else
      url="https://$host/dataPortal/$file"
    fi
    for try in 1 2; do
      if curl -fsSL --connect-timeout 30 --max-time 900 -o "$out.part" "$url" 2>/dev/null && valid "$out.part"; then
        mv "$out.part" "$out"
        echo "got     $file  [$(dims "$out"); $host]"
        return 0
      fi
      rm -f "$out.part"
      sleep 10
    done
  done
  if [ "$kind" = "cutout" ] && [ "$FULL" = "1" ]; then
    for host in almascience.eso.org almascience.nao.ac.jp almascience.nrao.edu; do
      if curl -fsSL --connect-timeout 30 --max-time 1800 -o "$out.part" "https://$host/dataPortal/$file" 2>/dev/null \
         && valid "$out.part"; then
        mv "$out.part" "$out"; echo "got     $file  [WHOLE CUBE $(dims "$out"); $host]"; return 0
      fi
      rm -f "$out.part"
    done
  fi
  echo "FAILED  $file  (run the script again later)"
}
export -f fetch

if [ "${1:-}" = "--test" ]; then
  n=$(awk -F'\t' 'NR > 1 && $2 == "cutout" { print NR; exit }' "$MANIFEST")
  IFS=$'\t' read -r file kind target ra dec radius prio <<< "$(sed -n "${n}p" "$MANIFEST")"
  echo "Test: one cut-out of $file"
  echo "      (radius $radius arcsec around $target at RA $ra, Dec $dec)"
  rm -f "$DEST/$file"
  fetch "$n"
  if valid "$DEST/$file"; then
    size=$(wc -c < "$DEST/$file" | tr -d ' ')
    nx=$(dims "$DEST/$file" | cut -d' ' -f1)
    echo "Received $((size / 1024)) kB, image $(dims "$DEST/$file") pixels (the whole cube is 62,400 kB, 256 x 256 x 238)."
    if [ "${nx:-999}" -lt 200 ]; then
      echo "The cut-out service works. Now run:  bash scripts/alma/get_cutouts.sh"
    else
      echo "The archive sent the WHOLE cube: it ignored the cut-out request (FULL=1 downloads whole cubes, ~26 GB)."
    fi
  else
    echo "The cut-out did not arrive. What the archive answered (first lines):"
    curl -sSL --max-time 120 "https://almascience.eso.org/soda/sync?ID=$file&CIRCLE=$ra%20$dec%20$(awk -v r="$radius" 'BEGIN { printf "%.6f", r / 3600 }')" 2>&1 | head -c 1500
    echo
    echo "Fallback: FULL=1 bash scripts/alma/get_cutouts.sh 1   (whole cubes, about 12 GB for priority 1)"
  fi
  exit 0
fi

MAXP="${1:-3}"
LINES=$(awk -F'\t' -v m="$MAXP" 'NR > 1 && $7 <= m { print NR }' "$MANIFEST")
TOTAL=$(echo "$LINES" | grep -c .)
echo "$TOTAL files with priority <= $MAXP; saving to $DEST"
echo "$LINES" | xargs -P 2 -n 1 bash -c 'fetch "$0"'

ok=0
for n in $LINES; do
  f=$(sed -n "${n}p" "$MANIFEST" | cut -f1)
  valid "$DEST/$f" && ok=$((ok + 1))
done
echo "complete: $ok of $TOTAL"
[ "$ok" -eq "$TOTAL" ] && echo "ALL FILES DOWNLOADED" || echo "Some files are still missing: run the script again."
