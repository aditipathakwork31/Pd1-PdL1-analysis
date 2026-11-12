#!/bin/bash
# Phase 1 — Structure Collection (corrected final set)
# Downloads all crystal structures + AlphaFold full-length models,
# then auto-prints TITLE/SOURCE lines for verification (don't skip reading these).
#
# Usage:
#   chmod +x phase1_download.sh
#   bash phase1_download.sh

BASE_DIR="$HOME/pd1_geometry/structures"
mkdir -p "$BASE_DIR"
cd "$BASE_DIR" || { echo "Could not cd into $BASE_DIR"; exit 1; }

echo "==> Saving files to: $BASE_DIR"
echo

PDB_LIST="
4ZQK:4ZQK_pd1_pdl1_complex.pdb
3RRQ:3RRQ_pd1_unbound.pdb
5B8C:5B8C_pd1_pembrolizumab.pdb
5WT9:5WT9_pd1_nivolumab.pdb
5X8L:5X8L_pdl1_atezolizumab.pdb
5X8M:5X8M_pdl1_durvalumab.pdb
5GRJ:5GRJ_pdl1_avelumab.pdb
"

for entry in $PDB_LIST; do
  id="${entry%%:*}"
  out="${entry##*:}"
  echo "Downloading $id -> $out"
  curl -fSL "https://files.rcsb.org/download/${id}.pdb" -o "$out"
  if [ $? -ne 0 ]; then
    echo "  !! FAILED: $id (check https://rcsb.org/structure/$id)"
  fi
  echo
done

echo "==> Crystal structures done."
echo

AF_LIST="
Q15116:PD1_alphafold_full.pdb
Q9NZQ7:PDL1_alphafold_full.pdb
"

for entry in $AF_LIST; do
  uid="${entry%%:*}"
  out="${entry##*:}"
  echo "Fetching AlphaFold metadata for $uid..."
  url=$(curl -fsSL "https://alphafold.ebi.ac.uk/api/prediction/${uid}" | grep -o '"pdbUrl" *: *"[^"]*"' | head -1 | sed 's/.*"pdbUrl" *: *"//;s/"$//')
  if [ -z "$url" ]; then
    echo "  !! Could not resolve URL for $uid — check https://alphafold.ebi.ac.uk/entry/$uid manually"
  else
    echo "  Found: $url"
    curl -fSL "$url" -o "$out" && echo "  OK: saved $out"
  fi
  echo
done

echo "==================================================="
echo "VERIFICATION — read every TITLE/SOURCE below."
echo "Confirm species = HOMO SAPIENS and partner matches"
echo "the master plan table before trusting any file."
echo "==================================================="
echo

for f in *.pdb; do
  echo "== $f =="
  grep -E '^(TITLE|SOURCE.*ORGANISM_COMMON)' "$f"
  echo
done

echo "==> Files in $BASE_DIR:"
ls -la "$BASE_DIR"
