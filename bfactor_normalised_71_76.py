#!/usr/bin/env python3
"""Item 13 (second part): chain-normalised B-factors for PD-1 residues 71-76 (and 78).
Needs: pip install biopython numpy. Downloads PDB files from RCSB if not in the working folder.
Usage: python3 bfactor_normalised_71_76.py
Numbering: UniProt numbering is derived from the DBREF record (offset = dbSeqBegin - seqBegin) and
sanity-checked (UniProt 74 must be ASN, 78 must be LYS). If a check fails, set MANUAL_OFFSET below."""
import os, urllib.request, warnings
import numpy as np
from Bio.PDB import PDBParser
warnings.filterwarnings("ignore")

TARGETS = {"4ZQK": ["B"], "3RRQ": ["A"], "5B8C": ["C", "F", "I", "L"], "5WT9": ["G"]}
MANUAL_OFFSET = {}          # e.g. {("5WT9", "G"): 0}   offset = UniProt number - PDB author number
RES = [71, 72, 73, 74, 75, 76, 78]
UNIPROT = "Q15116"

def get(pdb):
    fn = f"{pdb}.pdb"
    if not os.path.exists(fn):
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb}.pdb", fn)
    return fn

def dbref_offsets(fn):
    offs = {}
    for line in open(fn):
        if line.startswith("DBREF ") and line[33:41].strip() == UNIPROT:
            chain = line[12]
            offs[chain] = int(line[55:60]) - int(line[14:18])
    return offs

print("structure chain offset  | z-score (B vs chain mean) at 71 72 73 74 75 76 78 | mean raw B 71-76 (CA / all atoms) | checks")
for pdb, chains in TARGETS.items():
    fn = get(pdb)
    s = PDBParser(QUIET=True).get_structure(pdb, fn)[0]
    offs = dbref_offsets(fn)
    for ch in chains:
        off = MANUAL_OFFSET.get((pdb, ch), offs.get(ch))
        if off is None:
            print(f"{pdb} {ch}: no DBREF offset found; set MANUAL_OFFSET"); continue
        chain = s[ch]
        cab = np.array([r["CA"].get_bfactor() for r in chain if r.id[0] == " " and "CA" in r])
        mu, sd = cab.mean(), cab.std()
        by_uni = {r.id[1] + off: r for r in chain if r.id[0] == " "}
        z, notes = [], []
        for u in RES:
            r = by_uni.get(u)
            z.append("  n/a" if r is None or "CA" not in r else f"{(r['CA'].get_bfactor() - mu) / sd:5.2f}")
        ok74 = by_uni.get(74) is not None and by_uni[74].get_resname() == "ASN"
        ok78 = by_uni.get(78) is not None and by_uni[78].get_resname() == "LYS"
        notes.append("numbering OK" if ok74 and ok78 else "NUMBERING CHECK FAILED (74 should be ASN, 78 LYS)")
        loop = [by_uni[u] for u in range(71, 77) if u in by_uni and "CA" in by_uni[u]]
        ca_mean = np.mean([r["CA"].get_bfactor() for r in loop]) if loop else float("nan")
        all_mean = np.mean([a.get_bfactor() for r in loop for a in r]) if loop else float("nan")
        print(f"{pdb} {ch} {off:+d} | {' '.join(z)} | {ca_mean:5.1f} / {all_mean:5.1f} | {'; '.join(notes)}")
print("\nManuscript reports mean B of 71-76 as: 72 (4ZQK), 69 (3RRQ), 45-69 (5B8C copies), 70 (5WT9). Compare with the raw columns.")
print("z > ~1 means the segment is more mobile than that chain's average; around 0 or below means unremarkable.")
