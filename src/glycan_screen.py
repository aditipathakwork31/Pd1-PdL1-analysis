#!/usr/bin/env python3
"""Glycan SCREEN (not a glycan model).  For each N-glycosylation sequon Asn of the hybrid, report the distance from
the Asn ND2 (CB if ND2 absent) to the nearest heavy atom of the partner protein.
N-glycans extend well beyond the Asn (a GlcNAc2 core reaches roughly 1 nm; complex glycans more), so a site
within ~15 A of the partner is flagged as 'could contact partner'.  The 15 A figure is a heuristic screen, not a result.
Sequon positions are from UniProt (PD-1 N49,N58,N74,N116; PD-L1 N35,N192,N200,N219)."""
import argparse, pathlib
import numpy as np
from Bio.PDB import PDBParser, NeighborSearch

ROOT = pathlib.Path(__file__).resolve().parents[1]
SITES = {"P": [49, 58, 74, 116], "L": [35, 192, 200, 219]}
PARTNER = {"P": "L", "L": "P"}
FLAG_A = 15.0


def screen(model):
    rows = []
    for cid, sites in SITES.items():
        ns = NeighborSearch([x for r in model[PARTNER[cid]] if r.id[0] == " " for x in r])
        for n in sites:
            r = next((r for r in model[cid] if r.id[0] == " " and r.id[1] == n), None)
            if r is None:
                rows.append((cid, n, "absent", float("nan"), "")); continue
            at = r["ND2"] if "ND2" in r else r["CB"]
            near = ns.search(at.coord, 60.0, "A")
            d = min((float(np.linalg.norm(x.coord - at.coord)) for x in near), default=float("inf"))
            rows.append((cid, n, r.get_resname(), d, "could contact partner" if d <= FLAG_A else "clear of partner"))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hybrid", default=str(ROOT / "outputs/hybrid_crystal_core.pdb"))
    ap.add_argument("--out", default=str(ROOT / "outputs"))
    a = ap.parse_args()
    rows = screen(PDBParser(QUIET=True).get_structure("h", a.hybrid)[0])
    L = [f"{'protein':8}{'Asn':>5}{'residue':>9}{'ND2-to-partner (A)':>20}  verdict (flag <= {FLAG_A:.0f} A, heuristic)"]
    for cid, n, nm, d, v in rows:
        L.append(f"{'PD-1' if cid == 'P' else 'PD-L1':8}{n:>5}{nm:>9}{d:>20.1f}  {v}")
    txt = "\n".join(L); pathlib.Path(a.out, "glycan_screen.txt").write_text(txt + "\n"); print(txt)


if __name__ == "__main__":
    main()
