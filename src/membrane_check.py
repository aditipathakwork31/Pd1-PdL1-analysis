#!/usr/bin/env python3
"""All-atom membrane-plane check on the hybrid (chain P = PD-1, chain L = PD-L1), 4ZQK frame.
Membrane midplane point/normal per protein from geometry.py (TM C-alpha centroid; sign-fixed SVD axis).
delta = assumed distance from the midplane to the lipid-headgroup surface (scanned, A).
  own    : every heavy atom of the ectodomain must lie at height >= delta above its OWN membrane
  partner: every non-TM heavy atom of the other protein must lie at height >= delta on the ectodomain side of
           this membrane's plane (i.e. it does not cross into the other membrane).
Half-space tests on each plane separately.  The two planes are NOT parallel in this frame (see span output)."""
import argparse, pathlib, sys
import numpy as np
from Bio.PDB import PDBParser

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import geometry as g

SEL = {"P": dict(tm=(171, 191), ecto=(35, 145), own=(25, 145), nontm=(25, 170)),
       "L": dict(tm=(239, 259), ecto=(19, 127), own=(19, 225), nontm=(19, 238))}


def coords(chain, lo, hi, ca_only=False):
    return np.array([a.coord for r in chain if r.id[0] == " " and lo <= r.id[1] <= hi for a in r
                     if (a.get_id() == "CA" or not ca_only) and (a.element or "").strip().upper() != "H"])


def min_height(xyz, p, n):
    return float(((np.asarray(xyz) - p) @ n).min())


def analyse(model, deltas=(15.0, 20.0, 25.0)):
    anchors = {}
    for cid, s in SEL.items():
        tm = coords(model[cid], *s["tm"], ca_only=True)
        anchors[cid] = (g.centroid(tm), g.tm_normal(tm, coords(model[cid], *s["ecto"], ca_only=True)))
    res = {"own_P": min_height(coords(model["P"], *SEL["P"]["own"]), *anchors["P"]),
           "own_L": min_height(coords(model["L"], *SEL["L"]["own"]), *anchors["L"]),
           "L_vs_P_membrane": min_height(coords(model["L"], *SEL["L"]["nontm"]), *anchors["P"]),
           "P_vs_L_membrane": min_height(coords(model["P"], *SEL["P"]["nontm"]), *anchors["L"])}
    return res, {d: all(v >= d for v in res.values()) for d in deltas}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hybrid", default=str(ROOT / "outputs/hybrid_crystal_core.pdb"))
    ap.add_argument("--out", default=str(ROOT / "outputs"))
    a = ap.parse_args()
    res, ok = analyse(PDBParser(QUIET=True).get_structure("h", a.hybrid)[0])
    L = ["Minimum heights of heavy atoms above the membrane midplane (nm):"]
    L += [f"  {k:18} {v / 10:6.2f}" for k, v in res.items()]
    L += ["", "Pass (all minimum heights >= delta) for assumed midplane-to-surface offset delta:"]
    L += [f"  delta {d / 10:.1f} nm: {'PASS' if v else 'FAIL'}" for d, v in ok.items()]
    txt = "\n".join(L); pathlib.Path(a.out, "membrane_check.txt").write_text(txt + "\n"); print(txt)


if __name__ == "__main__":
    main()
