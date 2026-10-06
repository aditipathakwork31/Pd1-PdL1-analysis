#!/usr/bin/env python3
"""
Robustness checks for the PD-1:PD-L1 membrane-anchored hybrid model (phase4_grafted_model.pdb).

Reads the grafted model (chain X = PD-1 AlphaFold full-length, chain Y = PD-L1 AlphaFold full-length, both already
placed in the 4ZQK frame by phase4_span_calculator.py) and reports:
  1. membrane normals with a DEFINED sign (pointing from the TM midpoint toward the protein's own ectodomain),
     the angle between them, and the three span estimates plus a common-normal estimate;
  2. tilt sensitivity: spans when each membrane normal is tilted by 0/10/20/30 degrees in 12 azimuths;
  3. the smallest mismatch angle reachable with tilts <= 30 degrees;
  4. heavy-atom clashes (<2.0 A) between X and Y, split by region, and with signal peptides / cytoplasmic tails removed;
  5. mean pLDDT (B-factor column of AlphaFold models) by region;
  6. whether either protein crosses the other's membrane plane.
Deterministic (no random numbers). Usage:  python span_robustness.py phase4_grafted_model.pdb
"""
import sys, collections
import numpy as np
from Bio.PDB import PDBParser, NeighborSearch

PD1 = dict(chain="X", signal=(1, 24), ecto=(35, 145), stalk=(146, 170), tm=(171, 191), cyto=(192, 288))
PDL1 = dict(chain="Y", signal=(1, 18), ecto=(19, 127), ecto2=(133, 225), stalk=(226, 238), tm=(239, 259), cyto=(260, 290))
TILTS = (0, 10, 20, 30)
AZIMUTHS = np.linspace(0, 2 * np.pi, 12, endpoint=False)


def ca(chain, lo, hi):
    return np.array([r["CA"].coord for r in chain if lo <= r.id[1] <= hi and "CA" in r])


def rot_about(axis, th):
    a = axis / np.linalg.norm(axis)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * K @ K


def perp(v):
    t = np.array([1, 0, 0]) if abs(v[0]) < 0.9 else np.array([0, 1, 0])
    u = np.cross(v, t)
    return u / np.linalg.norm(u)


def tilted(n, theta_deg, az):
    return rot_about(rot_about(n, az) @ perp(n), np.radians(theta_deg)) @ n


def main(path):
    m = PDBParser(QUIET=True).get_structure("g", path)[0]
    X, Y = m[PD1["chain"]], m[PDL1["chain"]]

    # ---- 1. normals and spans
    tm = {"X": ca(X, *PD1["tm"]), "Y": ca(Y, *PDL1["tm"])}
    ecto = {"X": ca(X, *PD1["ecto"]), "Y": ca(Y, *PDL1["ecto"])}
    pt = {k: v.mean(0) for k, v in tm.items()}
    n = {}
    for k in "XY":
        axis = np.linalg.svd(tm[k] - tm[k].mean(0))[2][0]
        n[k] = axis * np.sign(axis @ (ecto[k].mean(0) - pt[k]))      # sign fixed: points toward the ectodomain
    n1, n2, d = n["X"], n["Y"], pt["Y"] - pt["X"]
    ang = np.degrees(np.arccos(np.clip(n1 @ n2, -1, 1)))
    nc = (n1 - n2) / np.linalg.norm(n1 - n2)
    print("1. NORMALS AND SPANS (normals point from TM midpoint toward own ectodomain)")
    print(f"   angle between normals {ang:.1f} deg (180 = facing, parallel membranes); deviation from antiparallel {180 - ang:.1f} deg")
    print(f"   span along PD-1 normal {abs(d @ n1) / 10:.2f} nm | along PD-L1 normal {abs(d @ n2) / 10:.2f} nm | "
          f"direct {np.linalg.norm(d) / 10:.2f} nm | along common normal {abs(d @ nc) / 10:.2f} nm")

    # ---- 2/3. tilt sensitivity
    print("\n2. TILT SENSITIVITY (each normal tilted by theta in 12 azimuths, all 144 combinations; nm)")
    for th in TILTS:
        s1, s2, sc, mis = [], [], [], []
        for a1 in AZIMUTHS:
            for a2 in AZIMUTHS:
                N1, N2 = tilted(n1, th, a1), tilted(n2, th, a2)
                c = (N1 - N2) / np.linalg.norm(N1 - N2)
                s1.append(abs(d @ N1) / 10); s2.append(abs(d @ N2) / 10); sc.append(abs(d @ c) / 10)
                mis.append(180 - np.degrees(np.arccos(np.clip(N1 @ N2, -1, 1))))
        print(f"   tilt {th:>2} deg: along PD-1 normal {min(s1):.1f}-{max(s1):.1f}; along PD-L1 normal {min(s2):.1f}-{max(s2):.1f}; "
              f"common normal {min(sc):.1f}-{max(sc):.1f}; mismatch {min(mis):.0f}-{max(mis):.0f} deg")
    best = min(180 - np.degrees(np.arccos(np.clip(tilted(n1, t1, a1) @ tilted(n2, t2, a2), -1, 1)))
               for t1 in range(0, 31, 2) for a1 in np.linspace(0, 2 * np.pi, 24, endpoint=False)
               for t2 in range(0, 31, 2) for a2 in np.linspace(0, 2 * np.pi, 24, endpoint=False))
    print(f"   smallest mismatch reachable with each normal tilted <= 30 deg: {best:.1f} deg")

    # ---- 4. clashes
    def region(cfg, i):
        for name in ("signal", "ecto", "ecto2", "stalk", "tm", "cyto"):
            if name in cfg and cfg[name][0] <= i <= cfg[name][1]:
                return name
        return "linker"
    heavy = lambda ch: [a for r in ch if r.id[0] == " " for a in r if a.element != "H"]
    ns = NeighborSearch(heavy(Y))
    pairs = [(a.parent.id[1], b.parent.id[1], a - b) for a in heavy(X) for b in ns.search(a.coord, 2.0)]
    print(f"\n4. CLASHES: {len(pairs)} heavy-atom pairs PD-1 vs PD-L1 closer than 2.0 A")
    for (rx, ry), c in sorted(collections.Counter((region(PD1, i), region(PDL1, j)) for i, j, _ in pairs).items(), key=lambda x: -x[1]):
        print(f"   PD-1 {rx:<7} vs PD-L1 {ry:<7}: {c}")
    mature = [(i, j, dd) for i, j, dd in pairs if region(PD1, i) not in ("signal", "cyto") and region(PDL1, j) not in ("signal", "cyto")]
    print(f"   after removing signal peptides and cytoplasmic tails: {len(mature)} pairs")
    if mature:
        print(f"   closest remaining: {min(x[2] for x in mature):.2f} A; residue pairs (PD-1, PD-L1): {sorted({(i, j) for i, j, _ in mature})}")

    # ---- 5. pLDDT
    def mbf(ch, rng):
        # per-residue mean: pLDDT of each residue (taken from its CA atom) averaged over the region
        v = [r["CA"].bfactor for r in ch if rng[0] <= r.id[1] <= rng[1] and "CA" in r]
        return np.mean(v) if v else float("nan")      # nan = region absent (e.g. trimmed model)
    print("\n5. MEAN pLDDT (B-factor column)")
    print("   PD-1 : " + ", ".join(f"{k} {mbf(X, PD1[k]):.0f}" for k in ("signal", "ecto", "stalk", "tm", "cyto")))
    print("   PD-L1: " + ", ".join(f"{k} {mbf(Y, PDL1[k]):.0f}" for k in ("signal", "ecto", "ecto2", "stalk", "tm", "cyto")))
    for nm, ch, cfg in (("PD-1", X, PD1), ("PD-L1", Y, PDL1)):
        v = [r["CA"].bfactor for r in ch if cfg["stalk"][0] <= r.id[1] <= cfg["stalk"][1]]
        print(f"   {nm} stalk CA pLDDT min/mean/max {min(v):.0f}/{np.mean(v):.0f}/{max(v):.0f}")

    # ---- 6. membrane-plane crossing
    yX = np.array([(r["CA"].coord - pt["X"]) @ n1 for r in Y if "CA" in r])
    xY = np.array([(r["CA"].coord - pt["Y"]) @ n2 for r in X if "CA" in r])
    print(f"\n6. PD-L1 CA heights above the PD-1 membrane midplane: {yX.min() / 10:.1f} to {yX.max() / 10:.1f} nm; "
          f"PD-1 CA heights above the PD-L1 midplane: {xY.min() / 10:.1f} to {xY.max() / 10:.1f} nm "
          "(positive throughout = no crossing)")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "phase4_grafted_model.pdb")
