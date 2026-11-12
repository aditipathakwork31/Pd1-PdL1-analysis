#!/usr/bin/env python3
"""Feasible intermembrane span for a trans PD-1:PD-L1 complex (replaces the unweighted tilt grid).

MODEL (all assumptions are explicit; see docs/geometry_note.md):
  * The ectodomain complex is rigid: PD-1 (chain P, 25-145 minus the AF loop) + PD-L1 (chain L, 19-225), 4ZQK frame.
  * The two membranes are PARALLEL planes with common unit normal N (N points from the PD-1 membrane to the PD-L1 one).
    Heights are z = x.N.  Membrane midplane heights Z1 (PD-1 cell) < Z2 (PD-L1 cell); surfaces are delta from the midplane.
  * Non-penetration: every ectodomain atom lies between the two surfaces:  Z1 <= zmin - delta,  Z2 >= zmax + delta.
  * Stalk reach: the ectodomain C-terminus (PD-1 CA145, PD-L1 CA225) must be within the stalk contour length of its
    membrane surface (vertical reach; the TM anchor can move laterally):
        z(CA145) - (Z1 + delta) <= L1,   (Z2 - delta) - z(CA225) <= L2,   L = (number of peptide steps) * step.
  * span = Z2 - Z1 (midplane to midplane), so for each N:
        S_min = (zmax - zmin) + 2 delta      S_max = z(CA225) - z(CA145) + L1 + L2 + 2 delta
    and N is feasible only if both stalks can reach: z(CA145) - zmin <= L1 and zmax - z(CA225) <= L2.
  * Orientation prior: UNIFORM over directions N (Fibonacci sphere).  AlphaFold gives no information on the TM
    orientation (PAE ~25-28 A), so no orientation is preferred.  Tilt of the TM helices relative to N is not constrained
    because the stalks (25 and 13 residues) can reorient them.
  * Optional IgC2 hinge uncertainty: PD-L1 133-225 rotated about CA132 by a random angle <= hinge_deg.
Outputs: outputs/feasible_span.csv and outputs/feasible_span.txt"""
import argparse, csv, pathlib, sys
import numpy as np
from scipy.spatial import ConvexHull
from scipy.spatial.transform import Rotation
from Bio.PDB import PDBParser

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import geometry as g

STEPS = {"P": 26, "L": 14}      # peptide steps from last Ig residue (145 / 225) to first TM residue (171 / 239)


def fib_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n); th = np.pi * (1 + 5 ** 0.5) * i
    return np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]


def scan(N, pts, c1, c2, L1, L2, delta):
    """N: (k,3) unit directions; pts: ectodomain atom coordinates.  Returns feasible mask, S_min, S_max (A)."""
    hull = pts[ConvexHull(pts).vertices]
    z = hull @ N.T
    zmin, zmax = z.min(0), z.max(0)
    z1, z2 = c1 @ N.T, c2 @ N.T
    feas = (z1 - zmin <= L1) & (zmax - z2 <= L2)
    return feas, zmax - zmin + 2 * delta, z2 - z1 + L1 + L2 + 2 * delta


def atoms(chain, ranges, skip=()):
    out = []
    for r in chain:
        if r.id[0] != " ":
            continue
        n = r.id[1]
        if any(lo <= n <= hi for lo, hi in ranges) and not any(lo <= n <= hi for lo, hi in skip):
            out += [a.coord for a in r if (a.element or "").strip().upper() != "H"]
    return np.array(out)


def _randvec(rng, max_deg):
    axis = rng.normal(size=3); axis /= np.linalg.norm(axis)
    return axis * np.radians(rng.uniform(0, max_deg))


def run(model, deltas=(15, 20, 25), steps=(3.3, 3.5, 3.8), hinges=(0, 20, 40), ndir=10000, nhinge=60,
        skip=((85, 92),), seed=1):
    rng = np.random.default_rng(seed)
    N = fib_sphere(ndir)
    p_ecto = atoms(model["P"], [(25, 145)], skip)
    l_igv = atoms(model["L"], [(19, 132)])
    l_igc = atoms(model["L"], [(133, 225)])
    ca = lambda ch, n: next(r["CA"].coord for r in model[ch] if r.id[0] == " " and r.id[1] == n)
    c1, c2, pivot = ca("P", 145), ca("L", 225), ca("L", 132)
    rows = []
    for delta in deltas:
        for step in steps:
            for hinge in hinges:
                F, SMIN, SMAX = [], [], []
                for k in range(nhinge if hinge else 1):
                    if hinge:
                        R = Rotation.from_rotvec(_randvec(rng, hinge))
                        igc = (l_igc - pivot) @ R.as_matrix().T + pivot
                        c2k = (c2 - pivot) @ R.as_matrix().T + pivot
                    else:
                        igc, c2k = l_igc, c2
                    f, a, b = scan(N, np.vstack([p_ecto, l_igv, igc]), c1, c2k, STEPS["P"] * step, STEPS["L"] * step, delta)
                    F.append(f); SMIN.append(a); SMAX.append(b)
                f, a, b = np.concatenate(F), np.concatenate(SMIN), np.concatenate(SMAX)
                rows.append(dict(delta_nm=delta / 10, step_A=step, hinge_deg=hinge, frac_directions_feasible=round(float(f.mean()), 3),
                                 Smin_lowest_nm=round(float(a[f].min()) / 10, 1) if f.any() else None,
                                 Smin_median_nm=round(float(np.median(a[f])) / 10, 1) if f.any() else None,
                                 Smax_median_nm=round(float(np.median(b[f])) / 10, 1) if f.any() else None,
                                 Smax_highest_nm=round(float(b[f].max()) / 10, 1) if f.any() else None))
    return rows, (p_ecto, l_igv, l_igc, c1, c2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hybrid", default=str(ROOT / "outputs/hybrid_crystal_core.pdb"))
    ap.add_argument("--out", default=str(ROOT / "outputs"))
    a = ap.parse_args()
    model = PDBParser(QUIET=True).get_structure("h", a.hybrid)[0]
    rows, (p_ecto, l_igv, l_igc, c1, c2) = run(model)
    od = pathlib.Path(a.out); od.mkdir(parents=True, exist_ok=True)
    with open(od / "feasible_span.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n"); w.writeheader(); w.writerows(rows)
    L = ["Feasible midplane-to-midplane span (nm); uniform prior over membrane-normal directions",
         f"{'delta':>6}{'step A':>8}{'hinge':>7}{'feasible':>10}{'Smin low':>10}{'Smin med':>10}{'Smax med':>10}{'Smax high':>11}"]
    for r in rows:
        L.append(f"{r['delta_nm']:>6}{r['step_A']:>8}{r['hinge_deg']:>7}{r['frac_directions_feasible']:>10}{r['Smin_lowest_nm']!s:>10}"
                 f"{r['Smin_median_nm']!s:>10}{r['Smax_median_nm']!s:>10}{r['Smax_highest_nm']!s:>11}")
    def cap(ch, lo, hi): return np.array([x["CA"].coord for x in model[ch] if x.id[0] == " " and lo <= x.id[1] <= hi])
    n1 = g.tm_normal(cap("P", 171, 191), cap("P", 35, 145)); n2 = g.tm_normal(cap("L", 239, 259), cap("L", 19, 127))
    nc = (n1 - n2) / np.linalg.norm(n1 - n2)
    L += ["", "At the three AF-derived directions (delta 2.0 nm, 3.5 A/residue, rigid):"]
    for name, n in (("PD-1 normal n1", n1), ("-PD-L1 normal", -n2), ("common normal", nc)):
        f, lo, hi = scan(n[None, :], np.vstack([p_ecto, l_igv, l_igc]), c1, c2, STEPS["P"] * 3.5, STEPS["L"] * 3.5, 20)
        L.append(f"  {name:16} feasible={bool(f[0])}  span range {lo[0] / 10:.1f}-{hi[0] / 10:.1f} nm")
    txt = "\n".join(L); (od / "feasible_span.txt").write_text(txt + "\n"); print(txt)


if __name__ == "__main__":
    main()
