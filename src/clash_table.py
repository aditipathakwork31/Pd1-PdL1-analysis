#!/usr/bin/env python3
"""Per-atom-pair clash table for the crystal-core hybrid (chain P = PD-1 vs chain L = PD-L1).

Reads outputs/hybrid_crystal_core.pdb and outputs/hybrid_residue_table.csv.
Writes outputs/clash_pairs.csv (every cross-chain heavy-atom pair closer than 3.0 A) and outputs/clash_summary.txt.

Definitions (state them in the Methods):
  * distance flags: lt2.0, lt2.5, lt3.0 (plain heavy-atom distance).
  * vdw_flag: approximate heavy-atom-only overlap = (r_i + r_j) - d, minus 0.6 A allowance if both atoms are N/O
    (possible H-bond); flagged when >= 0.4 A.  This is a proxy, NOT MolProbity (no hydrogens): cross-check with
    ChimeraX `clashes` / MolProbity.
  * pair_source: crystal-crystal pairs are real native contacts (reference); any pair involving an AF-derived
    residue is 'uncertain'.  Raw pLDDT and PAE-to-anchor are given so any threshold can be applied.
"""
import argparse, csv, pathlib
import numpy as np
from Bio.PDB import PDBParser, NeighborSearch

ROOT = pathlib.Path(__file__).resolve().parents[1]
RADII = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}


def element(a):
    e = (a.element or "").strip().upper()
    return e or a.get_id().strip()[0].upper()


def vdw_overlap(a, b, d):
    ea, eb = element(a), element(b)
    ov = RADII.get(ea, 1.7) + RADII.get(eb, 1.7) - d
    if ea in "NO" and eb in "NO":
        ov -= 0.6
    return ov


def find_pairs(chain_a, chain_b, cutoff=3.0):
    atoms_b = [x for r in chain_b if r.id[0] == " " for x in r if element(x) != "H"]
    ns = NeighborSearch(atoms_b)
    out = []
    for r in chain_a:
        if r.id[0] != " ":
            continue
        for a in r:
            if element(a) == "H":
                continue
            for b in ns.search(a.coord, cutoff, "A"):
                d = float(np.linalg.norm(a.coord - b.coord))
                out.append((r, a, b.get_parent(), b, d))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hybrid", default=str(ROOT / "outputs/hybrid_crystal_core.pdb"))
    ap.add_argument("--table", default=str(ROOT / "outputs/hybrid_residue_table.csv"))
    ap.add_argument("--out", default=str(ROOT / "outputs"))
    ap.add_argument("--exclude-pd1", default="85-92", help="PD-1 residue range (AF-placed loop) reported separately")
    a = ap.parse_args()
    xlo, xhi = (int(v) for v in a.exclude_pd1.split("-"))
    meta = {(r["chain"], int(r["resseq"])): r for r in csv.DictReader(open(a.table))}
    m = PDBParser(QUIET=True).get_structure("h", a.hybrid)[0]
    pairs = find_pairs(m["P"], m["L"], 3.0)
    rows = []
    for r1, a1, r2, a2, d in pairs:
        m1, m2 = meta[("P", r1.id[1])], meta[("L", r2.id[1])]
        s1, s2 = m1["source"], m2["source"]
        psrc = "crystal-crystal" if s1 == s2 == "crystal" else ("AF-AF" if s1 == s2 == "AF" else "crystal-AF")
        ov = vdw_overlap(a1, a2, d)
        rows.append(dict(
            pd1_res=r1.id[1], pd1_resname=r1.get_resname(), pd1_atom=a1.get_id(), pd1_region=m1["region"],
            pd1_source=s1, pd1_plddt=m1["af_plddt"], pd1_pae_to_anchor=m1["pae_to_anchor"],
            pdl1_res=r2.id[1], pdl1_resname=r2.get_resname(), pdl1_atom=a2.get_id(), pdl1_region=m2["region"],
            pdl1_source=s2, pdl1_plddt=m2["af_plddt"], pdl1_pae_to_anchor=m2["pae_to_anchor"],
            distance=round(d, 2), vdw_overlap=round(ov, 2), vdw_flag=int(ov >= 0.4),
            lt2_0=int(d < 2.0), lt2_5=int(d < 2.5), lt3_0=int(d < 3.0), pair_source=psrc,
            confidence="reliable (native crystal contact)" if psrc == "crystal-crystal" else "uncertain (AF-derived residue)"))
    rows.sort(key=lambda r: r["distance"])
    outdir = pathlib.Path(a.out); outdir.mkdir(parents=True, exist_ok=True)
    if rows:
        with open(outdir / "clash_pairs.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n"); w.writeheader(); w.writerows(rows)
    lines = [f"{len(rows)} cross-chain heavy-atom pairs closer than 3.0 A (PD-1 chain P vs PD-L1 chain L)", ""]
    lines.append(f"{'pair source':18}{'<2.0':>7}{'<2.5':>7}{'<3.0':>7}{'vdw_flag':>10}   (atom pairs)")
    for src in ("crystal-crystal", "crystal-AF", "AF-AF"):
        sub = [r for r in rows if r["pair_source"] == src]
        lines.append(f"{src:18}{sum(r['lt2_0'] for r in sub):>7}{sum(r['lt2_5'] for r in sub):>7}"
                     f"{len(sub):>7}{sum(r['vdw_flag'] for r in sub):>10}")
    keep = [r for r in rows if not (xlo <= r["pd1_res"] <= xhi)]
    lines += ["", f"Same table EXCLUDING PD-1 residues {xlo}-{xhi} (loop unresolved in the crystal):",
              f"{'pair source':18}{'<2.0':>7}{'<2.5':>7}{'<3.0':>7}{'vdw_flag':>10}"]
    for src in ("crystal-crystal", "crystal-AF", "AF-AF"):
        sub = [r for r in keep if r["pair_source"] == src]
        lines.append(f"{src:18}{sum(r['lt2_0'] for r in sub):>7}{sum(r['lt2_5'] for r in sub):>7}"
                     f"{len(sub):>7}{sum(r['vdw_flag'] for r in sub):>10}")
    lines += ["", "Residue pairs with any atom pair < 2.0 A:"]
    seen = {}
    for r in rows:
        if r["lt2_0"]:
            k = (r["pd1_res"], r["pdl1_res"]); seen.setdefault(k, (r["distance"], r["pair_source"], r["pd1_region"], r["pdl1_region"]))
    for (i, j), (d, ps, g1, g2) in sorted(seen.items()):
        lines.append(f"  PD-1 {i} ({g1}) - PD-L1 {j} ({g2}): {d} A  [{ps}]")
    if not seen:
        lines.append("  none")
    txt = "\n".join(lines)
    (outdir / "clash_summary.txt").write_text(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
