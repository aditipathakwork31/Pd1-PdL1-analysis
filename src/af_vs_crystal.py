#!/usr/bin/env python3
"""Where does the AlphaFold monomer disagree with the crystal, and is it a rigid-fit (alignment) problem?

For every residue resolved in the 4ZQK crystal (same residue name in AF), after the GLOBAL fit used by phase4:
  ca_dev (A), backbone RMSD (N,CA,C,O), side-chain RMSD (common heavy atoms).
Then re-fit AF onto the crystal using ONLY the Calpha atoms of interface residues (any atom within 4.5 A of the
partner in the crystal) and recompute.  If deviations at the clashing residues persist after the local fit, the
disagreement is real (conformational / side-chain), not an alignment artefact.

Reads outputs/phase4_grafted_model.pdb (chains A,B crystal; X,Y AlphaFold).  Writes outputs/af_vs_crystal.csv/.txt
"""
import argparse, csv, pathlib
import numpy as np
from Bio.PDB import PDBParser, NeighborSearch, Superimposer

ROOT = pathlib.Path(__file__).resolve().parents[1]
BB = ("N", "CA", "C", "O")
PAIRS = {"PD1": ("B", "X", "A"), "PDL1": ("A", "Y", "B")}   # crystal chain, AF chain, partner crystal chain


def get(chain, n):
    for r in chain:
        if r.id[0] == " " and r.id[1] == n:
            return r
    return None


def matched(crys, af):
    out = []
    for r in crys:
        if r.id[0] != " ":
            continue
        a = get(af, r.id[1])
        if a is not None and a.get_resname() == r.get_resname() and "CA" in r and "CA" in a:
            out.append((r, a))
    return out


def interface(chain, partner, cutoff=4.5):
    ns = NeighborSearch([x for r in partner if r.id[0] == " " for x in r])
    return {r.id[1] for r in chain if r.id[0] == " " and any(ns.search(x.coord, cutoff) for x in r)}


def rms(a, b):
    return float(np.sqrt(np.mean(np.sum((np.array(a) - np.array(b)) ** 2, axis=1)))) if len(a) else float("nan")


def residue_stats(r, a, rot=None, tran=None):
    def xyz(atom):
        c = atom.coord
        return c @ rot + tran if rot is not None else c
    ca = float(np.linalg.norm(r["CA"].coord - xyz(a["CA"])))
    bb = [n for n in BB if n in r and n in a]
    sc = [x.get_id() for x in r if x.get_id() not in BB and x.get_id() in a and not x.get_id().startswith("H")]
    return ca, rms([r[n].coord for n in bb], [xyz(a[n]) for n in bb]), \
        rms([r[n].coord for n in sc], [xyz(a[n]) for n in sc])


def analyse(model, focus):
    rows = []
    for prot, (cc, ac, pc) in PAIRS.items():
        mt = matched(model[cc], model[ac])
        iface = interface(model[cc], model[pc])
        fit = [(r, a) for r, a in mt if r.id[1] in iface]
        sup = Superimposer(); sup.set_atoms([r["CA"] for r, _ in fit], [a["CA"] for _, a in fit])
        rot, tran = sup.rotran
        for r, a in mt:
            g = residue_stats(r, a); l = residue_stats(r, a, rot, tran)
            rows.append(dict(protein=prot, resseq=r.id[1], resname=r.get_resname(), interface=int(r.id[1] in iface),
                             focus=int(r.id[1] in focus.get(prot, ())), plddt=round(a["CA"].get_bfactor(), 1),
                             ca_global=round(g[0], 2), bb_global=round(g[1], 2), sc_global=round(g[2], 2),
                             ca_local=round(l[0], 2), bb_local=round(l[1], 2), sc_local=round(l[2], 2)))
        rows.append(dict(protein=prot, resseq=-1, resname="LOCALFIT", interface=len(fit), focus=0, plddt=0,
                         ca_global=0, bb_global=0, sc_global=0, ca_local=round(sup.rms, 2), bb_local=0, sc_local=0))
    return rows


def summarize(rows):
    L = []
    for prot in PAIRS:
        rr = [r for r in rows if r["protein"] == prot and r["resseq"] > 0]
        lf = next(r for r in rows if r["protein"] == prot and r["resseq"] == -1)
        L.append(f"{prot}: {len(rr)} crystal-resolved residues; local fit on {lf['interface']} interface C-alpha, fit RMSD {lf['ca_local']} A")
        for name, sub in (("all", rr), ("interface", [r for r in rr if r["interface"]])):
            for k in ("ca_global", "ca_local", "sc_global", "sc_local"):
                v = np.array([r[k] for r in sub if not np.isnan(r[k])])
                if v.size == 0:
                    L.append(f"   {name:9} {k:10} n/a"); continue
                L.append(f"   {name:9} {k:10} median {np.median(v):.2f}  90th pct {np.percentile(v, 90):.2f}  max {v.max():.2f}")
        L.append("   focus residues (those in the original 21 clashing pairs):")
        L.append("   res name  pLDDT  CA_glob CA_loc  SC_glob SC_loc")
        for r in rr:
            if r["focus"]:
                L.append(f"   {r['resseq']:<4}{r['resname']:5}{r['plddt']:6}{r['ca_global']:8}{r['ca_local']:7}{r['sc_global']:8}{r['sc_local']:7}")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grafted", default=str(ROOT / "outputs/phase4_grafted_model.pdb"))
    ap.add_argument("--out", default=str(ROOT / "outputs"))
    ap.add_argument("--focus-pd1", default="88,131,132")
    ap.add_argument("--focus-pdl1", default="19,56,63,66")
    a = ap.parse_args()
    focus = {"PD1": {int(x) for x in a.focus_pd1.split(",")}, "PDL1": {int(x) for x in a.focus_pdl1.split(",")}}
    model = PDBParser(QUIET=True).get_structure("g", a.grafted)[0]
    rows = analyse(model, focus)
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    with open(out / "af_vs_crystal.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n"); w.writeheader(); w.writerows(rows)
    txt = summarize(rows)
    (out / "af_vs_crystal.txt").write_text(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
