#!/usr/bin/env python3
"""Build the crystal-core hybrid (replaces the AlphaFold-everywhere hybrid for clash analysis).

Reads the grafted model written by phase4 (chains A,B = 4ZQK crystal PD-L1/PD-1; chains X,Y = AlphaFold
PD-1/PD-L1 superimposed onto the crystal; all in the 4ZQK frame) and assembles, with UniProt numbering:

  PD-1  (new chain P): AF 25-32 | crystal B 33-84 | AF 85-92 (loop unresolved in crystal) | crystal B 93-146 | AF 147-191
  PD-L1 (new chain L): crystal A 18-132 | AF 133-259

Signal peptides (cleaved) and cytoplasmic tails are NOT included.  Waters/ions are dropped.
The hybrid is a rigid-body assembly of pieces, NOT a covalently continuous or energy-minimised model:
the C-N distance at every seam is reported in hybrid_seams.csv.

Outputs (in --out):  hybrid_crystal_core.pdb, hybrid_residue_table.csv, hybrid_seams.csv
"""
import argparse, csv, json, pathlib
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain

ROOT = pathlib.Path(__file__).resolve().parents[1]

REGIONS = {  # UniProt-based regions (Table 2 of the manuscript), inclusive ranges
    "PD1": [(1, 24, "signal"), (25, 34, "linker"), (35, 145, "ecto"), (146, 170, "stalk"),
            (171, 191, "tm"), (192, 288, "cyto")],
    "PDL1": [(1, 18, "signal"), (19, 127, "ecto"), (128, 132, "linker"), (133, 225, "ecto2"),
             (226, 238, "stalk"), (239, 259, "tm"), (260, 290, "cyto")],
}
PIECES = {  # (source chain in grafted file, first, last, source label, note)
    "PD1": [("X", 25, 32, "AF", ""), ("B", 33, 84, "crystal", ""),
            ("X", 85, 92, "AF", "unresolved in crystal; AF-placed loop"),
            ("B", 93, 146, "crystal", "residue 93 is Ser in crystal (C93S), Cys in AF"), ("X", 147, 191, "AF", "")],
    "PDL1": [("A", 18, 132, "crystal", ""), ("Y", 133, 259, "AF", "")],
}
NEW_CHAIN = {"PD1": "P", "PDL1": "L"}
AF_CHAIN = {"PD1": "X", "PDL1": "Y"}
ANCHOR = {"PD1": (35, 145), "PDL1": (19, 127)}   # domain whose crystal structure fixes the frame


def region(protein, n):
    for lo, hi, name in REGIONS[protein]:
        if lo <= n <= hi:
            return name
    return "?"


def get_res(chain, n):
    for r in chain:
        if r.id[0] == " " and r.id[1] == n:
            return r
    raise KeyError(f"residue {n} missing in chain {chain.id}")


def pae_to_anchor(pae, protein, n):
    """Mean expected position error of residue n when the AF model is aligned on the anchor domain.
    AF convention: pae[i][j] = error at residue j when aligned on residue i (1-based residue = index+1)."""
    lo, hi = ANCHOR[protein]
    return float(np.mean(pae[lo - 1:hi, n - 1]))


def load_pae(path):
    d = json.load(open(path))
    d = d[0] if isinstance(d, list) else d
    return np.array(d["predicted_aligned_error"], float)


def build(grafted, pae_paths, out):
    out = pathlib.Path(out); out.mkdir(parents=True, exist_ok=True)
    model = PDBParser(QUIET=True).get_structure("g", str(grafted))[0]
    paes = {p: load_pae(pae_paths[p]) for p in PIECES}
    st = Structure.Structure("hybrid"); md = Model.Model(0); st.add(md)
    rows, seams, mism = [], [], []
    for prot, pieces in PIECES.items():
        ch = Chain.Chain(NEW_CHAIN[prot]); md.add(ch)
        af = model[AF_CHAIN[prot]]
        for src, lo, hi, label, note in pieces:
            for n in range(lo, hi + 1):
                r = get_res(model[src], n)
                afr = get_res(af, n)
                if label == "crystal" and r.get_resname() != afr.get_resname():
                    mism.append((prot, n, r.get_resname(), afr.get_resname()))
                ch.add(r.copy())
                rows.append(dict(protein=prot, chain=NEW_CHAIN[prot], resseq=n, resname=r.get_resname(),
                                 source=label, note=note, region=region(prot, n),
                                 af_plddt=round(afr["CA"].get_bfactor(), 1),
                                 pae_to_anchor=round(pae_to_anchor(paes[prot], prot, n), 2),
                                 crystal_resolved=int(label == "crystal")))
        res = [r for r in ch]
        for a, b in zip(res, res[1:]):
            d = float(np.linalg.norm(a["C"].coord - b["N"].coord))
            sa = next(x["source"] for x in rows if x["chain"] == ch.id and x["resseq"] == a.id[1])
            sb = next(x["source"] for x in rows if x["chain"] == ch.id and x["resseq"] == b.id[1])
            gap = b.id[1] != a.id[1] + 1
            if sa != sb or d > 1.6 or gap:
                seams.append(dict(chain=ch.id, res_i=a.id[1], res_j=b.id[1], source_i=sa, source_j=sb,
                                  c_n_distance=round(d, 2), bonded=int(d <= 1.6 and not gap)))
    io = PDBIO(); io.set_structure(st); io.save(str(out / "hybrid_crystal_core.pdb"))
    for name, data in (("hybrid_residue_table.csv", rows), ("hybrid_seams.csv", seams)):
        with open(out / name, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(data[0].keys()), lineterminator="\n"); w.writeheader(); w.writerows(data)
    return rows, seams, mism


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grafted", default=str(ROOT / "outputs/phase4_grafted_model.pdb"))
    ap.add_argument("--pae-pd1", default=str(ROOT / "data/raw/alphafold/AF-Q15116-F1-predicted_aligned_error_v6.json"))
    ap.add_argument("--pae-pdl1", default=str(ROOT / "data/raw/alphafold/AF-Q9NZQ7-F1-predicted_aligned_error_v6.json"))
    ap.add_argument("--out", default=str(ROOT / "outputs"))
    a = ap.parse_args()
    rows, seams, mism = build(a.grafted, {"PD1": a.pae_pd1, "PDL1": a.pae_pdl1}, a.out)
    n1 = sum(r["protein"] == "PD1" for r in rows); n2 = len(rows) - n1
    print(f"PD-1 chain P: {n1} residues; PD-L1 chain L: {n2} residues")
    print("crystal vs AF residue-name mismatches at the same UniProt number (expect only PD-1 93):", mism)
    print("seams / breaks (C-N distance, A):")
    for s in seams:
        print(f"  chain {s['chain']} {s['res_i']}->{s['res_j']} {s['source_i']}->{s['source_j']}  {s['c_n_distance']}  bonded={s['bonded']}")


if __name__ == "__main__":
    main()
