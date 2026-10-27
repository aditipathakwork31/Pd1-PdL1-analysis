#!/usr/bin/env python3
"""Confidence summary by protein / region / source from outputs/hybrid_residue_table.csv.
pLDDT = AlphaFold local confidence.  PAE-to-anchor = mean expected position error of the residue when the AF model
is aligned on the crystal-anchored domain (AF convention: pae[i][j] = error at j aligned on i), i.e. how well
AlphaFold itself places this residue RELATIVE to the anchored domain.  Writes outputs/region_confidence.csv/.txt"""
import argparse, csv, pathlib
from collections import defaultdict
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]


def summarize(rows):
    g = defaultdict(list)
    for r in rows:
        g[(r["protein"], r["region"], r["source"])].append(r)
    out = []
    for (p, reg, src), rs in sorted(g.items()):
        pl = np.array([float(r["af_plddt"]) for r in rs]); pae = np.array([float(r["pae_to_anchor"]) for r in rs])
        out.append(dict(protein=p, region=reg, source=src, n=len(rs), first=min(int(r["resseq"]) for r in rs),
                        last=max(int(r["resseq"]) for r in rs), plddt_mean=round(pl.mean(), 1), plddt_min=round(pl.min(), 1),
                        frac_plddt_lt70=round(float((pl < 70).mean()), 2), pae_mean=round(pae.mean(), 1),
                        pae_max=round(pae.max(), 1)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", default=str(ROOT / "outputs/hybrid_residue_table.csv"))
    ap.add_argument("--out", default=str(ROOT / "outputs"))
    a = ap.parse_args()
    out = summarize(list(csv.DictReader(open(a.table))))
    od = pathlib.Path(a.out); od.mkdir(parents=True, exist_ok=True)
    with open(od / "region_confidence.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()), lineterminator="\n"); w.writeheader(); w.writerows(out)
    L = [f"{'protein':6} {'region':7} {'source':8} {'n':>4} {'res':>9} {'pLDDT mean/min':>15} {'<70':>5} {'PAE mean/max (A)':>17}"]
    for r in out:
        L.append(f"{r['protein']:6} {r['region']:7} {r['source']:8} {r['n']:>4} {str(r['first'])+'-'+str(r['last']):>9} "
                 f"{str(r['plddt_mean'])+'/'+str(r['plddt_min']):>15} {r['frac_plddt_lt70']:>5} {str(r['pae_mean'])+'/'+str(r['pae_max']):>17}")
    txt = "\n".join(L); (od / "region_confidence.txt").write_text(txt + "\n"); print(txt)


if __name__ == "__main__":
    main()
