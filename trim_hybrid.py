#!/usr/bin/env python3
"""
Trim the AlphaFold chains of the grafted hybrid model: remove signal peptides (cleaved in the mature proteins) and
cytoplasmic tails (intracellular, unconstrained in the monomer prediction). Coordinates of the retained residues are
unchanged, so the TM midpoints, membrane normals and every span value are unaffected; only the clash count changes.

Usage:  python trim_hybrid.py phase4_grafted_model.pdb phase4_grafted_model_trimmed.pdb
Keeps PD-1 (chain X) residues 25-191 and PD-L1 (chain Y) residues 19-259; crystal chains A and B are copied untouched.
"""
import sys
from Bio.PDB import PDBParser, PDBIO, Select

KEEP = {"X": (25, 191), "Y": (19, 259)}


class Trim(Select):
    def accept_residue(self, res):
        ch = res.get_parent().id
        if ch not in KEEP:
            return True
        lo, hi = KEEP[ch]
        return lo <= res.id[1] <= hi


def main(src, dst):
    s = PDBParser(QUIET=True).get_structure("g", src)
    io = PDBIO(); io.set_structure(s); io.save(dst, Trim())
    t = PDBParser(QUIET=True).get_structure("t", dst)[0]
    for c in t:
        r = [x for x in c if x.id[0] == " "]
        print(f"chain {c.id}: {len(r)} residues ({r[0].id[1]}-{r[-1].id[1]})")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
