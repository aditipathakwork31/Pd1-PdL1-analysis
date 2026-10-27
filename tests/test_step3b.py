import csv, pathlib, sys
import numpy as np
from Bio.PDB import Structure, Model, Chain, Residue, Atom
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import region_confidence as rc, membrane_check as mc, glycan_screen as gs


def test_region_confidence_groups():
    rows = [dict(protein="PD1", region="stalk", source="AF", resseq=r, af_plddt=p, pae_to_anchor=a)
            for r, p, a in ((146, 60, 10.0), (147, 80, 20.0))]
    s = rc.summarize(rows)[0]
    assert s["n"] == 2 and s["plddt_mean"] == 70.0 and s["frac_plddt_lt70"] == 0.5 and s["pae_max"] == 20.0


def test_min_height():
    xyz = [[0, 0, 30.0], [5, 5, 45.0], [1, 1, 20.5]]
    assert abs(mc.min_height(xyz, np.array([0, 0, 0.0]), np.array([0, 0, 1.0])) - 20.5) < 1e-9
    assert abs(mc.min_height(xyz, np.array([0, 0, 10.0]), np.array([0, 0, 1.0])) - 10.5) < 1e-9


def test_glycan_screen_distance_and_flag():
    m = Model.Model(0)
    for cid, n, x in (("P", 49, 0.0), ("L", 35, 12.0)):
        ch = Chain.Chain(cid); r = Residue.Residue((" ", n, " "), "ASN", " ")
        for nm, el, p in (("CB", "C", (x, 0, 0)), ("ND2", "N", (x, 1, 0))):
            r.add(Atom.Atom(nm, np.array(p, float), 1, 1, " ", nm.ljust(4), 1, el))
        ch.add(r); m.add(ch)
    rows = gs.screen(m)
    p49 = next(r for r in rows if r[0] == "P" and r[1] == 49)
    assert abs(p49[3] - 12.0) < 1e-9 and p49[4] == "could contact partner"
    assert next(r for r in rows if r[0] == "P" and r[1] == 58)[2] == "absent"
