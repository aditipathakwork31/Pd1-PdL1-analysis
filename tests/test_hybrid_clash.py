import csv, json, pathlib, sys
import numpy as np
import pytest
from Bio.PDB import Structure, Model, Chain, Residue, Atom, PDBIO

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import build_hybrid as bh
import clash_table as ct


def mk_res(n, x, y, name="GLY"):
    r = Residue.Residue((" ", n, " "), name, " ")
    for nm, el, dx in (("N", "N", -1.2), ("CA", "C", 0.0), ("C", "C", 1.2), ("O", "O", 1.2)):
        r.add(Atom.Atom(nm, np.array([x + dx, y + (0.5 if nm == "O" else 0), 0.0]), 50.0, 1.0, " ", nm.ljust(4), n * 4 + len(nm), el))
    return r


def mk_chain(cid, lo, hi, y):
    c = Chain.Chain(cid)
    for n in range(lo, hi + 1):
        c.add(mk_res(n, n * 3.8, y))
    return c


@pytest.fixture
def synth(tmp_path):
    s = Structure.Structure("g"); m = Model.Model(0); s.add(m)
    for cid, lo, hi, y in (("A", 18, 132, 100.0), ("B", 33, 146, 0.0), ("X", 1, 288, 0.0), ("Y", 1, 290, 100.0)):
        m.add(mk_chain(cid, lo, hi, y))
    io = PDBIO(); io.set_structure(s); io.save(str(tmp_path / "g.pdb"))
    for nm, n in (("p1", 288), ("p2", 290)):
        json.dump([{"predicted_aligned_error": np.full((n, n), 3.0).tolist(), "max_predicted_aligned_error": 31.75}],
                  open(tmp_path / f"{nm}.json", "w"))
    return tmp_path


def test_region_boundaries():
    assert bh.region("PD1", 145) == "ecto" and bh.region("PD1", 146) == "stalk"
    assert bh.region("PD1", 171) == "tm" and bh.region("PD1", 191) == "tm"
    assert bh.region("PDL1", 132) == "linker" and bh.region("PDL1", 133) == "ecto2"
    assert bh.region("PDL1", 239) == "tm"


def test_build_counts_and_sources(synth):
    rows, seams, mism = bh.build(synth / "g.pdb", {"PD1": synth / "p1.json", "PDL1": synth / "p2.json"}, synth)
    p = [r for r in rows if r["chain"] == "P"]; l = [r for r in rows if r["chain"] == "L"]
    assert len(p) == 191 - 25 + 1 and len(l) == 259 - 18 + 1
    assert [r["resseq"] for r in p] == list(range(25, 192))          # contiguous UniProt numbering
    assert sum(r["source"] == "crystal" for r in p) == (84 - 33 + 1) + (146 - 93 + 1)
    assert all(r["source"] == "AF" for r in p if 85 <= r["resseq"] <= 92)
    assert all(abs(r["pae_to_anchor"] - 3.0) < 1e-9 for r in rows)
    assert mism == []
    assert all(s["bonded"] == 1 for s in seams)                       # synthetic pieces share coordinates


def test_pair_finder_and_vdw():
    a = Chain.Chain("P"); b = Chain.Chain("L")
    ra = Residue.Residue((" ", 1, " "), "GLY", " "); rb = Residue.Residue((" ", 1, " "), "GLY", " ")
    ra.add(Atom.Atom("N", np.array([0.0, 0, 0]), 1, 1, " ", " N  ", 1, "N"))
    rb.add(Atom.Atom("O", np.array([2.0, 0, 0]), 1, 1, " ", " O  ", 2, "O"))
    rb.add(Atom.Atom("CA", np.array([0.0, 2.9, 0]), 1, 1, " ", " CA ", 3, "C"))
    a.add(ra); b.add(rb)
    pairs = ct.find_pairs(a, b, 3.0)
    assert len(pairs) == 2
    d = {p[3].get_id(): p[4] for p in pairs}
    assert abs(d["O"] - 2.0) < 1e-9 and abs(d["CA"] - 2.9) < 1e-9
    n = ra["N"]
    assert abs(ct.vdw_overlap(n, rb["O"], 2.0) - (1.55 + 1.52 - 2.0 - 0.6)) < 1e-9     # N-O, H-bond allowance: 0.47 -> flagged
    assert ct.vdw_overlap(n, rb["O"], 2.3) < 0.4                                         # 0.17 -> not flagged
    assert abs(ct.vdw_overlap(n, rb["CA"], 2.9) - (1.55 + 1.70 - 2.9)) < 1e-9           # N-C: 0.35 -> not flagged


def test_af_vs_crystal_detects_sidechain_difference():
    import af_vs_crystal as av
    def res(n, x, cb_dy):
        r = Residue.Residue((" ", n, " "), "ALA", " ")
        for nm, el, p in (("N", "N", (x - 1.2, 0, 0)), ("CA", "C", (x, 0, 0)), ("C", "C", (x + 1.2, 0, 0)),
                          ("O", "O", (x + 1.2, 1, 0)), ("CB", "C", (x, 1.5 + cb_dy, 0))):
            r.add(Atom.Atom(nm, np.array(p, float), 80.0, 1.0, " ", nm.ljust(4), 1, el))
        return r
    m = Model.Model(0)
    for cid, dy in (("B", 0.0), ("X", 0.0), ("A", 0.0), ("Y", 0.0)):
        c = Chain.Chain(cid)
        for n in range(1, 7):
            c.add(res(n, n * 3.8 + (0 if cid in "BX" else 0), dy if n != 3 else (3.0 if cid in "XY" else 0.0)))
        m.add(c)
    rows = av.analyse(m, {"PD1": {3}, "PDL1": set()})
    r3 = next(r for r in rows if r["protein"] == "PD1" and r["resseq"] == 3)
    r4 = next(r for r in rows if r["protein"] == "PD1" and r["resseq"] == 4)
    assert r3["ca_global"] < 1e-6 and r3["bb_global"] < 1e-6        # backbone identical
    assert abs(r3["sc_global"] - 3.0) < 1e-6 and r4["sc_global"] < 1e-6   # only residue 3's CB moved, by 3.0 A
