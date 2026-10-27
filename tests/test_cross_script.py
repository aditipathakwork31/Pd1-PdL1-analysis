"""Integration test: phase3 -> phase4 -> span_robustness vs geometry.py, on the pinned inputs.

Runs the original scripts in a temporary copy of the repo's structures/, then checks
  (1) phase4's printed spans == geometry.py spans computed from the grafted file,
  (2) span_robustness.py's printed spans == geometry.py spans,
  (3) the numbers equal the values reported in the manuscript/README.
Skipped if structures/ has not been created (run `make link`).
"""
import pathlib, re, shutil, subprocess, sys
import numpy as np
import pytest
from Bio.PDB import PDBParser

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import geometry as g

PD1 = dict(chain="X", ecto=(35, 145), tm=(171, 191))
PDL1 = dict(chain="Y", ecto=(19, 127), tm=(239, 259))

pytestmark = pytest.mark.skipif(not (ROOT / "structures" / "PD1_alphafold_full.pdb").exists(),
                                reason="structures/ missing: run `make link`")


def ca(chain, lo, hi):
    return np.array([r["CA"].coord for r in chain if lo <= r.id[1] <= hi and "CA" in r])


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    w = tmp_path_factory.mktemp("xs")
    shutil.copytree(ROOT / "structures", w / "structures")
    out = {}
    for name in ("phase3_manual_orient.py", "phase4_span_calculator.py"):
        r = subprocess.run([sys.executable, str(ROOT / "src" / name)], cwd=w, capture_output=True, text=True)
        assert r.returncode == 0, r.stdout + r.stderr
        out[name] = r.stdout
    pdb = w / "phase4_grafted_model.pdb"
    assert pdb.exists()
    r = subprocess.run([sys.executable, str(ROOT / "src" / "span_robustness.py"), str(pdb)],
                       cwd=w, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    out["span_robustness.py"] = r.stdout
    m = PDBParser(QUIET=True).get_structure("g", str(pdb))[0]
    X, Y = m["X"], m["Y"]
    tm1, tm2 = ca(X, *PD1["tm"]), ca(Y, *PDL1["tm"])
    n1 = g.tm_normal(tm1, ca(X, *PD1["ecto"])); n2 = g.tm_normal(tm2, ca(Y, *PDL1["ecto"]))
    out["geom"] = g.spans(g.centroid(tm1), n1, g.centroid(tm2), n2)
    return out


def test_phase4_matches_geometry(run):
    t = run["phase4_span_calculator.py"]
    p1 = float(re.search(r"along PD-1 normal\):\s+([\d.]+) ", t).group(1))
    p2 = float(re.search(r"along PD-L1 normal\):\s+([\d.]+) ", t).group(1))
    pd = float(re.search(r"Direct point distance:\s+([\d.]+) ", t).group(1))
    s = run["geom"]
    assert abs(abs(s["along_n1"]) - p1) < 0.06 and abs(abs(s["along_n2"]) - p2) < 0.06 and abs(s["direct"] - pd) < 0.06


def test_span_robustness_matches_geometry(run):
    t = run["span_robustness.py"]
    m = re.search(r"along PD-1 normal ([\d.]+) nm \| along PD-L1 normal ([\d.]+) nm \| direct ([\d.]+) nm \| along common normal ([\d.]+) nm", t)
    assert m, t
    a, b, d, c = (float(x) for x in m.groups())
    s = run["geom"]
    assert abs(abs(s["along_n1"]) / 10 - a) < 0.006 and abs(abs(s["along_n2"]) / 10 - b) < 0.006
    assert abs(s["direct"] / 10 - d) < 0.006 and abs(abs(s["along_common"]) / 10 - c) < 0.006


def test_signs_are_those_of_facing_membranes(run):
    assert g.facing_ok(run["geom"]), run["geom"]


def test_reproduces_reported_values(run):
    s = run["geom"]
    assert abs(s["along_n1"] / 10 - 15.22) < 0.01 and abs(abs(s["along_n2"]) / 10 - 17.39) < 0.01
    assert abs(s["direct"] / 10 - 18.10) < 0.01 and abs(s["along_common"] / 10 - 17.78) < 0.01
    assert abs(s["angle_between_normals"] - 133.0) < 0.1
