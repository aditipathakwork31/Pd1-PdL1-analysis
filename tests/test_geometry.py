import sys, pathlib
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
import geometry as g


def helix(center, axis, n=21, rise=1.5, radius=2.3):
    """Ideal alpha-helix C-alpha trace along `axis` through `center`."""
    axis = np.asarray(axis, float); axis /= np.linalg.norm(axis)
    t = np.cross(axis, [1, 0, 0]) if abs(axis[0]) < 0.9 else np.cross(axis, [0, 1, 0])
    t /= np.linalg.norm(t); u = np.cross(axis, t)
    k = np.arange(n) - (n - 1) / 2
    ang = np.radians(100) * k
    pts = np.array(center, float) + np.outer(k * rise, axis) + radius * (np.outer(np.cos(ang), t) + np.outer(np.sin(ang), u))
    return pts


def test_hand_checkable_100A_is_10nm():
    # PD-1 TM at origin, ectodomain above; PD-L1 TM at z=100, its ectodomain below it
    tm1, tm2 = helix([0, 0, 0], [0, 0, 1]), helix([0, 0, 100], [0, 0, 1])
    n1 = g.tm_normal(tm1, [[0, 0, 50]])
    n2 = g.tm_normal(tm2, [[0, 0, 50]])
    s = g.spans(g.centroid(tm1), n1, g.centroid(tm2), n2)   # reference points are TM centroids by definition
    assert abs(s["along_n1"] - 100) < 0.5 and abs(s["direct"] - 100) < 1e-9
    assert abs(s["along_n2"] + 100) < 0.5
    assert abs(s["along_common"] - 100) < 0.5
    # The SVD axis of a 21-residue C-alpha trace is not exactly the ideal helix axis (non-integer number
    # of turns); it is off by ~1 degree here, i.e. ~2.5 A of apparent lateral offset over 100 A.
    # This is a real, small uncertainty in the normal that the geometry note should quote.
    assert s["mismatch_from_antiparallel"] < 0.5 and s["lateral_offset_n1"] < 5.0
    assert g.facing_ok(s)


def test_sign_independent_of_point_order():
    tm = helix([3, -2, 7], [1, 2, 3]); ecto = [[10, 20, 30]]
    a = g.tm_normal(tm, ecto); b = g.tm_normal(tm[::-1], ecto)
    assert np.allclose(a, b, atol=1e-9)
    assert a @ (np.mean(ecto, 0) - tm.mean(0)) > 0


def test_wrong_side_is_detected():
    n1 = np.array([0, 0, 1.0]); n2 = np.array([0, 0, -1.0])
    s = g.spans([0, 0, 0], n1, [0, 0, -100], n2)   # PD-L1 membrane BELOW PD-1 membrane
    assert not g.facing_ok(s)


def test_lateral_offset():
    s = g.spans([0, 0, 0], np.array([0, 0, 1.0]), [30, 0, 100], np.array([0, 0, -1.0]))
    assert abs(s["lateral_offset_n1"] - 30) < 1e-9


def test_tilt_zero_alpha():
    lo, hi = g.tilt_range(100.0, 100.0, 30)
    assert abs(lo - 86.6025) < 1e-3 and abs(hi - 86.6025) < 1e-3     # 100*cos(30 deg)


def test_reproduces_manuscript_table4_within_grid_error():
    # Table 4 (nm): numbers produced by the 12-azimuth grid; analytic exact values must agree within 0.2 nm
    d = 18.10
    table = {("pd1", 10): (13.3, 16.6), ("pd1", 20): (11.0, 17.6), ("pd1", 30): (8.4, 18.0),
             ("pdl1", 10): (16.3, 18.0), ("pdl1", 20): (14.7, 18.0), ("pdl1", 30): (12.6, 17.5)}
    along = {"pd1": 15.22, "pdl1": 17.39}
    for (k, th), (lo_t, hi_t) in table.items():
        lo, hi = g.tilt_range(d, along[k], th)
        assert abs(lo - lo_t) < 0.2 and abs(hi - hi_t) < 0.2, (k, th, lo, hi)
