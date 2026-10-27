"""Single definition of the membrane-anchor geometry used by every script.

Conventions (all coordinates in angstrom, in the shared 4ZQK frame):
  * Membrane reference point of a protein = centroid of the C-alpha atoms of its annotated TM helix.
    ASSUMPTION: this point lies on the bilayer midplane.
  * Membrane normal = first principal axis of those TM C-alpha atoms (SVD), SIGN FIXED so it points
    from the TM centroid toward the protein's own ectodomain.
  * d = p_PDL1 - p_PD1 (from the PD-1 reference point to the PD-L1 reference point).
  * For two facing membranes: d.n1 > 0 (PD-L1 lies on the ectodomain side of PD-1's membrane) and
    d.n2 < 0, and n1.n2 = -1 (antiparallel).  Spans below are SIGNED so violations are visible.
"""
import numpy as np


def centroid(xyz):
    return np.asarray(xyz, float).mean(axis=0)


def tm_normal(tm_ca, ecto_ca):
    """Sign-fixed membrane normal (unit vector) from TM C-alpha coordinates."""
    tm_ca = np.asarray(tm_ca, float)
    c = centroid(tm_ca)
    axis = np.linalg.svd(tm_ca - c)[2][0]
    toward_ecto = centroid(ecto_ca) - c
    return axis if axis @ toward_ecto >= 0 else -axis


def angle_deg(u, v):
    return float(np.degrees(np.arccos(np.clip(u @ v / (np.linalg.norm(u) * np.linalg.norm(v)), -1, 1))))


def spans(p1, n1, p2, n2):
    """All span estimates, signed, in angstrom (divide by 10 for nm)."""
    d = np.asarray(p2, float) - np.asarray(p1, float)
    nc = (n1 - n2) / np.linalg.norm(n1 - n2)
    direct = float(np.linalg.norm(d))
    out = {
        "along_n1": float(d @ n1),
        "along_n2": float(d @ n2),
        "along_common": float(d @ nc),
        "direct": direct,
        "angle_between_normals": angle_deg(n1, n2),
    }
    out["mismatch_from_antiparallel"] = 180.0 - out["angle_between_normals"]
    out["lateral_offset_n1"] = float(np.sqrt(max(direct**2 - out["along_n1"] ** 2, 0)))
    out["lateral_offset_n2"] = float(np.sqrt(max(direct**2 - out["along_n2"] ** 2, 0)))
    return out


def facing_ok(s):
    """True if signs are consistent with two facing membranes (see module docstring)."""
    return s["along_n1"] > 0 and s["along_n2"] < 0


def tilt_range(direct, along, theta_deg):
    """Exact analytic range of |d.N| when the normal is tilted by EXACTLY theta from its
    original direction (all azimuths).  alpha = angle between d and the original normal.
    Returns (min, max) in the units of `direct`."""
    alpha = np.degrees(np.arccos(np.clip(abs(along) / direct, -1, 1)))
    lo_angle = min(alpha + theta_deg, 180.0)
    hi_angle = abs(alpha - theta_deg)
    return float(direct * np.cos(np.radians(lo_angle))), float(direct * np.cos(np.radians(hi_angle)))
