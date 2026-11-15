#!/usr/bin/env python3
"""Figure 1: modelled PD-1:PD-L1 separation against reference values.

The script that drew the manuscript's first version of this figure was not preserved; this script redraws it.
Every number is defined once below, with its source.

Usage (from the repository root):  python src/make_figure1.py
Output: outputs/figure1_span_comparison.png
"""
import pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "figure1_span_comparison.png"

# Computed here (src/span_robustness.py on the full hybrid; see expected/span_robustness.txt):
# midplane-to-midplane span at zero tilt, lowest and highest of the four definitions (nm).
MODEL_MIN, MODEL_MAX = 15.22, 18.10
# Illustration only, not computed: assumed total bilayer offset between midplane-to-midplane
# and surface-to-surface conventions (nm).
SURFACE_OFFSET = 4.0

# Reference value (secondary source, not re-derived here):
#  - Immune receptor-ligand complexes (TCR:pMHC, KIR:pMHC, NKG2D:MICA) have an extracellular span of about 10-15 nm
#    from crystal structures, against about 40 nm for integrin complexes (Alakoskela et al. 2011, Biophys J 100:2865-2874,
#    PMC3123984). That paper cites Groves 2007 (review), Fan et al. 2001 (KIR2DL1-HLA-Cw4) and Li et al. 2001
#    (NKG2D-MICA). It does not state the measurement convention (end-to-end or membrane-normal projected length).
#  - A separate TCR:pMHC value (~13 nm) and an "IgSF pairs ~15 nm" value (draft attribution "Strazza 2024" not found)
#    were removed because no primary source was verified.
RL_MIN, RL_MAX = 10.0, 15.0  # nm, range as reported

# Verified: Cordoba et al. 2013 (Blood 121:4295-4302) reports about 28 nm (smallest isoform) to about 50 nm (largest)
# for the CD45 ectodomain, from low-resolution electron microscopy (Woollett et al. 1985; McCall et al. 1992).
# The manuscript draft used ~25 nm.
CD45_SMALLEST = 28.0  # nm


def main():
    fig, ax = plt.subplots(figsize=(8, 5), dpi=200)
    lo, hi = MODEL_MIN, MODEL_MAX
    ax.bar(0, hi - lo, bottom=lo, width=0.55, color="#4c72b0")
    ax.text(0, hi + 0.4, f"{lo:.1f}-{hi:.1f}", ha="center", va="bottom", fontsize=9)

    slo, shi = lo - SURFACE_OFFSET, hi - SURFACE_OFFSET
    ax.bar(1, shi - slo, bottom=slo, width=0.55, color="#9db3d8", edgecolor="#4c72b0", hatch="/")
    ax.text(1, shi + 0.4, f"~{slo:.1f}-{shi:.1f}", ha="center", va="bottom", fontsize=9)

    ax.bar(2, RL_MAX - RL_MIN, bottom=RL_MIN, width=0.55, color="#8c8c8c")
    ax.text(2, RL_MAX + 0.4, f"~{RL_MIN:.0f}-{RL_MAX:.0f}", ha="center", va="bottom", fontsize=9)

    ax.axhline(CD45_SMALLEST, color="firebrick", linestyle="--", linewidth=1.2)
    ax.text(2.35, CD45_SMALLEST + 0.5,
            f"CD45 smallest-isoform ectodomain, ~{CD45_SMALLEST:.0f} nm (range 28-50 nm; sterically excluded, not a bound span)",
            ha="right", va="bottom", fontsize=7.5, color="firebrick")

    ax.set_xticks(range(3))
    ax.set_xticklabels(["PD-1:PD-L1\n(this work)\nmidplane-to-midplane",
                        "PD-1:PD-L1\n(this work)\napprox. surface-to-surface*",
                        "Immune receptor-ligand\ncomplexes, crystal structures\n(Alakoskela et al. 2011)"], fontsize=7.5)
    ax.set_ylabel("Intermembrane separation (nm)")
    ax.set_ylim(0, 33)
    ax.set_title("Single-conformer model span vs reference values", fontsize=10)
    fig.text(0.01, 0.005, f"*Assumes a ~{SURFACE_OFFSET:.0f} nm total offset between conventions; illustrative, not computed.",
             fontsize=6.5, ha="left", va="bottom")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    OUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUT)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
