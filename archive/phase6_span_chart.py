#!/usr/bin/env python3
import matplotlib.pyplot as plt

PD1_PDL1_SPAN_NM = 16.5  # range 15.2-18.1 nm across the three estimates, see Phase 4 output

RECEPTOR_PAIRS = {
    "PD-1:PD-L1":     {"span_nm": PD1_PDL1_SPAN_NM, "source": "This work (Phase 4)"},
    "TCR:pMHC":       {"span_nm": 13.0, "source": "Cordoba et al. 2013, Blood 121(21):4295"},
    "IgSF pairs (general)": {"span_nm": 15.0, "source": "Strazza et al. 2024, Immunology (review)"},
    "CD28:CD86":      {"span_nm": None, "source": "Not found with a distinct sourced value"},
    "CTLA-4:CD80":    {"span_nm": None, "source": "Not found with a distinct sourced value"},
    "LAG-3:MHC-II":   {"span_nm": None, "source": "Not found with a distinct sourced value"},
}

CD45_R0_HEIGHT_NM = 25.0
CD45_RABC_HEIGHT_NM = 40.0
CD45_HEIGHT_NM = CD45_R0_HEIGHT_NM


def main():
    labels = []
    values = []
    colors = []

    for name, info in RECEPTOR_PAIRS.items():
        if info["span_nm"] is None:
            print(f"  Skipping '{name}' -- no value filled in yet ({info['source']})")
            continue
        labels.append(name)
        values.append(info["span_nm"])
        colors.append("#4C72B0" if name == "PD-1:PD-L1" else "#8C8C8C")

    if not labels:
        print("\nNo values filled in yet.")
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(labels, values, color=colors)

    if CD45_HEIGHT_NM is not None:
        ax.axhline(CD45_HEIGHT_NM, color="firebrick", linestyle="--", linewidth=1.5,
                    label=f"CD45 ectodomain height (~{CD45_HEIGHT_NM} nm) -- sterically excluded, not a bound span")
        ax.legend(fontsize=8)

    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, val + 0.2, f"{val:.1f}",
                 ha="center", va="bottom", fontsize=9)

    ax.set_ylabel("Intermembrane span (nm)")
    ax.set_title("Immune synapse receptor:ligand intermembrane spans")
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    plt.savefig("phase6_span_comparison.png", dpi=300)
    print("Saved phase6_span_comparison.png")


if __name__ == "__main__":
    main()
