#!/usr/bin/env python3
import json
import textwrap
import urllib.request

PROTEINS = {
    "PD1":  {"uniprot": "Q15116", "name": "Programmed cell death protein 1 (PDCD1)"},
    "PDL1": {"uniprot": "Q9NZQ7", "name": "Programmed cell death 1 ligand 1 (CD274)"},
}

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{}.json"


def fetch_uniprot(acc):
    url = UNIPROT_URL.format(acc)
    with urllib.request.urlopen(url, timeout=30) as resp:
        return json.loads(resp.read().decode())


def get_sequence(entry):
    return entry["sequence"]["value"]


def summarize_features(entry):
    wanted_types = {
        "Signal": "SIGNAL PEPTIDE",
        "Transmembrane": "TM HELIX",
        "Topological domain": "TOPOLOGICAL DOMAIN",
        "Domain": "DOMAIN (Ig-fold etc.)",
        "Glycosylation": "N-GLYCOSYLATION SITE",
    }

    results = {label: [] for label in wanted_types.values()}

    for feat in entry.get("features", []):
        ftype = feat.get("type")
        if ftype not in wanted_types:
            continue
        label = wanted_types[ftype]
        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        desc = feat.get("description", "")
        if start == end:
            span = f"{start}"
        else:
            span = f"{start}-{end}"
        results[label].append(f"{span}: {desc}".strip(": "))

    return results


def write_fasta(name, header, seq, path):
    with open(path, "w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")


def main():
    report_lines = []
    report_lines.append("=" * 70)
    report_lines.append("PD-1 / PD-L1 DOMAIN ARCHITECTURE  (live UniProt pull)")
    report_lines.append("=" * 70)
    report_lines.append("")

    for key, info in PROTEINS.items():
        acc = info["uniprot"]
        print(f"Fetching UniProt entry {acc} ({key})...")
        try:
            entry = fetch_uniprot(acc)
        except Exception as e:
            print(f"  !! FAILED to fetch {acc}: {e}")
            continue

        seq = get_sequence(entry)
        seq_len = entry["sequence"]["length"]
        feats = summarize_features(entry)

        write_fasta(key, f"{acc}|{key}|{info['name']}", seq, f"{key}.fasta")

        report_lines.append(f"--- {key}  (UniProt {acc})  {info['name']} ---")
        report_lines.append(f"Full length: {seq_len} aa")
        report_lines.append("")

        for label, items in feats.items():
            report_lines.append(f"  {label}:")
            if not items:
                report_lines.append("    (none annotated)")
            for item in items:
                report_lines.append(f"    - {item}")
            report_lines.append("")

        report_lines.append("-" * 70)
        report_lines.append("")

    report_lines.append(
        "NOTES:\n"
        "- TM HELIX here is UniProt's curated annotation. Cross-validate it\n"
        "  independently via TMHMM 2.0 or DeepTMHMM if you want a second\n"
        "  opinion before trusting it for membrane orientation.\n"
        "- N-glycosylation sites won't be modeled in the final geometry but\n"
        "  matter for interpretation near the stalk/ectodomain junction.\n"
        "- Use TOPOLOGICAL DOMAIN entries (not just TM boundaries) to define\n"
        "  ectodomain vs stalk vs TM vs cytoplasmic cut points.\n"
    )

    out_path = "domain_architecture.txt"
    with open(out_path, "w") as f:
        f.write("\n".join(report_lines))

    print()
    print(f"Done. Wrote {out_path}, PD1.fasta, PDL1.fasta")
    print()
    print("\n".join(report_lines))


if __name__ == "__main__":
    main()
