#!/usr/bin/env python3
"""Download every external input ONCE into data/raw/ and record provenance in data/MANIFEST.tsv.

Usage:
  python src/fetch_inputs.py fetch [--force]   # download missing files, append manifest rows
  python src/fetch_inputs.py verify            # re-hash every file and compare with the manifest
  python src/fetch_inputs.py link              # copy pinned files to structures/ with the names the scripts expect
  python src/fetch_inputs.py compare --dir D   # check your existing structure files against the pinned downloads

After `fetch`, no other script may use the network; they read data/raw/ only.
"""
import argparse, csv, datetime, hashlib, json, pathlib, sys
import requests

ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
MANIFEST = ROOT / "data" / "MANIFEST.tsv"
FIELDS = ["file", "url", "utc_download", "sha256", "sha1", "bytes", "note"]

PDB_IDS = ["4ZQK", "3RRQ", "5B8C", "5WT9", "5X8L", "5X8M", "5GRJ", "5IUS", "5JXE"]
UNIPROT = {"PD1": "Q15116", "PDL1": "Q9NZQ7"}
# Names the phase*/revision scripts expect inside structures/
ALIASES = {
    "pdb/4ZQK.pdb": "4ZQK_pd1_pdl1_complex.pdb",
    "pdb/3RRQ.pdb": "3RRQ_pd1_unbound.pdb",
    "pdb/5B8C.pdb": "5B8C_pd1_pembrolizumab.pdb",
    "pdb/5WT9.pdb": "5WT9_pd1_nivolumab.pdb",
    "pdb/5X8L.pdb": "5X8L_pdl1_atezolizumab.pdb",
    "pdb/5X8M.pdb": "5X8M_pdl1_durvalumab.pdb",
    "pdb/5GRJ.pdb": "5GRJ_pdl1_avelumab.pdb",
    "alphafold/AF-Q15116-F1-model_v6.pdb": "PD1_alphafold_full.pdb",
    "alphafold/AF-Q9NZQ7-F1-model_v6.pdb": "PDL1_alphafold_full.pdb",
}
AF_VERSION = "v6"
# SHA-1 values stated in the manuscript for the AlphaFold files actually used
AF_EXPECTED_SHA1 = {
    "AF-Q15116-F1-model_v6.pdb": "8c23c10098430fc80e89287adab0aeecf7d8096b",
    "AF-Q9NZQ7-F1-model_v6.pdb": "068df4065f2f49baade0e4f7f7a7225728c69ebc",
}


def targets():
    t = []
    for p in PDB_IDS:
        t.append((f"pdb/{p}.cif", f"https://files.rcsb.org/download/{p}.cif"))
        t.append((f"rcsb_entry/{p}.json", f"https://data.rcsb.org/rest/v1/core/entry/{p}"))
        t.append((f"sifts/{p}.json", f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{p.lower()}"))
    for rel in ALIASES:  # PDB-format copies used by the phase scripts
        if rel.startswith("pdb/"):
            pid = rel[4:8]
            t.append((rel, f"https://files.rcsb.org/download/{pid}.pdb"))
    for acc in UNIPROT.values():
        t.append((f"uniprot/{acc}.fasta", f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"))
        t.append((f"uniprot/{acc}.json", f"https://rest.uniprot.org/uniprotkb/{acc}.json"))
        t.append((f"alphafold/AF-{acc}-F1-model_{AF_VERSION}.pdb",
                  f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_{AF_VERSION}.pdb"))
        t.append((f"alphafold/AF-{acc}-F1-predicted_aligned_error_{AF_VERSION}.json",
                  f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-predicted_aligned_error_{AF_VERSION}.json"))
    return t


def digest(path):
    b = pathlib.Path(path).read_bytes()
    return hashlib.sha256(b).hexdigest(), hashlib.sha1(b).hexdigest(), len(b)


def read_manifest():
    if not MANIFEST.exists():
        return {}
    with open(MANIFEST, newline="") as f:
        return {r["file"]: r for r in csv.DictReader(f, delimiter="\t")}


def write_manifest(rows):
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with open(MANIFEST, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        for k in sorted(rows):
            w.writerow(rows[k])


def fetch(force=False):
    rows = read_manifest()
    failed = []
    for rel, url in targets():
        dest = RAW / rel
        if dest.exists() and rel in rows and not force:
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        try:
            r = requests.get(url, timeout=120)
            r.raise_for_status()
        except Exception as e:  # keep going; report at the end
            failed.append((rel, str(e)))
            continue
        dest.write_bytes(r.content)
        sha256, sha1, n = digest(dest)
        note = ""
        if rel.startswith("uniprot/") and rel.endswith(".json"):
            note = "UniProt release header: " + r.headers.get("X-UniProt-Release", "NOT PROVIDED")
        if rel.startswith("rcsb_entry/"):
            rev = json.loads(r.content).get("rcsb_accession_info", {}).get("revision_date", "?")
            note = f"RCSB revision_date {rev}"
        name = dest.name
        if name in AF_EXPECTED_SHA1:
            note = "matches manuscript SHA-1" if sha1 == AF_EXPECTED_SHA1[name] else \
                   f"!! SHA-1 DIFFERS from manuscript ({AF_EXPECTED_SHA1[name]})"
        rows[rel] = dict(file=rel, url=url,
                         utc_download=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
                         sha256=sha256, sha1=sha1, bytes=n, note=note)
        print(f"fetched {rel} ({n} bytes) {note}")
    write_manifest(rows)
    if failed:
        print("\nFAILED downloads:")
        for rel, e in failed:
            print(" ", rel, "->", e)
        sys.exit(1)


def verify():
    rows = read_manifest()
    if not rows:
        sys.exit("MANIFEST.tsv missing or empty: run `fetch` first")
    bad = 0
    for rel, r in rows.items():
        p = RAW / rel
        if not p.exists():
            print("MISSING ", rel); bad += 1; continue
        sha256, _, _ = digest(p)
        if sha256 != r["sha256"]:
            print("CHANGED ", rel); bad += 1
    wanted = {rel for rel, _ in targets()}
    for rel in sorted(wanted - set(rows)):
        print("NOT IN MANIFEST", rel); bad += 1
    print("verify OK" if not bad else f"verify FAILED: {bad} problem(s)")
    sys.exit(1 if bad else 0)


def link():
    """Copy raw files to structures/ under the names the existing scripts expect."""
    import shutil
    out = ROOT / "structures"
    out.mkdir(exist_ok=True)
    for rel, alias in ALIASES.items():
        shutil.copyfile(RAW / rel, out / alias)
        print("structures/" + alias)


def _atoms(path):
    keep = set()
    for line in open(path, errors="replace"):
        if line.startswith(("ATOM", "HETATM")):
            keep.add(line[12:54])  # atom name, resname, chain, resnum, xyz
    return keep


def compare(local_dir):
    """Is each of your EXISTING structure files identical to the pinned download?"""
    local = pathlib.Path(local_dir).expanduser()
    for rel, alias in ALIASES.items():
        a, b = local / alias, RAW / rel
        if not a.exists():
            print(f"{alias:34} not in {local}"); continue
        if digest(a)[0] == digest(b)[0]:
            print(f"{alias:34} IDENTICAL (byte-for-byte)"); continue
        A, B = _atoms(a), _atoms(b)
        verdict = "same coordinates, header/other differences" if A == B else \
                  f"COORDINATES DIFFER (only local {len(A-B)}, only pinned {len(B-A)} atom records)"
        print(f"{alias:34} {verdict}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["fetch", "verify", "link", "compare"])
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--dir", help="folder of your existing structure files (compare)")
    a = ap.parse_args()
    if a.cmd == "fetch": fetch(a.force)
    elif a.cmd == "verify": verify()
    elif a.cmd == "link": link()
    else: compare(a.dir or sys.exit("compare needs --dir <folder>"))
