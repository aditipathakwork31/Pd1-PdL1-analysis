#!/usr/bin/env python3
"""Compare the ChimeraX output files from item12_epitope_crosscheck.cxc with the manuscript values.
Usage: python3 compare_chimerax_results.py [folder_with_txt_files]   (default: ~/Desktop)"""
import re, sys
from pathlib import Path

folder = Path(sys.argv[1]).expanduser() if len(sys.argv) > 1 else Path("~/Desktop").expanduser()

D45 = [23,26,56,58,59,60,61,111,112,113,115,121,122,123,124,125,127]   # 5X8M, 4.5 A (17)
A45 = [54,56,58,59,60,61,62,63,66,69,73,75,76,78,113,115,117,121,123,125]  # 5GRJ, 4.5 A (20)
# (file, expected_low, expected_high, expected residue list or None, label)
CHECKS = [
 ("5x8m_4.0.txt", 15, 17, None, "5X8M durvalumab 4.0 A (manuscript range 15-17; pasted plan said 15)"),
 ("5x8m_4.5.txt", 17, 17, D45, "5X8M durvalumab 4.5 A"),
 ("5x8m_5.0.txt", 17, 17, None, "5X8M durvalumab 5.0 A"),
 ("5grj_4.0.txt", 17, 17, None, "5GRJ avelumab 4.0 A"),
 ("5grj_4.5.txt", 20, 20, A45, "5GRJ avelumab 4.5 A"),
 ("5grj_5.0.txt", 22, 22, None, "5GRJ avelumab 5.0 A"),
 ("5b8c_C_own.txt", 27, 27, None, "5B8C chain C, own Fv (A,B)"),
 ("5b8c_C_all.txt", 45, 45, None, "5B8C chain C, all 8 antibody chains (original method)"),
 ("5b8c_I_own.txt", 23, 23, None, "5B8C chain I, own Fv (G,H)"),
 ("5x8l_A_all.txt", 40, 40, None, "5X8L chain A, all 10 antibody chains (original method)"),
 ("5x8l_A_own.txt", 24, 30, None, "5X8L chain A, own pair (F,K); manuscript consensus 24 / union 30"),
]

def read_residues(path):
    res = set()
    for line in path.read_text().splitlines():
        m = re.search(r"residue id /(\S+?):(-?\d+)", line)
        if m:
            res.add((m.group(1), int(m.group(2))))
    return res

print(f"Reading from {folder}\n")
for fname, lo, hi, exp_list, label in CHECKS:
    p = folder / fname
    if not p.exists():
        print(f"MISSING  {label}: {fname} not found"); continue
    res = read_residues(p)
    n = len(res)
    if lo <= n <= hi: status = "MATCH  "
    elif lo - 2 <= n <= hi + 2: status = "CLOSE  "   # boundary effects: report, do not change anything
    else: status = "DIFFERS"
    exp = f"{lo}" if lo == hi else f"{lo}-{hi}"
    print(f"{status}  {label}: found {n}, expected {exp}")
    if exp_list is not None:
        got = sorted(r[1] for r in res)
        miss = sorted(set(exp_list) - set(got)); extra = sorted(set(got) - set(exp_list))
        print(f"          residue list: missing {miss or 'none'}; extra {extra or 'none'}")
print("\nCLOSE = off by 1-2 (possible 4.5 A boundary difference): send me the lists before changing the manuscript.")
print("DIFFERS = off by more: do not submit until the cause is found.")
