#!/usr/bin/env python3
import csv
from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa

CUTOFF = 4.5

STRUCTURES = [
    {
        "file": "structures/5B8C_pd1_pembrolizumab.pdb",
        "label": "PD-1 + pembrolizumab",
        "antigen_chain": "C",
        "antibody_chains": ["A", "B", "D", "E", "G", "H", "J", "K"],
    },
    {
        "file": "structures/5WT9_pd1_nivolumab.pdb",
        "label": "PD-1 + nivolumab",
        "antigen_chain": "G",
        "antibody_chains": ["H", "L"],
    },
    {
        "file": "structures/5X8L_pdl1_atezolizumab.pdb",
        "label": "PD-L1 + atezolizumab",
        "antigen_chain": "A",
        "antibody_chains": ["F", "G", "H", "S", "J", "K", "L", "M", "N", "O"],
    },
    {
        "file": "structures/5X8M_pdl1_durvalumab.pdb",
        "label": "PD-L1 + durvalumab",
        "antigen_chain": "A",
        "antibody_chains": ["B", "C"],
    },
    {
        "file": "structures/5GRJ_pdl1_avelumab.pdb",
        "label": "PD-L1 + avelumab",
        "antigen_chain": "A",
        "antibody_chains": ["H", "L"],
    },
]


def heavy_atoms(chain):
    atoms = []
    for res in chain:
        if not is_aa(res, standard=True):
            continue
        for atom in res:
            if atom.element != "H":
                atoms.append(atom)
    return atoms


def find_epitope(structure_path, antigen_chain_id, antibody_chain_ids):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("x", structure_path)[0]

    antigen_chain = structure[antigen_chain_id]
    antibody_atoms = []
    for cid in antibody_chain_ids:
        antibody_atoms.extend(heavy_atoms(structure[cid]))

    ns = NeighborSearch(antibody_atoms)

    epitope_residues = set()
    for res in antigen_chain:
        if not is_aa(res, standard=True):
            continue
        for atom in res:
            if atom.element == "H":
                continue
            nearby = ns.search(atom.coord, CUTOFF)
            if nearby:
                resnum = res.get_id()[1]
                resname = res.get_resname()
                epitope_residues.add((resnum, resname))
                break

    return sorted(epitope_residues)


def main():
    all_results = []

    for cfg in STRUCTURES:
        print(f"\n=== {cfg['label']} ({cfg['file']}) ===")
        try:
            epitope = find_epitope(cfg["file"], cfg["antigen_chain"], cfg["antibody_chains"])
        except KeyError as e:
            print(f"  !! Chain {e} not found -- check:")
            print(f"     grep COMPND {cfg['file']}")
            continue
        except FileNotFoundError:
            print(f"  !! File not found: {cfg['file']} -- run Phase 1 first.")
            continue

        if not epitope:
            print("  !! Zero contact residues found -- almost certainly wrong chain IDs.")
            continue

        print(f"  {len(epitope)} epitope residues (within {CUTOFF} A heavy-atom cutoff):")
        residue_str = ", ".join(f"{resname}{resnum}" for resnum, resname in epitope)
        print(f"  {residue_str}")

        for resnum, resname in epitope:
            all_results.append({
                "structure": cfg["label"],
                "antibody": cfg["label"].split("+")[-1].strip(),
                "residue_number": resnum,
                "residue_name": resname,
            })

    if all_results:
        with open("phase5_epitope_residues.csv", "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["structure", "antibody", "residue_number", "residue_name"])
            writer.writeheader()
            writer.writerows(all_results)
        print(f"\nWrote phase5_epitope_residues.csv ({len(all_results)} rows total)")


if __name__ == "__main__":
    main()
