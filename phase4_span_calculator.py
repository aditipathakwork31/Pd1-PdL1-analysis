#!/usr/bin/env python3
import sys
import numpy as np
from Bio.PDB import PDBParser, Superimposer, PDBIO
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner

CRYSTAL_FILE = "structures/4ZQK_pd1_pdl1_complex.pdb"
CRYSTAL_PD1_CHAIN = "B"
CRYSTAL_PDL1_CHAIN = "A"

PD1_PPM_FILE = "structures/PD1_ppm_oriented.pdb"
PDL1_PPM_FILE = "structures/PDL1_ppm_oriented.pdb"

OUTPUT_GRAFTED_PDB = "phase4_grafted_model.pdb"
MIN_MATCH_LENGTH = 30


def get_chain_seq_and_residues(chain):
    seq = []
    residues = []
    for res in chain:
        if not is_aa(res, standard=True):
            continue
        if "CA" not in res:
            continue
        try:
            aa = seq1(res.get_resname())
        except KeyError:
            continue
        seq.append(aa)
        residues.append(res)
    return "".join(seq), residues


def align_and_get_transform(fixed_chain, moving_chain, label):
    fixed_seq, fixed_res = get_chain_seq_and_residues(fixed_chain)
    moving_seq, moving_res = get_chain_seq_and_residues(moving_chain)

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.match_score = 2
    aligner.mismatch_score = -1

    alignment = aligner.align(fixed_seq, moving_seq)[0]
    fixed_aligned, moving_aligned = alignment.aligned

    fixed_atoms = []
    moving_atoms = []
    n_identical = 0

    for (f_start, f_end), (m_start, m_end) in zip(fixed_aligned, moving_aligned):
        for offset in range(f_end - f_start):
            fi = f_start + offset
            mi = m_start + offset
            if fixed_seq[fi] == moving_seq[mi]:
                fixed_atoms.append(fixed_res[fi]["CA"])
                moving_atoms.append(moving_res[mi]["CA"])
                n_identical += 1

    if n_identical < MIN_MATCH_LENGTH:
        sys.exit(
            f"[{label}] Only {n_identical} identical-residue CA pairs found -- aborting."
        )

    sup = Superimposer()
    sup.set_atoms(fixed_atoms, moving_atoms)
    rot, tran = sup.rotran
    print(f"[{label}] Aligned on {n_identical} identical CA pairs, RMSD = {sup.rms:.2f} Å")

    return rot, tran, sup.rms, n_identical


def transform_structure(structure, rot, tran):
    for atom in structure.get_atoms():
        atom.coord = np.dot(atom.coord, rot) + tran


def membrane_plane_in_new_frame(rot):
    normal = np.dot(np.array([0.0, 0.0, 1.0]), rot)
    normal = normal / np.linalg.norm(normal)
    return normal


def main():
    parser = PDBParser(QUIET=True)

    print("Loading structures...")
    crystal = parser.get_structure("crystal", CRYSTAL_FILE)[0]
    pd1_ppm = parser.get_structure("pd1_ppm", PD1_PPM_FILE)[0]
    pdl1_ppm = parser.get_structure("pdl1_ppm", PDL1_PPM_FILE)[0]

    try:
        crystal_pd1_chain = crystal[CRYSTAL_PD1_CHAIN]
        crystal_pdl1_chain = crystal[CRYSTAL_PDL1_CHAIN]
    except KeyError as e:
        sys.exit(f"Chain {e} not found in {CRYSTAL_FILE}.")

    pd1_ppm_chain = next(iter(pd1_ppm))
    pdl1_ppm_chain = next(iter(pdl1_ppm))

    print("\n--- Aligning PD-1 ---")
    rot1, tran1, rmsd1, n1 = align_and_get_transform(crystal_pd1_chain, pd1_ppm_chain, "PD-1")

    print("\n--- Aligning PD-L1 ---")
    rot2, tran2, rmsd2, n2 = align_and_get_transform(crystal_pdl1_chain, pdl1_ppm_chain, "PD-L1")

    transform_structure(pd1_ppm, rot1, tran1)
    transform_structure(pdl1_ppm, rot2, tran2)

    pd1_normal = membrane_plane_in_new_frame(rot1)
    pd1_point = tran1
    pdl1_normal = membrane_plane_in_new_frame(rot2)
    pdl1_point = tran2

    cos_angle = np.clip(np.dot(pd1_normal, -pdl1_normal), -1.0, 1.0)
    angle_deg = np.degrees(np.arccos(cos_angle))

    diff = pdl1_point - pd1_point
    span_along_pd1_normal = abs(np.dot(diff, pd1_normal))
    span_along_pdl1_normal = abs(np.dot(diff, pdl1_normal))
    span_direct = np.linalg.norm(diff)

    print("\n" + "=" * 60)
    print("RESULTS")
    print("=" * 60)
    print(f"PD-1 alignment:  {n1} matched residues, RMSD {rmsd1:.2f} Å")
    print(f"PD-L1 alignment: {n2} matched residues, RMSD {rmsd2:.2f} Å")
    print(f"Membrane normal parallelism (should be <~20°): {angle_deg:.1f}°")
    if angle_deg > 20:
        print("  !! WARNING: check geometry in ChimeraX before trusting the span.")
    print()
    print(f"Span (along PD-1 normal):  {span_along_pd1_normal:.1f} Å ({span_along_pd1_normal/10:.2f} nm)")
    print(f"Span (along PD-L1 normal): {span_along_pdl1_normal:.1f} Å ({span_along_pdl1_normal/10:.2f} nm)")
    print(f"Direct point distance:     {span_direct:.1f} Å ({span_direct/10:.2f} nm)")

    io = PDBIO()
    from Bio.PDB import Structure, Model
    combined = Structure.Structure("grafted")
    model = Model.Model(0)
    combined.add(model)
    for chain in crystal:
        model.add(chain.copy())
    for chain in pd1_ppm:
        chain.id = "X"
        model.add(chain)
    for chain in pdl1_ppm:
        chain.id = "Y"
        model.add(chain)
    io.set_structure(combined)
    io.save(OUTPUT_GRAFTED_PDB)
    print(f"\nWrote {OUTPUT_GRAFTED_PDB} -- open in ChimeraX to cross-check the span manually.")


if __name__ == "__main__":
    main()
