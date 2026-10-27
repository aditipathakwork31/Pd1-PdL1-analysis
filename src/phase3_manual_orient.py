#!/usr/bin/env python3
import numpy as np
from Bio.PDB import PDBParser, PDBIO

CONFIG = [
    {
        "input": "structures/PD1_alphafold_full.pdb",
        "output": "structures/PD1_ppm_oriented.pdb",
        "tm_start": 171,
        "tm_end": 191,
    },
    {
        "input": "structures/PDL1_alphafold_full.pdb",
        "output": "structures/PDL1_ppm_oriented.pdb",
        "tm_start": 239,
        "tm_end": 259,
    },
]


def get_tm_ca_coords(chain, tm_start, tm_end):
    coords = []
    for res in chain:
        resnum = res.get_id()[1]
        if tm_start <= resnum <= tm_end and "CA" in res:
            coords.append(res["CA"].coord)
    return np.array(coords)


def orient_structure(input_path, output_path, tm_start, tm_end):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("x", input_path)
    model = structure[0]
    chain = next(iter(model))

    tm_coords = get_tm_ca_coords(chain, tm_start, tm_end)
    if len(tm_coords) < 5:
        raise ValueError(
            f"Only found {len(tm_coords)} CA atoms in TM range {tm_start}-{tm_end} "
            f"for {input_path} -- check residue numbering."
        )

    centroid = tm_coords.mean(axis=0)
    centered = tm_coords - centroid
    _, _, vh = np.linalg.svd(centered)
    helix_axis = vh[0]

    z_axis = np.array([0.0, 0.0, 1.0])
    v = np.cross(helix_axis, z_axis)
    s = np.linalg.norm(v)
    c = np.dot(helix_axis, z_axis)

    if s < 1e-8:
        rot = np.eye(3) if c > 0 else -np.eye(3)
    else:
        vx = np.array([
            [0, -v[2], v[1]],
            [v[2], 0, -v[0]],
            [-v[1], v[0], 0],
        ])
        rot = np.eye(3) + vx + vx @ vx * ((1 - c) / (s ** 2))

    for atom in structure.get_atoms():
        atom.coord = np.dot(rot, atom.coord - centroid)

    rotated_tm = np.array([np.dot(rot, c - centroid) for c in tm_coords])
    z_mid = rotated_tm[:, 2].mean()
    for atom in structure.get_atoms():
        atom.coord[2] -= z_mid

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_path)
    print(f"  Wrote {output_path}  (TM axis aligned to Z, midpoint at Z=0)")


def main():
    print("Manual membrane orientation (PPM fallback)")
    print("=" * 60)
    for cfg in CONFIG:
        print(f"\nProcessing {cfg['input']} (TM helix {cfg['tm_start']}-{cfg['tm_end']})...")
        try:
            orient_structure(cfg["input"], cfg["output"], cfg["tm_start"], cfg["tm_end"])
        except FileNotFoundError:
            print(f"  !! File not found: {cfg['input']} -- run Phase 1 first.")
        except ValueError as e:
            print(f"  !! {e}")

    print("\nDone. These outputs use the SAME filenames Phase 4 expects.")


if __name__ == "__main__":
    main()
