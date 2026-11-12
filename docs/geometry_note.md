# Feasible intermembrane span: model and assumptions

This note describes `src/feasible_span.py`, which replaces the unweighted tilt grid in `span_robustness.py` for the question
"how far apart can the two membranes be?". It was written from the script's docstring and outputs and checked against the code.

## Model

1. **Rigid ectodomain complex.** PD-1 (chain P, residues 25-145, excluding the AlphaFold loop 85-92) and PD-L1 (chain L, residues 19-225) are held in the 4ZQK frame. The input is the crystal-core hybrid (`outputs/hybrid_crystal_core.pdb`, built by `src/build_hybrid.py`). Only non-hydrogen atoms of standard residues are used, so waters, ligands and glycans are not included; the extent of the complex along N is taken from the convex hull of these atoms.
2. **Parallel membranes.** The two membranes are parallel planes with a common unit normal N pointing from the PD-1 membrane to the PD-L1 membrane. Heights are z = x.N. The midplanes are at heights Z1 < Z2, and each membrane surface lies a distance delta from its midplane (about half a bilayer).
3. **Non-penetration.** Every ectodomain atom lies between the two surfaces: Z1 <= zmin - delta and Z2 >= zmax + delta.
4. **Stalk reach.** The ectodomain C-terminus (PD-1 CA145, PD-L1 CA225) must be within the stalk contour length of its membrane surface (vertical reach; the TM anchor can move sideways). Contour length L = number of peptide steps x step length, with 26 steps for PD-1 (145 to 171) and 14 for PD-L1 (225 to 239).
5. **Span.** span = Z2 - Z1 (midplane to midplane). For each direction N:
   S_min = (zmax - zmin) + 2 delta, and S_max = z(CA225) - z(CA145) + L1 + L2 + 2 delta.
   N is feasible only if both stalks can reach: z(CA145) - zmin <= L1 and zmax - z(CA225) <= L2.
6. **Orientation prior.** Uniform over directions N (Fibonacci sphere). AlphaFold gives no information on TM orientation (PAE about 25-28 A), so none is preferred. The tilt of the TM helices relative to N is not constrained, because the stalks (25 and 13 residues) can reorient them. The number of directions is 10,000.
7. **Optional hinge.** PD-L1 residues 133-225 (IgC) are rotated about the CA of residue 132 by a rotation about a random axis through an angle drawn uniformly between 0 and the stated hinge angle. Each hinge setting above 0 uses 60 such random rotations; hinge 0 is a single rigid evaluation.

## Parameters scanned

Half-bilayer thickness delta 1.5, 2.0 and 2.5 nm; step length 3.3, 3.5 and 3.8 A per residue; hinge 0, 20 and 40 degrees.
Directions: 10,000 (Fibonacci sphere). Hinge rotations: 60 per hinge setting above 0. Random seed: 1 (`numpy.random.default_rng(1)`), so a re-run gives identical output files.

## Reading the output

| Column | Meaning |
|---|---|
| feasible | fraction of directions (pooled over the hinge rotations of that setting) for which both stalks can reach the membranes |
| Smin low | lowest S_min over all feasible directions and hinge rotations |
| Smin med | median S_min over the same feasible set |
| Smax med | median S_max over the same feasible set |
| Smax high | highest S_max over the same feasible set |

"Low" and "high" are extremes, not percentiles or confidence limits.

## Central result

For delta = 2.0 nm, 3.5 A per residue and a rigid complex: 71.7% of directions feasible; S_min 8.1 nm (lowest) and 9.8 nm (median); S_max 20.0 nm (median) and 25.0 nm (highest).

## How to interpret it

- Differences between the hinge rows are small (for example 0.692, 0.693 and 0.696 feasible at 3.3 A) and may partly reflect sampling noise from 60 random rotations per setting with one seed, so no conclusion about the hinge effect should be drawn from them.
- It is a geometric feasibility bound, not a probability distribution of the real separation.
- The maximum assumes fully extended stalks, so it is an upper bound.
- Glycans and lipid detail are not modelled, so real minimum separations would be larger.
- The 15.2-18.1 nm spans from the AlphaFold-oriented hybrid lie inside the feasible range, consistent with the conclusion that the span is not tightly constrained by this model.
