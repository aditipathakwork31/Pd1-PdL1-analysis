[requirements.txt](https://github.com/user-attachments/files/33874924/requirements.txt)
[run_log.txt](https://github.com/user-attachments/files/33874838/run_log.txt)
# PD-1:PD-L1 geometric model and therapeutic antibody epitope mapping

Analysis scripts for "A Membrane-Anchored Geometric Model of PD-1:PD-L1 and Therapeutic Antibody Epitope Mapping".
All results come from public structures (RCSB, AlphaFold DB, UniProt). Nothing here is experimental.

## Read this first: superseded scripts

| Script | Status |
|---|---|
| `phase5_epitope_mapper.py` | **Superseded.** Pairs ONE antigen copy with the antibody chains of ALL copies, so crystal-packing contacts are counted. It gives 45 residues for pembrolizumab (5B8C chain C) and 40 for atezolizumab (5X8L chain A). The corrected counts (20 to 29 and 24 to 30) come from `pd1_pdl1_revision_checks.py epitope`. |
| `phase6_span_chart.py` | **Superseded.** Draws a single 16.5 nm bar. Use `make_figure1.py`, which shows the range and the stated assumptions. |

## Environment

Python 3.14.7, Biopython 1.88, NumPy 2.5.3, requests 2.34.2, Matplotlib 3.11.2 (`requirements.txt`).
UCSF ChimeraX 1.13.dev202609011810 (2026-09-01), used for Figure 2 and the independent checks.
AlphaFold database version 6 (AF-Q15116-F1 and AF-Q9NZQ7-F1, model date 2025-08-01).
SHA-1 of the AlphaFold files used: PD-1 `8c23c10098430fc80e89287adab0aeecf7d8096b`, PD-L1 `068df4065f2f49baade0e4f7f7a7225728c69ebc`.
UniProt release 2026_02 (effective 10 June 2026); records pulled 1 September 2026, one day before release 2026_03 went live, so 2026_02 is the release these files reflect. Entry versions at pull: PD-1 (Q15116) 207, PD-L1 (Q9NZQ7) 206.
RCSB/PDBe queries were run in October 2026. Database content can change, so re-running later may differ slightly.

## Inputs

The `phase*` scripts read local files with fixed names from `structures/`. They do not download anything.
`pd1_pdl1_revision_checks.py` downloads from RCSB, PDBe and UniProt each time it runs.

creates `structures/` with these names: `4ZQK_pd1_pdl1_complex.pdb`, `5B8C_pd1_pembrolizumab.pdb`, `5WT9_pd1_nivolumab.pdb`,
`5X8L_pdl1_atezolizumab.pdb`, `5X8M_pdl1_durvalumab.pdb`, `5GRJ_pdl1_avelumab.pdb`, `PD1_alphafold_full.pdb`, `PDL1_alphafold_full.pdb`.
It also checks the two AlphaFold files against the SHA-1 values above and prints the chains of each crystal structure.

Check that the inputs are right: 4ZQK chain A must be PD-L1 and chain B PD-1; `specs.json` lists the antigen, heavy and light chains for the antibody structures.

## Run order
Place the input structures in `structures/` (see Inputs), then:

1. `python pd1_pdl1_revision_checks.py selftest`, then `meta`, `search`, `mapping --residue 78`, `epitope --specs specs.json`, `apo`, and `interface --pdb 5IUS --pd1-chain A --pdl1-chain C` (and B/D). These download from RCSB/PDBe/UniProt and reproduce the epitope, apo and interface results.
2. `python phase3_manual_orient.py` orients the AlphaFold models from the TM helix axis.
3. `python phase4_span_calculator.py` builds the hybrid model (`phase4_grafted_model.pdb`) and prints three spans.
4. `python span_robustness.py phase4_grafted_model.pdb` gives span under different normals, tilt sensitivity, clash count (195, falling to 21 after trimming) and pLDDT.
5. `python trim_hybrid.py phase4_grafted_model.pdb phase4_grafted_model_trimmed.pdb` removes signal peptides and cytoplasmic tails.

The full console output of steps 2-5 is in `run_log.txt`.
## Expected results (use these to check a re-run)

- `selftest` prints `selftest OK`. `mapping` prints UniProt Q15116 residue 78 = K and chain B author residue 78 = LYS in 4ZQK.
- `search` returns exactly two entries, 4ZQK and 5IUS. The first line of `epitope` output gives 18 residues for the 4ZQK PD-L1 interface at 4.5 A.
- `epitope` at 4.5 A: nivolumab 14 residues; pembrolizumab 20 (consensus) to 29 (union), per copy 27, 27, 23, 22;
  atezolizumab 24 to 30, per copy 27, 25, 26, 26, 25; durvalumab 17; avelumab 20.
- `span_robustness.py` on the full hybrid: angle between the ectodomain-pointing normals 133.0 degrees (47.0 degrees from antiparallel);
  spans 15.22 nm (along the PD-1 normal), 17.39 nm (PD-L1 normal), 18.10 nm (direct), 17.78 nm (common normal);
  195 heavy-atom pairs closer than 2.0 A between the two AlphaFold chains. On the trimmed hybrid the spans are unchanged and the pairs fall to 21.
- `interface` for 5IUS: pairs A-C and B-D are the real contacts. PD-L1 interface Jaccard 0.80 with 4ZQK in both; PD-1 interface Jaccard 0.42 (A-C) and 0.46 (B-D).
- `apo` (3RRQ against 4ZQK PD-1): 105 matched residues, all-atom fit 1.52 A, core fit 0.63 A, interface 2.93 A under the core fit;
  residues 74 and 75 deviate by 8.1 A and residue 76 by 5.1 A.

## Scripts

| File | Purpose |
|---|---|
| `phase2_domain_mapping.py` | UniProt domain architecture (writes `domain_architecture.txt` and the FASTA files) |
| `phase3_manual_orient.py`, `phase4_span_calculator.py` | Membrane orientation from the TM helix axis; grafting of full-length AlphaFold chains onto 4ZQK; span estimates |
| `span_robustness.py` | Normals with a defined direction, tilt sensitivity, clash check, pLDDT |
| `trim_hybrid.py` | Removes signal peptides and cytoplasmic tails from the hybrid |
| `pd1_pdl1_revision_checks.py` (with `specs.json`) | Metadata, RCSB search, numbering, epitopes with antibody pairing, 5IUS interface, bound versus unbound PD-1, PD-1 loop ensemble |
| `make_figure1.py`, `render_figure2.cxc` | Figures 1 and 2 |
| `item12_epitope_crosscheck.cxc`, `compare_chimerax_results.py` | Independent epitope cross-check in ChimeraX |
| `item13_3rrq_crystal_contacts.cxc`, `bfactor_normalised_71_76.py` | Loop 71 to 76 checks |
| `item10_provenance.sh` | Records versions and compares UniProt sequences |

## Limitations

- The membrane normal is taken from the TM helix axis. Tilting it by up to 30 degrees changes the span along the PD-1 normal to 8.4 to 18.0 nm, so the 15 to 18 nm bracket is not a constrained result.
- The hybrid model contains single static AlphaFold stalks (low confidence) and no lipids or glycans.
- Epitopes are defined by a heavy-atom distance cutoff, not by energy.
- The computed PD-1:PD-L1 span (15.2-18.1 nm; 8.4-18.0 nm under tilt sensitivity) is of the same order as the ~10-15 nm extracellular span widely cited for immune receptor-ligand pairs, originating in crystal structures of KIR2DL1-HLA-Cw4 and NKG2D-MICA (Fan et al. 2001; Li et al. 2001; Boyington et al. 2000) and propagated through subsequent reviews (Alakoskela et al. 2011; Strazza et al. 2024). The exact measurement convention underlying the original ~15 nm figure (e.g., whether end-to-end or membrane-normal projected length) is not stated in the reviewing literature and was not independently re-derived here; the comparison is therefore qualitative, not a matched metric. The CD45 ectodomain height (~28-50 nm across splice isoforms; Woollett et al. 1985, McCall et al. 1992, as cited in Cordoba et al. 2013) is well above this range, consistent with CD45's size-based exclusion from close-contact zones under the kinetic-segregation model.
