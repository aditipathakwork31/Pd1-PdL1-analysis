# PD-1:PD-L1 geometric model and therapeutic antibody epitope mapping

Analysis scripts for "A Membrane-Anchored Geometric Model of PD-1:PD-L1 and Therapeutic Antibody Epitope Mapping".
All results come from public structures (RCSB, AlphaFold DB, UniProt). Nothing here is experimental.
Status: manuscript in preparation (not yet submitted).

**Research question.** If full-length AlphaFold chains of PD-1 and PD-L1 are placed on the 4ZQK crystal complex and anchored in a membrane, how far apart can the two membranes be, how sensitive is that estimate to membrane orientation, and which residues do five approved therapeutic antibodies contact once crystal-packing contacts are excluded?

## Quick start

```bash
make env        # create .venv, install requirements.txt, record versions in outputs/versions.txt
make fetch      # download every external input once into data/raw/ and record hashes in data/MANIFEST.tsv
make verify     # re-hash the downloads and compare with the manifest
make link       # copy the pinned files into structures/ under the names the scripts expect
make test       # unit and cross-script tests
make geometry   # TM-axis orientation, grafted hybrid, span robustness (writes outputs/)
make hybrid     # crystal-core hybrid, clash table, AF-vs-crystal deviations, confidence, membrane check, glycan screen, feasible span
```

After `make fetch`, no script uses the network except the RCSB/PDBe/UniProt steps of `src/pd1_pdl1_revision_checks.py`.
Not yet part of `make`: the epitope, interface and bound-versus-unbound steps, `src/item10_provenance.sh`, the ChimeraX checks, and the figures.
Run those directly, for example `python src/pd1_pdl1_revision_checks.py selftest`, then `meta`, `search`, `mapping --residue 78`, `epitope --specs data/specs.json`, `apo`, and `interface --pdb 5IUS --pd1-chain A --pdl1-chain C` (and B/D). TODO: confirm each command and the `specs.json` path on a clean checkout, and add `diagnose` and `pd1ensemble` if they are part of the workflow.

## Repository layout

| Path | Contents |
|---|---|
| `src/` | Analysis scripts (see Scripts) |
| `tests/` | `test_geometry.py`, `test_cross_script.py`, `test_hybrid_clash.py`, `test_step3b.py` |
| `data/` | `raw/` (pinned downloads), `MANIFEST.tsv` (URL, UTC time, SHA-256 and SHA-1 of each file), `specs.json` |
| `expected/` | Saved baseline output (`span_robustness.txt`) to compare a re-run against |
| `outputs/` | Outputs of the run sequence |
| `docs/` | `run_log_original.txt` (console output of the original run) and `geometry_note.md` (assumptions of the feasible-span model) |
| `archive/` | Superseded scripts, kept for the record |
| `PD1.fasta`, `PDL1.fasta`, `domain_architecture.txt` | UniProt sequences and domain annotation (pulled 1 September 2026) |

## Inputs and provenance

- **Original run (1 September 2026):** the input structures were downloaded with `phase1_download.sh` (seven RCSB entries and the two AlphaFold models, written to `structures/` with the file names the scripts expect). TODO: add this script to `archive/` unchanged and say here that it is the original.
- **Pinned re-download (26-27 October 2026):** `src/fetch_inputs.py` downloads every external input once into `data/raw/`, records URL, UTC time and hashes in `data/MANIFEST.tsv`, checks the AlphaFold files against the SHA-1 values below, and `link` copies them into `structures/`. It was written after the original run, so it reproduces the inputs and does not claim to be the script that produced them.
- **Do the original files match the pinned ones?** Run `make fetch` and then `make compare DIR=$HOME/pd1_geometry/structures`. TODO: record the result here (identical, or which files differ and how). Entries on the RCSB and AlphaFold can be revised after the original download.

Check that the inputs are right: 4ZQK chain A must be PD-L1 and chain B PD-1; `data/specs.json` lists the antigen, heavy and light chains for the antibody structures.

## Read this first: superseded scripts

| Script | Status |
|---|---|
| `archive/phase5_epitope_mapper.py` | **Superseded.** Pairs ONE antigen copy with the antibody chains of ALL copies, so crystal-packing contacts are counted. It gives 45 residues for pembrolizumab (5B8C chain C) and 40 for atezolizumab (5X8L chain A). The corrected counts (20 to 29 and 24 to 30) come from `src/pd1_pdl1_revision_checks.py epitope`. |
| `archive/phase6_span_chart.py` | **Superseded.** Draws a single 16.5 nm bar. |

TODO: Figure 1 (the span range with its stated assumptions) was meant to come from `make_figure1.py`, which is not in the repository. Either add the script that produced it or state how the figure was made.

## Environment

Python 3.14.7, Biopython 1.88, NumPy 2.5.3, requests 2.34.2, Matplotlib 3.11.2 (`requirements.txt`). The revision script ran in a virtual environment (which adds requests and Matplotlib); the original phase scripts are assumed to have used the system python3 with the same Biopython and NumPy (not separately recorded). UCSF ChimeraX 1.13.dev202609011810 (2026-09-01), used for Figure 2 and the independent checks. AlphaFold database version 6 (AF-Q15116-F1 and AF-Q9NZQ7-F1, model date 2025-08-01). SHA-1 of the AlphaFold files used: PD-1 `8c23c10098430fc80e89287adab0aeecf7d8096b`, PD-L1 `068df4065f2f49baade0e4f7f7a7225728c69ebc`. UniProt release 2026_02 (effective 10 June 2026); records pulled on 1 September 2026 (taken from the modification time of the saved files; the release number was not recorded at the time, so 2026_02 is inferred from the release dates). TODO: confirm the 2026_02 / 2026_03 release dates on UniProt's release page. Entry versions at pull: PD-1 (Q15116) 207, PD-L1 (Q9NZQ7) 206. RCSB/PDBe queries were run in October 2026. Database content can change, so re-running later may differ slightly.

## Expected results (use these to check a re-run)

- `python src/pd1_pdl1_revision_checks.py selftest` prints `selftest OK`. `mapping` prints UniProt Q15116 residue 78 = K and chain B author residue 78 = LYS in 4ZQK.
- `search` returns exactly two entries, 4ZQK and 5IUS. The first line of `epitope` output gives 18 residues for the 4ZQK PD-L1 interface at 4.5 A.
- `epitope` at 4.5 A: nivolumab 14 residues; pembrolizumab 20 (consensus) to 29 (union), per copy 27, 27, 23, 22; atezolizumab 24 to 30, per copy 27, 25, 26, 26, 25; durvalumab 17; avelumab 20. The three anti-PD-L1 antibodies (atezolizumab, durvalumab, avelumab) share ten residues: 56, 58-61, 113, 115, 121, 123 and 125; seven of them (56, 58, 113, 115, 121, 123, 125) are in the 4ZQK PD-1-binding interface. Consensus epitopes overlap the 18-residue 4ZQK PD-L1 interface by 13 of 24 residues (atezolizumab), 11 of 17 (durvalumab) and 11 of 20 (avelumab). TODO: confirm these against your `epitope` output.
- `src/phase4_span_calculator.py`: PD-1 alignment on 105 identical CA pairs, RMSD 0.90 A; PD-L1 alignment on 115 identical CA pairs, RMSD 0.75 A. It prints a warning because the membrane normals are 47.0 degrees from antiparallel (script threshold about 20 degrees).
- `src/span_robustness.py` on the full hybrid (saved baseline: `expected/span_robustness.txt`; `make geometry` should reproduce it): angle between the ectodomain-pointing normals 133.0 degrees (47.0 degrees from antiparallel); spans 15.22 nm (along the PD-1 normal), 17.39 nm (PD-L1 normal), 18.10 nm (direct), 17.78 nm (common normal); 195 heavy-atom pairs closer than 2.0 A between the two AlphaFold chains. On the trimmed hybrid the spans are unchanged and the pairs fall to 21. Tilting each normal by 30 degrees gives 8.4-18.0 nm along the PD-1 normal, 12.6-17.5 nm along the PD-L1 normal and 13.4-18.1 nm along the common normal; the smallest mismatch reachable with each normal tilted by at most 30 degrees is 0.3 degrees.
- `src/feasible_span.py` (deterministic: two runs gave identical output files). With the membrane surface 2.0 nm from the midplane, 3.5 A per residue and a rigid complex (10,000 directions, seed 1), 71.7% of membrane-normal directions are feasible; the lowest minimum span is 8.1 nm (median 9.8 nm) and the highest maximum span is 25.0 nm (median 20.0 nm). "Low" and "high" are the extremes over feasible directions, not confidence limits. At the three AlphaFold-derived directions the feasible ranges are 12.6-23.2 nm (PD-1 normal), 13.1-25.0 nm (minus the PD-L1 normal) and 13.5-24.6 nm (common normal). The full grid (half-bilayer 1.5-2.5 nm, 3.3-3.8 A per residue, hinge 0-40 degrees) is in `outputs/feasible_span.txt`.
- Crystal-core hybrid, `outputs/clash_summary_console.txt` (cross-chain heavy-atom pairs, PD-1 chain P against PD-L1 chain L): 61 pairs closer than 3.0 A. Crystal-crystal: none below 2.5 A and 13 between 2.5 and 3.0 A, consistent with hydrogen-bond-length contacts at the crystal interface. Crystal-AlphaFold: 15 pairs below 2.0 A, 29 below 2.5 A, 48 below 3.0 A. AlphaFold-AlphaFold: none. Excluding PD-1 residues 85-92 (the AlphaFold loop, unresolved in the crystal), there are no crystal-AlphaFold or AlphaFold-AlphaFold pairs below 3.0 A, so every clash below 2.0 A involves that loop: PD-1 residue 88 with PD-L1 residue 18 (0.62 A) and with residue 19 (1.27 A). The `vdw_flag` column in `src/clash_table.py` marks a heavy-atom pair when (r_i + r_j) - d is at least 0.4 A, with radii C 1.70, N 1.55, O 1.52 and S 1.80 A and a 0.6 A allowance when both atoms are N or O; it flags 41 pairs (1 crystal-crystal, 40 crystal-AlphaFold). ChimeraX, run on the single model with `clashes /P restrict /L overlapCutoff 0.6 hbondAllowance 0.4`, reports 40 clashes. The two lists share 39 pairs and differ in three: (a) ChimeraX only: PD-1 Gln88 NE2 with PD-L1 Thr20 N (2.19 A). Both atoms are amide donors, so ChimeraX applies no hydrogen-bond allowance (overlap 1.087 A), whereas the script subtracts its blanket 0.6 A N/O allowance (0.91 A becomes 0.31 A) and does not flag it. (b) Script only: PD-1 Gln88 C with PD-L1 Ala18 CB (2.95 A) and PD-1 Gln88 CD with PD-L1 Thr20 CG2 (3.00 A). The script's overlaps are 0.45 and 0.40 A; ChimeraX's overlaps for the same pairs are 0.540 and 0.491 A, below its 0.6 A cutoff, so both programs see an overlap and only the script's lower cutoff flags it. For these two pairs ChimeraX's implied radius sum is 3.49 A against the script's 3.40 A; other atom types differ more (3.76 A for Gln88 CA with Ala18 CB). Apart from the PD-1 Glu136 / PD-L1 Arg113 contact (2.74 A, a salt-bridge-length polar contact flagged by both programs), every flagged pair involves PD-1 residues 87-88 in the AlphaFold loop, against PD-L1 residues 18, 19, 20 and 119. The ChimeraX list is in `outputs/chimerax_clashes.txt`.
- `make test` runs 17 tests; all passed on 4 November 2026.
- `interface` for 5IUS: pairs A-C and B-D are the real contacts. PD-L1 interface Jaccard 0.80 with 4ZQK in both; PD-1 interface Jaccard 0.42 (A-C) and 0.46 (B-D).
- `apo` (3RRQ against 4ZQK PD-1): 105 matched residues, all-atom fit 1.52 A, core fit 0.63 A, interface 2.93 A under the core fit; residues 74 and 75 deviate by 8.1 A and residue 76 by 5.1 A.

## Scripts

| File (in `src/`) | Purpose |
|---|---|
| `fetch_inputs.py` | Pinned downloads, manifest, verification, linking, comparison with earlier files |
| `phase2_domain_mapping.py` | UniProt domain architecture (writes `domain_architecture.txt` and the FASTA files) |
| `phase3_manual_orient.py`, `phase4_span_calculator.py` | Membrane orientation from the TM helix axis (outputs are named `*_ppm_oriented.pdb`, but the PPM server was not used); grafting of full-length AlphaFold chains onto 4ZQK; span estimates |
| `span_robustness.py` | Normals with a defined direction, tilt sensitivity, clash check, pLDDT |
| `trim_hybrid.py` | Removes signal peptides and cytoplasmic tails from the hybrid |
| `geometry.py` | Shared geometry functions, covered by `tests/test_geometry.py` and `tests/test_cross_script.py` |
| `build_hybrid.py`, `clash_table.py`, `af_vs_crystal.py` | Crystal-core hybrid, per-pair clash table, AlphaFold-versus-crystal deviation analysis |
| `region_confidence.py`, `membrane_check.py`, `glycan_screen.py` | Confidence summary, membrane check, glycan screen |
| `feasible_span.py` | Feasible midplane-to-midplane span for a rigid PD-1:PD-L1 ectodomain complex between parallel membranes, scanning membrane-normal directions uniformly and requiring both stalks to reach their membranes (assumptions in `docs/geometry_note.md`) |
| `figure_clashes.cxc` | ChimeraX figure of PD-1/PD-L1 heavy-atom clashes in the crystal-core hybrid. Start ChimeraX, run `cd <repository folder>` in its command line, then `open src/figure_clashes.cxc` |
| `pd1_pdl1_revision_checks.py` (with `data/specs.json`) | Metadata, RCSB search, numbering, epitopes with antibody pairing, 5IUS interface, bound versus unbound PD-1, PD-1 loop ensemble |
| `render_figure2.cxc` | Figure 2 |
| `item12_epitope_crosscheck.cxc`, `compare_chimerax_results.py` | Independent epitope cross-check in ChimeraX |
| `item13_3rrq_crystal_contacts.cxc`, `bfactor_normalised_71_76.py` | Loop 71 to 76 checks |
| `item10_provenance.sh` | Records versions and compares UniProt sequences |

TODO: add the expected results of `af_vs_crystal.py`, `region_confidence.py`, `membrane_check.py` and `glycan_screen.py` from `outputs/`.

## Limitations

- The membrane normal is taken from the TM helix axis. Tilting it by up to 30 degrees changes the span along the PD-1 normal to 8.4 to 18.0 nm, so the 15 to 18 nm bracket is not a constrained result.
- The 4ZQK-based hybrid places the two membrane normals 47 degrees from antiparallel, and `phase4_span_calculator.py` warns about this. Parallel membranes are reachable by tilting each normal by at most 30 degrees (smallest mismatch 0.3 degrees), but the tilt needed is not constrained by the data.
- The hybrid model contains single static AlphaFold stalks and no lipids or glycans (N-glycosylation was not modelled). The stalks are low confidence (mean pLDDT 62 for PD-1 residues 146-170 and 67 for PD-L1 residues 226-238, against 93-96 for the Ig domains), and the span depends on them.
- In the AlphaFold-derived hybrid (`phase4`), after trimming, 21 Ig-domain heavy-atom pairs remain closer than 2.0 A (closest 0.98 A; PD-1 residues 88, 131 and 132 against PD-L1 residues 19, 56, 63 and 66). They are not a superposition error: in the crystal the corresponding residue pairs are 3.4 to 6.6 A apart. AlphaFold places PD-1 residues 131 and 132 differently from the crystal (Calpha 2.9 and 1.5 A) and places the 85-92 loop, which is unresolved in 4ZQK, over PD-L1 residues 18-19. That hybrid is therefore not steric-clean, which is why the crystal-core hybrid (`src/build_hybrid.py`) takes the Ig domains from the crystal structure and uses AlphaFold only for the stalks, TM helices and the 85-92 loop. The spans depend only on the TM helices and are unchanged by trimming.
- PD-1 residues 85-92 are unresolved in 4ZQK, so contacts in that segment (including part of the pembrolizumab epitope) have no counterpart in 4ZQK.
- The second human complex, 5IUS, has an engineered high-affinity PD-1 and lower resolution (2.89 A). Its PD-L1 interface agrees closely with 4ZQK (16 of 18 residues shared) but its PD-1 interface less so (13 of 21), so the PD-1 interface is less well established.
- Epitope sizes vary between copies and with the cutoff. Contacts with symmetry-related molecules outside the asymmetric unit were not examined, and 5X8L (3.10 A) and 5GRJ (3.21 A) are low resolution.
- The feasible-span model assumes a rigid ectodomain complex, parallel membranes, a uniform prior over membrane-normal directions and fully extended stalks for the maximum, and it ignores glycans and lipid detail. Its 8-25 nm range is a geometric feasibility bound, not a prediction of the actual separation; the 15.2-18.1 nm hybrid values lie inside it.
- Epitopes are defined by a heavy-atom distance cutoff, not by energy.
- There is no molecular dynamics: all contacts come from static structures, so no contact frequencies over time are reported. Comparisons across structures (Jaccard, RMSD) are descriptive and carry no uncertainty estimates.
- TODO: the next paragraph names specific citations and values, but the manuscript still lists the literature spans and their measurement conventions as unverified (Cordoba et al. 2013; Strazza et al. 2024), and gives the CD45 height as about 25 nm for the smallest isoform. Verify each citation and value, then reconcile the CD45 figure with the 28-50 nm stated here, or soften the paragraph.
- The computed PD-1:PD-L1 span (15.2-18.1 nm; 8.4-18.0 nm under tilt sensitivity) is of the same order as the ~10-15 nm extracellular span widely cited for immune receptor-ligand pairs, originating in crystal structures of KIR2DL1-HLA-Cw4 and NKG2D-MICA (Fan et al. 2001; Li et al. 2001; Boyington et al. 2000) and propagated through subsequent reviews (Alakoskela et al. 2011; Strazza et al. 2024). The exact measurement convention underlying the original ~15 nm figure (e.g., whether end-to-end or membrane-normal projected length) is not stated in the reviewing literature and was not independently re-derived here; the comparison is therefore qualitative, not a matched metric. The CD45 ectodomain height (~28-50 nm across splice isoforms; Woollett et al. 1985, McCall et al. 1992, as cited in Cordoba et al. 2013) is well above this range, consistent with CD45's size-based exclusion from close-contact zones under the kinetic-segregation model.

## Open items

- Verify the literature span values and their measurement conventions against the original papers.
- Independent epitope cross-check in ChimeraX (`src/item12_epitope_crosscheck.cxc`) and crystal contacts around PD-1 residues 71-76 in 3RRQ with chain-normalised B-factors (`src/item13_3rrq_crystal_contacts.cxc`). TODO: record whether these have been run and what they gave.
- Resolve the placement of the AlphaFold 85-92 loop in the crystal-core hybrid, which is the only source of clashes below 2.0 A (against PD-L1 residues 18-19, and ChimeraX also flags contact with PD-L1 residue 119).
- Record the UniProt release number directly.

## License

MIT (add a `LICENSE` file).
