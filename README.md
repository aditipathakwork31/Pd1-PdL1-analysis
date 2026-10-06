# PD-1:PD-L1 geometric model and therapeutic antibody epitope mapping

Analysis scripts for the manuscript "A Membrane-Anchored Geometric Model of PD-1:PD-L1 and Therapeutic Antibody Epitope Mapping".

## Environment

Python 3.14.7, Biopython 1.88, NumPy 2.5.3, requests 2.34.2, Matplotlib 3.11.2. UCSF ChimeraX: version 1.13.dev202609011810 (2026-09-01).

AlphaFold models: database version 6 (AF-Q15116-F1 and AF-Q9NZQ7-F1, model date 2025-08-01).

Download the AlphaFold files separately. SHA-1 of the files used: PD-1 8c23c10098430fc80e89287adab0aeecf7d8096b, PD-L1 068df4065f2f49baade0e4f7f7a7225728c69ebc.

## Scripts

* phase2_domain_mapping.py: UniProt domain architecture (writes domain_architecture.txt and the FASTA files)
* phase3_manual_orient.py, phase4_span_calculator.py, phase6_span_chart.py: membrane orientation and span
* phase5_epitope_mapper.py: original epitope mapper
* pd1_pdl1_revision_checks.py (with specs.json): metadata, RCSB search, numbering, epitopes, 5IUS interface, bound versus unbound PD-1
* span_robustness.py: normals, tilt sensitivity, clashes, pLDDT
* trim_hybrid.py: removes signal peptides and cytoplasmic tails from the hybrid
* render_figure2.cxc: Figure 2 in ChimeraX
* item12_epitope_crosscheck.cxc and compare_chimerax_results.py: independent epitope cross-check in ChimeraX
* item13_3rrq_crystal_contacts.cxc and bfactor_normalised_71_76.py: loop 71-76 checks
* item10_provenance.sh: records versions and compares UniProt sequences

PDB files are downloaded from the RCSB at run time and are not included.

