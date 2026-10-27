#!/usr/bin/env python3
"""
Re-analysis checks for the PD-1:PD-L1 membrane-geometry / epitope study.

Steps (run any subset):
  meta      Fetch resolution, title, organism, chains for each PDB entry (RCSB REST).
  mapping   Show PDB-author <-> UniProt residue mapping (PDBe SIFTS) and check a residue
            (e.g. "Lys78") against the UniProt sequence AND the PDB numbering.
  epitope   Epitope residues with (a) cognate-copy filtering, (b) cutoff sensitivity,
            (c) overlap with the 4ZQK PD-L1 interface. Output in UniProt numbering.
  apo       Calpha RMSD of unbound human PD-1 (3RRQ) vs PD-1 in 4ZQK, matched by UniProt number.
  versions  Print software versions for the Methods section.

Requires: biopython, numpy, requests. Needs internet access to rcsb.org, ebi.ac.uk, uniprot.org.
NOTE: written for the revision of the manuscript; run it and record the output.
Unit-tested only on synthetic geometry (see __main__ --selftest), not yet on live data.
"""
import argparse, io, sys, json, itertools
from collections import defaultdict

import numpy as np
import requests
from Bio.PDB import MMCIFParser, NeighborSearch, Superimposer
from Bio.PDB.Polypeptide import is_aa

PDB_IDS = ["4ZQK", "3RRQ", "5B8C", "5WT9", "5X8L", "5X8M", "5GRJ"]
UNIPROT = {"PD1": "Q15116", "PDL1": "Q9NZQ7"}
CUTOFFS = (4.0, 4.5, 5.0)
MIN_PAIR_CONTACTS = 20   # atom-atom contacts (<= cutoff) for an antibody chain to count as cognate


def get(url, **kw):
    r = requests.get(url, timeout=60, **kw)
    r.raise_for_status()
    return r


# ---------------------------------------------------------------- meta
def step_meta(ids=None):
    """Resolution, title, and per-entity description / organism / UniProt accession / chains for each PDB id."""
    ids = [i.upper() for i in (ids or PDB_IDS)]
    print("PDB\tres(A)\tentities (description [organism] UniProt chains)\ttitle")
    for pid in ids:
        e = get(f"https://data.rcsb.org/rest/v1/core/entry/{pid}").json()
        res = e.get("rcsb_entry_info", {}).get("resolution_combined", ["NA"])[0]
        title = e["struct"]["title"]
        n_ent = int(e["rcsb_entry_info"]["polymer_entity_count"])
        ents = []
        for i in range(1, n_ent + 1):
            pe = get(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pid}/{i}").json()
            org = ",".join(o.get("ncbi_scientific_name", "?") for o in pe.get("rcsb_entity_source_organism", []))
            cid = pe["rcsb_polymer_entity_container_identifiers"]
            ch = ",".join(cid["auth_asym_ids"])
            acc = ",".join(r.get("database_accession", "?") for r in cid.get("reference_sequence_identifiers", [])
                           if r.get("database_name") == "UniProt") or "-"
            ents.append(f"{pe['rcsb_polymer_entity']['pdbx_description']} [{org}] {acc} chains {ch}")
        print(f"{pid}\t{res}\t" + " | ".join(ents) + f"\t{title}")


def rcsb_search_both(acc1="Q15116", acc2="Q9NZQ7", rows=100):
    """All PDB entries that contain BOTH UniProt accessions (RCSB Search API v2, entry level). HTTP 204 = no hits."""
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    attr = "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession"
    node = lambda acc: {"type": "terminal", "service": "text",
                        "parameters": {"attribute": attr, "operator": "exact_match", "value": acc}}
    ids, start = [], 0
    while True:
        body = {"query": {"type": "group", "logical_operator": "and", "nodes": [node(acc1), node(acc2)]},
                "return_type": "entry", "request_options": {"paginate": {"start": start, "rows": rows}}}
        r = requests.post(url, json=body, timeout=60)
        if r.status_code == 204:
            break
        r.raise_for_status()
        d = r.json()
        ids += [x["identifier"] for x in d.get("result_set", [])]
        start += rows
        if start >= d.get("total_count", 0):
            break
    return ids


def step_search():
    """Definitive list of PDB entries containing both human PD-1 (Q15116) and human PD-L1 (Q9NZQ7)."""
    ids = rcsb_search_both(UNIPROT["PD1"], UNIPROT["PDL1"])
    print(f"RCSB entries containing both {UNIPROT['PD1']} and {UNIPROT['PDL1']}: {len(ids)}")
    print(" ".join(ids))
    if ids:
        print("\nDetails (review each: complex type, engineered mutants, antibody or small-molecule partners):")
        step_meta(ids)


# ---------------------------------------------------------------- numbering
def sifts_segments(pid, uniprot):
    """Return {chain: [(auth_start, auth_end, unp_start, unp_end), ...]} for one UniProt accession."""
    d = get(f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pid.lower()}").json()[pid.lower()]
    segs = defaultdict(list)
    m = d["UniProt"].get(uniprot)
    if not m:
        return segs
    for s in m["mappings"]:
        a0, a1 = s["start"]["author_residue_number"], s["end"]["author_residue_number"]
        u0, u1 = s["unp_start"], s["unp_end"]
        # SIFTS returns null author numbers where the end (or start) residues are unobserved
        # (e.g. 4ZQK PD-1 residues 147-150). Infer them from the constant author->UniProt offset.
        if a1 is None and a0 is not None:
            a1 = a0 + (u1 - u0)
        if a0 is None and a1 is not None:
            a0 = a1 - (u1 - u0)
        if a0 is None or a1 is None:
            continue  # cannot place this segment; its residues stay unmapped (excluded, not guessed)
        segs[s["chain_id"]].append((a0, a1, u0, u1))
    return segs


def auth_to_unp(segs_chain, resnum):
    for a0, a1, u0, u1 in segs_chain:
        if a0 <= resnum <= a1 and (a1 - a0) == (u1 - u0):
            return u0 + (resnum - a0)
    return None  # not mappable (insertion codes / engineered residues / gaps)


def uniprot_seq(acc):
    fa = get(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta").text.splitlines()
    return "".join(fa[1:])


def step_mapping(pdb="4ZQK", residue=78):
    one3 = {"K": "LYS", "R": "ARG"}
    seq = uniprot_seq(UNIPROT["PD1"])
    print(f"UniProt {UNIPROT['PD1']} residue {residue} = {seq[residue-1]}")
    print("Neighbourhood (UniProt numbering):",
          " ".join(f"{i}{seq[i-1]}" for i in range(residue - 8, residue + 9)))
    segs = sifts_segments(pdb, UNIPROT["PD1"])
    print(f"SIFTS segments for {pdb} (chain: auth_start-auth_end -> UniProt start-end):")
    for ch, lst in segs.items():
        print("  ", ch, lst)
    cif = MMCIFParser(QUIET=True).get_structure(pdb, io.StringIO(get(f"https://files.rcsb.org/download/{pdb}.cif").text))
    for ch, lst in segs.items():
        for res in cif[0][ch]:
            if res.id[0] == " " and auth_to_unp(lst, res.id[1]) == residue:
                print(f"{pdb} chain {ch}: UniProt {residue} == author resnum {res.id[1]} ({res.get_resname()})")
    print("Decision rule: if UniProt residue is not Lys, the paper's 'Lys78' uses a different numbering; "
          "find the Lys in the neighbourhood above and record the offset.")


# ---------------------------------------------------------------- epitopes
def load(pid):
    return MMCIFParser(QUIET=True).get_structure(pid, io.StringIO(get(f"https://files.rcsb.org/download/{pid}.cif").text))[0]


def heavy_atoms(chain):
    return [a for r in chain if r.id[0] == " " and is_aa(r, standard=True) for a in r if a.element != "H"]


def pair_contacts(ag_chain, ab_chain, cutoff):
    ns = NeighborSearch(heavy_atoms(ab_chain))
    n = 0
    for a in heavy_atoms(ag_chain):
        n += len(ns.search(a.coord, cutoff))
    return n


def epitope_residues(ag_chain, ab_chains, cutoff):
    ns = NeighborSearch([a for c in ab_chains for a in heavy_atoms(c)])
    out = set()
    for r in ag_chain:
        if r.id[0] != " " or not is_aa(r, standard=True):
            continue
        if any(ns.search(a.coord, cutoff) for a in r if a.element != "H"):
            out.add(r.id[1])
    return out


def cognate_epitope(model, ag_chain_ids, ab_chain_ids, cutoff, min_contacts=MIN_PAIR_CONTACTS):
    """For each antigen chain, use only antibody chains with >= min_contacts atom contacts to it
    (the cognate Fab chains), so contacts with other copies in the asymmetric unit are not counted."""
    result = {}
    for ag in ag_chain_ids:
        cog = [model[ab] for ab in ab_chain_ids if pair_contacts(model[ag], model[ab], cutoff) >= min_contacts]
        result[ag] = (sorted(c.id for c in cog), epitope_residues(model[ag], cog, cutoff) if cog else set())
    return result


def consensus_unp(per_copy, segs):
    """Convert each copy's author-numbered epitope to UniProt numbering; return per-copy sets and the
    union / intersection across copies."""
    sets = []
    for ag, (cog, res) in per_copy.items():
        unp = {auth_to_unp(segs[ag], r) for r in res}
        unp.discard(None)
        sets.append(unp)
    if not sets:
        return [], set(), set()
    return sets, set().union(*sets), set.intersection(*sets)



def fallback_segments(chain, unp_seq, min_identity=0.8):
    """Sequence-alignment mapping (author number -> UniProt number) for chains SIFTS cannot map.
    Returns segments in the same (auth_start, auth_end, unp_start, unp_end) format, plus the identity."""
    from Bio.Align import PairwiseAligner
    from Bio.Data.IUPACData import protein_letters_3to1
    res = [r for r in chain if r.id[0] == " " and r.id[2] == " " and is_aa(r, standard=True)]
    seq = "".join(protein_letters_3to1[r.get_resname().capitalize()] for r in res)
    al = PairwiseAligner()
    al.mode = "local"; al.match_score = 2; al.mismatch_score = -1
    al.open_gap_score = -5; al.extend_gap_score = -1
    aln = al.align(unp_seq, seq)[0]
    pairs = []   # (author number, UniProt number)
    for (u0, u1), (s0, s1) in zip(*aln.aligned):
        for k in range(u1 - u0):
            pairs.append((int(res[s0 + k].id[1]), int(u0 + k + 1), unp_seq[u0 + k] == seq[s0 + k]))
    if not pairs:
        raise ValueError("sequence alignment found no matching region")
    identity = sum(p[2] for p in pairs) / len(pairs)
    if identity < min_identity:
        raise ValueError(f"alignment identity {identity:.2f} below {min_identity}")
    segs, run = [], []
    for a, u, _ in pairs:
        if run and a == run[-1][0] + 1 and u == run[-1][1] + 1:
            run.append((a, u))
        else:
            if run: segs.append((run[0][0], run[-1][0], run[0][1], run[-1][1]))
            run = [(a, u)]
    if run: segs.append((run[0][0], run[-1][0], run[0][1], run[-1][1]))
    return segs, identity


def segments_for(pid, model, acc, chains):
    """SIFTS segments per chain; where SIFTS has none, fall back to sequence alignment (and say so)."""
    segs = dict(sifts_segments(pid, acc))
    seq = None
    for ch in chains:
        if not segs.get(ch):
            seq = seq or uniprot_seq(acc)
            segs[ch], ident = fallback_segments(model[ch], seq)
            print(f"   [{pid} chain {ch}: no SIFTS mapping; used sequence alignment to {acc}, identity {ident:.2f}]")
    return segs


def pair_fabs(model, heavy, light, cutoff):
    """One-to-one heavy<->light pairing by heavy-atom contacts (greedy, strongest first)."""
    scored = sorted(((pair_contacts(model[h], model[l], cutoff), h, l) for h in heavy for l in light), reverse=True)
    used_h, used_l, groups = set(), set(), []
    for n, h, l in scored:
        if h not in used_h and l not in used_l and n > 0:
            groups.append((h, l)); used_h.add(h); used_l.add(l)
    return groups


def assign_antigens(model, ag_chains, groups, cutoff):
    """One-to-one antigen copy <-> Fab/Fv group by total contacts (greedy, strongest first)."""
    scored = sorted(((pair_contacts(model[ag], model[h], cutoff) + pair_contacts(model[ag], model[l], cutoff), ag, (h, l))
                     for ag in ag_chains for (h, l) in groups), reverse=True)
    used_ag, used_g, out = set(), set(), {}
    for n, ag, g in scored:
        if ag not in used_ag and g not in used_g and n > 0:
            out[ag] = g; used_ag.add(ag); used_g.add(g)
    return out


def paired_epitopes(model, sp, cutoff):
    """Epitope per antigen copy using ONLY its own paired antibody chains. Returns
    {antigen chain: (group, residues, contacts_with_group, contacts_with_other_antibody_chains)}."""
    groups = pair_fabs(model, sp["heavy"], sp["light"], cutoff)
    pairing = assign_antigens(model, sp["ag_chains"], groups, cutoff)
    out = {}
    for ag in sp["ag_chains"]:
        if ag not in pairing:
            out[ag] = (None, set(), 0, 0); continue
        h, l = pairing[ag]
        own = pair_contacts(model[ag], model[h], cutoff) + pair_contacts(model[ag], model[l], cutoff)
        other = sum(pair_contacts(model[ag], model[c], cutoff) for c in sp["heavy"] + sp["light"] if c not in (h, l))
        out[ag] = ((h, l), epitope_residues(model[ag], [model[h], model[l]], cutoff), own, other)
    return out


def step_epitope(specs):
    """specs: JSON file, per antibody structure: {"antigen": "PD1"|"PDL1", "ag_chains": [...], "heavy": [...],
    "light": [...], "ab_chains": [...]}. Chain IDs MUST be confirmed from the entity records (run `meta`).
    Each antigen copy is paired one-to-one with a heavy/light group by contact count; only that group's
    chains define its epitope, so contacts with neighbouring copies (crystal packing) are excluded."""
    spec = json.load(open(specs))
    # PD-L1 interface in 4ZQK (UniProt numbering) for overlap
    m4 = load("4ZQK")
    pd1_ch, pdl1_ch = spec["4ZQK"]["PD1_chain"], spec["4ZQK"]["PDL1_chain"]
    seg_pdl1 = sifts_segments("4ZQK", UNIPROT["PDL1"])[pdl1_ch]
    iface = {auth_to_unp(seg_pdl1, r) for r in epitope_residues(m4[pdl1_ch], [m4[pd1_ch]], 4.5)} - {None}
    print(f"4ZQK PD-L1 interface residues (UniProt, 4.5 A): n={len(iface)}")
    for pid, sp in spec.items():
        if pid == "4ZQK":
            continue
        model = load(pid)
        segs = segments_for(pid, model, UNIPROT[sp["antigen"]], sp["ag_chains"])
        for cutoff in CUTOFFS:
            per = paired_epitopes(model, sp, cutoff)
            sets = []
            for ag, (grp, res, own, other) in per.items():
                unp = {auth_to_unp(segs[ag], r) for r in res}; unp.discard(None)
                sets.append(unp)
                if cutoff == 4.5:
                    print(f"   {pid} antigen {ag} <-> (heavy {grp[0]}, light {grp[1]})" if grp else f"   {pid} antigen {ag}: NO PAIRING FOUND",
                          f"| own contacts {own}, contacts with other antibody chains (packing) {other}")
            union = set().union(*sets) if sets else set()
            inter = set.intersection(*sets) if sets else set()
            line = f"{pid} cutoff {cutoff}: per-copy sizes {[len(x) for x in sets]}; union {len(union)}; consensus {len(inter)}"
            if sp["antigen"] == "PDL1" and inter:
                ov = inter & iface
                line += f"; overlap with 4ZQK interface {len(ov)}/{len(inter)} ({100*len(ov)/len(inter):.0f}%)"
            print(line)
            if cutoff == 4.5:
                print("   consensus residues (UniProt):", sorted(inter))
                print("   union residues (UniProt):", sorted(union))


# ---------------------------------------------------------------- diagnose
def step_diagnose(specs, only=None):
    """Print, for each antibody structure: chains present, SIFTS segments, and the antigen x antibody
    heavy-atom contact-count matrix at 4.5 A. Use it to see which antibody chains truly pair with which
    antigen copy, and to find why an epitope came out empty."""
    spec = json.load(open(specs))
    for pid, sp in spec.items():
        if pid == "4ZQK" or (only and pid != only):
            continue
        model = load(pid)
        print(f"===== {pid} (antigen {sp['antigen']})")
        print("chains in file:", {c.id: sum(1 for r in c if r.id[0] == ' ') for c in model})
        segs = dict(sifts_segments(pid, UNIPROT[sp["antigen"]]))
        print("SIFTS segments (auth_start, auth_end, unp_start, unp_end):", segs)
        missing = [c for c in sp["ag_chains"] + sp["ab_chains"] if c not in model]
        if missing:
            print("!! chain IDs in specs.json not found in file:", missing)
            continue
        print("contact counts (rows = antigen chain, columns = antibody chain):")
        print("      " + "".join(f"{ab:>6}" for ab in sp["ab_chains"]))
        for ag in sp["ag_chains"]:
            row = [pair_contacts(model[ag], model[ab], 4.5) for ab in sp["ab_chains"]]
            print(f"  {ag:>3} " + "".join(f"{n:>6}" for n in row))
        for ag in sp["ag_chains"]:
            raw = sorted(epitope_residues(model[ag], [model[ab] for ab in sp["ab_chains"]], 4.5))
            mapped = [auth_to_unp(segs.get(ag, []), r) for r in raw]
            print(f"  antigen {ag}: {len(raw)} raw author-numbered contact residues (any antibody chain); "
                  f"{sum(m is not None for m in mapped)} map to UniProt; first residues {raw[:12]}")

# ---------------------------------------------------------------- apo vs bound
def kabsch(mob, ref):
    """Least-squares superposition of mob onto ref (both Nx3). Returns (aligned mob, rmsd)."""
    mc, rc = mob.mean(0), ref.mean(0)
    H = (mob - mc).T @ (ref - rc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    out = (mob - mc) @ R.T + rc
    return out, float(np.sqrt(((out - ref) ** 2).sum(1).mean()))


def step_apo(bound_chain="B", apo_chain="A", core_cut=2.0):
    """Unbound human PD-1 (3RRQ) vs PD-1 in the 4ZQK complex, residues matched by UniProt number.
    Reports Calpha RMSD for all matched residues, for the 4ZQK interface residues, and for a core fit
    (residues within core_cut A after a first fit), plus the residues that move most."""
    bound, apo = load("4ZQK"), load("3RRQ")
    seg_b = segments_for("4ZQK", bound, UNIPROT["PD1"], [bound_chain])
    seg_a = segments_for("3RRQ", apo, UNIPROT["PD1"], [apo_chain])
    def ca_by_unp(model, ch, segs):
        d = {}
        for r in model[ch]:
            if r.id[0] == " " and "CA" in r:
                u = auth_to_unp(segs[ch], r.id[1])
                if u is not None:
                    d[u] = r["CA"].coord
        return d
    cb, ca = ca_by_unp(bound, bound_chain, seg_b), ca_by_unp(apo, apo_chain, seg_a)
    common = sorted(set(cb) & set(ca))
    print(f"4ZQK chain {bound_chain}: {len(cb)} residues with CA; 3RRQ chain {apo_chain}: {len(ca)}; matched by UniProt number: {len(common)}")
    rng = range(min(min(cb), min(ca)), max(max(cb), max(ca)) + 1)
    print("  residues in 33-150 missing in 4ZQK:", [u for u in range(33, 151) if u not in cb])
    print("  residues in 33-150 missing in 3RRQ:", [u for u in range(33, 151) if u not in ca])
    B = np.array([cb[u] for u in common]); A = np.array([ca[u] for u in common])
    A1, rms_all = kabsch(A, B)
    dev = np.sqrt(((A1 - B) ** 2).sum(1))
    iface = {auth_to_unp(seg_b[bound_chain], r) for r in epitope_residues(bound[bound_chain], [bound["A"]], 4.5)} - {None}
    idx_if = [i for i, u in enumerate(common) if u in iface]
    rms_if = float(np.sqrt((dev[idx_if] ** 2).mean())) if idx_if else float("nan")
    print(f"  fit on all {len(common)} matched CA: RMSD {rms_all:.2f} A")
    print(f"  interface residues present in both ({len(idx_if)} of {len(iface)}): RMSD {rms_if:.2f} A under that fit")
    core = [i for i in range(len(common)) if dev[i] < core_cut]
    A2, rms_core = kabsch(A[core], B[core])
    # apply the core fit to all residues: refit with core only, then measure everything
    mc, rc = A[core].mean(0), B[core].mean(0)
    H = (A[core] - mc).T @ (B[core] - rc); U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))]) @ U.T
    dev2 = np.sqrt((((A - mc) @ R.T + rc - B) ** 2).sum(1))
    print(f"  core fit ({len(core)} residues within {core_cut} A after the first fit): core RMSD {rms_core:.2f} A")
    if idx_if:
        print(f"  interface RMSD under core fit: {float(np.sqrt((dev2[idx_if] ** 2).mean())):.2f} A")
    far = sorted(((dev2[i], common[i]) for i in range(len(common)) if dev2[i] >= core_cut), reverse=True)
    print(f"  residues deviating >= {core_cut} A under the core fit ({len(far)}):", [(u, round(float(d), 1)) for d, u in far][:15])
    print("  note: 4ZQK has C93S and 3RRQ has A132L (engineered); positions are matched, side chains differ.")


def step_pd1ensemble(loop=range(71, 77), probe=(78,), core_cut=2.0):
    """Compare the PD-1 loop (UniProt 71-76) and Lys78 across every PD-1 copy against 4ZQK chain B,
    after a core superposition, to test whether the 3RRQ/4ZQK difference is a bound-vs-unbound
    split or just scatter between independent copies. Also prints mean B-factors of the loop."""
    ref_pid, ref_ch = "4ZQK", "B"
    entries = [("3RRQ", "A", "unbound"), ("5B8C", "C", "pembrolizumab"), ("5B8C", "F", "pembrolizumab"),
               ("5B8C", "I", "pembrolizumab"), ("5B8C", "L", "pembrolizumab"), ("5WT9", "G", "nivolumab")]
    models = {}
    def maps(pid, ch):
        if pid not in models:
            models[pid] = load(pid)
        model = models[pid]
        segs = segments_for(pid, model, UNIPROT["PD1"], [ch])
        ca, bf = {}, {}
        for r in model[ch]:
            if r.id[0] == " " and "CA" in r:
                u = auth_to_unp(segs[ch], r.id[1])
                if u is not None:
                    ca[u] = r["CA"].coord
                    bf[u] = float(np.mean([a.bfactor for a in r]))
        return ca, bf
    ref, ref_bf = maps(ref_pid, ref_ch)
    cols = list(loop) + list(probe)
    print("deviation (A) of Calpha from 4ZQK chain B after a core fit; '-' = residue not modelled")
    print(f"{'structure/chain':<22}{'core n':>7}{'coreRMSD':>9}  " + "".join(f"{c:>6}" for c in cols) + "   mean B of loop")
    print(f"{'4ZQK B (reference)':<22}{'':>7}{'':>9}  " + "".join(f"{'0.0':>6}" if c in ref else f"{'-':>6}" for c in cols)
          + f"   {np.mean([ref_bf[c] for c in loop if c in ref_bf]):.0f}")
    for pid, ch, label in entries:
        ca, bf = maps(pid, ch)
        common = sorted(set(ca) & set(ref))
        A = np.array([ca[u] for u in common]); B = np.array([ref[u] for u in common])
        A1, _ = kabsch(A, B)
        dev = np.sqrt(((A1 - B) ** 2).sum(1))
        core = [i for i in range(len(common)) if dev[i] < core_cut]
        mc, rc = A[core].mean(0), B[core].mean(0)
        H = (A[core] - mc).T @ (B[core] - rc); U, S, Vt = np.linalg.svd(H)
        R = Vt.T @ np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))]) @ U.T
        d2 = dict(zip(common, np.sqrt((((A - mc) @ R.T + rc - B) ** 2).sum(1))))
        core_rms = float(np.sqrt((np.array([d2[common[i]] for i in core]) ** 2).mean()))
        cells = "".join(f"{d2[c]:>6.1f}" if c in d2 else f"{'-':>6}" for c in cols)
        lb = [bf[c] for c in loop if c in bf]
        print(f"{pid+' '+ch+' '+label:<22}{len(core):>7}{core_rms:>9.2f}  {cells}   {np.mean(lb):.0f}" if lb else f"{pid} {ch}: loop not modelled")
    print("Read: if 3RRQ alone differs at residues 74-75 while all the antibody-bound copies sit with 4ZQK, the loop "
          "is bound-like in complexes; if the copies scatter widely, the loop is intrinsically flexible.")


def interface_sets(model, pd1_chain, pdl1_chain, pid, cutoff=4.5):
    """Interface residues of each chain (UniProt numbering) for one PD-1/PD-L1 pair."""
    sg1 = segments_for(pid, model, UNIPROT["PD1"], [pd1_chain])[pd1_chain]
    sg2 = segments_for(pid, model, UNIPROT["PDL1"], [pdl1_chain])[pdl1_chain]
    to_unp = lambda segs, res: {auth_to_unp(segs, r) for r in res} - {None}
    return (to_unp(sg1, epitope_residues(model[pd1_chain], [model[pdl1_chain]], cutoff)),
            to_unp(sg2, epitope_residues(model[pdl1_chain], [model[pd1_chain]], cutoff)))


def compare_sets(name, a, b, a_label, b_label):
    inter = a & b
    print(f"  {name}: {a_label} {len(a)}, {b_label} {len(b)}, shared {len(inter)}, Jaccard {len(inter) / len(a | b):.2f}")
    print(f"    only in {a_label}: {sorted(a - b)}")
    print(f"    only in {b_label}: {sorted(b - a)}")


def step_interface(pdb, pd1_chain, pdl1_chain, cutoff=4.5):
    """Cross-validate the 4ZQK interface against another human PD-1:PD-L1 entry (e.g. 5IUS).
    Chain letters must be read from the 'meta' output, not guessed."""
    ref = load("4ZQK"); other = load(pdb.upper())
    r1, r2 = interface_sets(ref, "B", "A", "4ZQK", cutoff)
    o1, o2 = interface_sets(other, pd1_chain, pdl1_chain, pdb.upper(), cutoff)
    print(f"Interface residues at {cutoff} A (UniProt numbering): 4ZQK vs {pdb.upper()} (PD-1 chain {pd1_chain}, PD-L1 chain {pdl1_chain})")
    compare_sets("PD-1", r1, o1, "4ZQK", pdb.upper())
    compare_sets("PD-L1", r2, o2, "4ZQK", pdb.upper())
    print(f"  Lys78 in PD-1 interface: 4ZQK {78 in r1}, {pdb.upper()} {78 in o1}")
    print(f"  PD-L1 groove residues 56/115/122/123 in interface: 4ZQK {[x in r2 for x in (56, 115, 122, 123)]}, "
          f"{pdb.upper()} {[x in o2 for x in (56, 115, 122, 123)]}")


def step_versions():
    import Bio, numpy, platform
    print("python", platform.python_version(), "| biopython", Bio.__version__, "| numpy", numpy.__version__)


# ---------------------------------------------------------------- self-test (synthetic geometry only)
def _selftest():
    from Bio.PDB import Structure, Model, Chain, Residue, Atom
    def mk(chain_id, coords):
        c = Chain.Chain(chain_id)
        for i, xyz in enumerate(coords, 1):
            r = Residue.Residue((" ", i, " "), "ALA", " ")
            for name in ("N", "CA", "C", "O"):
                r.add(Atom.Atom(name, np.array(xyz, float), 1.0, 1.0, " ", name, 0, name[0]))
            c.add(r)
        return c
    m = Model.Model(0)
    m.add(mk("A", [(0, 0, 0), (10, 0, 0), (20, 0, 0)]))            # antigen copy
    m.add(mk("H", [(0, 3, 0), (10, 3, 0)]))                        # cognate Ab: contacts res 1,2
    m.add(mk("X", [(20, 3, 0)] ))                                  # other copy's Ab: contacts res 3 only
    r = cognate_epitope(m, ["A"], ["H", "X"], 4.0, min_contacts=20)  # H: 32 contacts, X: 16
    assert r["A"][1] == {1, 2}, r          # residue 3 must NOT be called (non-cognate chain filtered out)
    assert r["A"][0] == ["H"], r
    assert auth_to_unp([(10, 20, 33, 43)], 15) == 38 and auth_to_unp([(10, 20, 33, 43)], 25) is None
    # regression test for the 4ZQK case: SIFTS segment with a null author end (unobserved C-terminus)
    import unittest.mock as um
    fake = {"4zqk": {"UniProt": {"Q15116": {"mappings": [
        {"chain_id": "B", "start": {"author_residue_number": 33}, "end": {"author_residue_number": None},
         "unp_start": 33, "unp_end": 150}]}}}}
    class R:  # minimal stand-in for a requests response
        def json(self): return fake
    with um.patch(__name__ + ".get", lambda url, **kw: R()):
        sg = sifts_segments("4ZQK", "Q15116")
    assert sg["B"] == [(33, 150, 33, 150)], sg
    assert auth_to_unp(sg["B"], 78) == 78

    # ---- Fab pairing must exclude a planted packing contact
    def mk2(chain_id, coords):
        c = Chain.Chain(chain_id)
        for i, xyz in enumerate(coords, 1):
            r = Residue.Residue((" ", i, " "), "ALA", " ")
            for name in ("N", "CA", "C", "O"):
                r.add(Atom.Atom(name, np.array(xyz, float), 1.0, 1.0, " ", name, 0, name[0]))
            c.add(r)
        return c
    m2 = Model.Model(0)
    m2.add(mk2("A", [(1.5, 3, 0), (1.5, 3, 20)]))            # antigen copy 1: res1 near Fab1, res2 only near L2's extra residue
    m2.add(mk2("B", [(1.5, 103, 0)]))                        # antigen copy 2
    m2.add(mk2("H", [(0, 0, 0)]));  m2.add(mk2("L", [(3, 0, 0)]))          # Fab 1
    m2.add(mk2("X", [(0, 100, 0)])); m2.add(mk2("Y", [(3, 100, 0), (1.5, 6.5, 20)]))  # Fab 2 (Y has a packing residue near A res2)
    sp = {"ag_chains": ["A", "B"], "heavy": ["H", "X"], "light": ["L", "Y"]}
    pe = paired_epitopes(m2, sp, 4.5)
    assert pe["A"][0] == ("H", "L") and pe["B"][0] == ("X", "Y"), pe
    assert pe["A"][1] == {1}, pe["A"]          # residue 2 touches only the NON-cognate chain Y -> excluded
    assert pe["A"][3] > 0                      # and the packing contact is reported
    # ---- sequence-alignment fallback mapping
    unp = "ACDEFGHIKLMNPQRSTVWY"
    names = {"A":"ALA","C":"CYS","D":"ASP","E":"GLU","F":"PHE","G":"GLY","H":"HIS","I":"ILE","K":"LYS","L":"LEU","M":"MET","N":"ASN"}
    ch = Chain.Chain("Z")
    for k, aa in enumerate(unp[2:12]):                       # UniProt residues 3..12 -> author 103..112
        ch.add(Residue.Residue((" ", 103 + k, " "), names[aa], " "))
    sg, ident = fallback_segments(ch, unp)
    assert auth_to_unp(sg, 103) == 3 and auth_to_unp(sg, 112) == 12 and ident == 1.0, (sg, ident)

    # ---- Kabsch: a rotated + translated copy must superpose with ~0 RMSD
    rng = np.random.default_rng(1)
    P = rng.normal(size=(30, 3))
    th = 0.7
    Rz = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
    _, rms = kabsch(P @ Rz.T + np.array([5.0, -2.0, 3.0]), P)
    assert rms < 1e-8, rms

    # ---- RCSB search paging and the empty (204) response, with the network mocked
    import unittest.mock as um
    class Resp:
        def __init__(self, code, data=None): self.status_code, self._d = code, data
        def json(self): return self._d
        def raise_for_status(self): pass
    pages = [Resp(200, {"total_count": 3, "result_set": [{"identifier": "AAAA"}, {"identifier": "BBBB"}]}),
             Resp(200, {"total_count": 3, "result_set": [{"identifier": "CCCC"}]})]
    with um.patch.object(requests, "post", side_effect=pages):
        assert rcsb_search_both(rows=2) == ["AAAA", "BBBB", "CCCC"]
    with um.patch.object(requests, "post", return_value=Resp(204)):
        assert rcsb_search_both() == []
    print("selftest OK")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["meta", "search", "mapping", "epitope", "diagnose", "interface", "apo", "pd1ensemble", "versions", "selftest"])
    ap.add_argument("--specs", help="JSON chain specification for the epitope step")
    ap.add_argument("--residue", type=int, default=78)
    ap.add_argument("--pdb", help="PDB id (diagnose: restrict to one entry; interface: the entry to compare with 4ZQK)")
    ap.add_argument("--ids", nargs="*", help="PDB ids for the meta step (default: the seven entries in the study)")
    ap.add_argument("--pd1-chain", help="PD-1 chain letter (interface step)")
    ap.add_argument("--pdl1-chain", help="PD-L1 chain letter (interface step)")
    a = ap.parse_args()
    if a.step == "meta": step_meta(a.ids)
    elif a.step == "search": step_search()
    elif a.step == "mapping": step_mapping(residue=a.residue)
    elif a.step == "epitope": step_epitope(a.specs)
    elif a.step == "diagnose": step_diagnose(a.specs, a.pdb)
    elif a.step == "interface": step_interface(a.pdb, a.pd1_chain, a.pdl1_chain)
    elif a.step == "apo": step_apo()
    elif a.step == "pd1ensemble": step_pd1ensemble()
    elif a.step == "versions": step_versions()
    else: _selftest()
