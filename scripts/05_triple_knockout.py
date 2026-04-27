"""
triple_knockout.py
==================
Triple-gene knockout analysis — targeted synthetic lethal triples.
Run order: single_knockout.py → double_knockout.py → triple_knockout.py

Outputs
-------
  triple_ko_conditionally_essential.csv -- CE triples per substrate
  triple_ko_rescue.csv                  -- rescue reactions per triple (deduped)
  triple_ko_rescue_raw.csv              -- raw rescue table
"""

import sys, time, multiprocessing, pandas as pd
from itertools import combinations
from collections import defaultdict
from src.knockout_utils import (load_model, load_off_reactions, load_carbon_sources,
                             prepare_condition, load_reaction_library,
                             find_rescue_genes_from_library,
                             deduplicate_rescues, ESS_THRESH, OUT)


def load_precomputed_essentials(cs_names: list, glc_ess: set):
    """
    Load single and double KO essentials from previously saved CSV files.

    Returns
    -------
    s_ess_by_cs : {cs_name: set of gene_id strings}
        Genes essential on each substrate (and NOT essential on glucose).
        These are the CE^l single genes.

    d_ess_by_cs : {cs_name: set of frozensets of 2 gene_ids}
        Synthetic lethal pairs on each substrate (neither gene in the pair
        is essential alone on that substrate or on glucose).
    """
    # Single CE genes per substrate
    s_ess_by_cs: dict[str, set] = {}
    try:
        se_df = pd.read_csv(OUT / "single_ko_conditionally_essential.csv")
        for cs_name in cs_names:
            rows = se_df[se_df.substrate == cs_name]
            # CE genes: essential on cs but not on glucose (already filtered in single_knockout.py)
            s_ess_by_cs[cs_name] = set(rows["gene_id"].astype(str))
        print(f"  Loaded single CE: "
              f"{sum(len(v) for v in s_ess_by_cs.values())} gene-substrate pairs")
    except FileNotFoundError:
        print("  WARNING: single_ko_conditionally_essential.csv not found. "
              "Run single_knockout.py first.")
        s_ess_by_cs = {cs: set() for cs in cs_names}

    # Double CE pairs per substrate
    d_ess_by_cs: dict[str, set] = {}
    try:
        de_df = pd.read_csv(OUT / "double_ko_conditionally_essential.csv")
        for cs_name in cs_names:
            rows = de_df[de_df.substrate == cs_name]
            pairs = {frozenset([r.gene1, r.gene2]) for _, r in rows.iterrows()}
            d_ess_by_cs[cs_name] = pairs
        print(f"  Loaded double CE: "
              f"{sum(len(v) for v in d_ess_by_cs.values())} pair-substrate pairs")
    except FileNotFoundError:
        print("  WARNING: double_ko_conditionally_essential.csv not found. "
              "Run double_knockout.py first.")
        d_ess_by_cs = {cs: set() for cs in cs_names}

    return s_ess_by_cs, d_ess_by_cs


def main():
    print("Loading model, data, and reaction library...")
    model      = load_model()
    off_rxns   = load_off_reactions()
    cs_map     = load_carbon_sources()
    lib_df, mi = load_reaction_library()
    external_model = None
    cs_names   = list(cs_map.keys())
    all_genes  = [g.id for g in model.genes if g.reactions]

    # ── Glucose baseline (fast — just to know which genes are glc-essential) ──
    # We only need the glc_ess set to define the screen space.
    # The expensive per-substrate single/double essentials come from saved CSVs.
    print("Computing glucose essentials baseline...")
    glc_off = off_rxns.get("Glucose", set())
    with model:
        prepare_condition(model, cs_map["Glucose"], glc_off)
        wt_glc  = model.slim_optimize()
        # Load from single KO file if available to avoid recomputation
        try:
            se_df = pd.read_csv(OUT / "single_ko_es_matrix.csv")
            # Genes essential on glucose: in the ES matrix, glucose row = 1
            # Actually glucose was the baseline, not in the matrix.
            # Use the fact that screen genes = all model genes minus glc_ess.
            # Load glc_ess from single_ko es_matrix if it has a glucose column
            if "Glucose" in se_df.columns:
                glc_ess = set(se_df.columns[se_df.loc[se_df.substrate=="Glucose"].iloc[0] == 1].tolist()) \
                    if "substrate" in se_df.columns else set()
            # Fallback: recompute (fast, ~30s)
            from knockout_utils import get_essential_genes
            glc_ess = get_essential_genes(model, wt_glc, processes=4)
        except Exception:
            from knockout_utils import get_essential_genes
            glc_ess = get_essential_genes(model, wt_glc, processes=4)

    screen     = set(all_genes) - glc_ess
    screen_lst = sorted(screen)
    print(f"  Glucose-essential: {len(glc_ess)} | Screen: {len(screen_lst)}")

    # ── Load pre-computed single + double essentials ───────────────────────────
    print("\nLoading pre-computed single and double KO results...")
    s_ess_by_cs, d_ess_by_cs = load_precomputed_essentials(cs_names, glc_ess)

    # ── Build reaction-neighbourhood candidate triples ─────────────────────────
    print("\nBuilding targeted triple candidates via reaction neighbourhood...")
    gene_rxns = defaultdict(set)
    rxn_genes = defaultdict(set)
    for rxn in model.reactions:
        for g in rxn.genes:
            gene_rxns[g.id].add(rxn.id)
            rxn_genes[rxn.id].add(g.id)

    candidates = set()

    # Isozyme groups — genes sharing the same reaction (≥3 genes)
    for rid, genes in rxn_genes.items():
        g_list = sorted(genes & screen)
        if len(g_list) >= 3:
            for triple in combinations(g_list, 3):
                candidates.add(triple)

    # Pairs sharing a reaction + one reaction-neighbouring gene
    for rid, genes in rxn_genes.items():
        g_list = sorted(genes & screen)
        if len(g_list) == 2:
            g1, g2 = g_list
            nbr_rxns = (gene_rxns[g1] | gene_rxns[g2]) - {rid}
            third = set()
            for r2 in nbr_rxns:
                third |= rxn_genes[r2] & screen
            for g3 in third - {g1, g2}:
                candidates.add(tuple(sorted([g1, g2, g3])))

    candidates = list(candidates)
    print(f"  Candidate triples to test: {len(candidates):,}")

    # ── Per-substrate triple KO loop ───────────────────────────────────────────
    print("\nRunning triple-gene knockout screening...")
    ce_rows, rescue_rows = [], []
    t_total = time.time()

    for cs_name, ex_id in cs_map.items():
        if cs_name == "Glucose":
            continue
        off = off_rxns.get(cs_name, set())
        t0  = time.time()

        # WT biomass for this substrate
        with model:
            prepare_condition(model, ex_id, off)
            wt = model.slim_optimize()
        if wt < 1e-6:
            continue

        threshold = wt * ESS_THRESH

        # Load pre-computed essentials for this substrate
        s_ess  = s_ess_by_cs.get(cs_name, set())   # genes essential alone
        d_ess  = d_ess_by_cs.get(cs_name, set())   # synthetic lethal pairs

        found = 0
        with model:
            prepare_condition(model, ex_id, off)

            for (g1, g2, g3) in candidates:
                # ── Filter 1: skip if any single gene is already essential ────
                if any(g in s_ess or g in glc_ess for g in [g1, g2, g3]):
                    continue
                # ── Filter 2: skip if any pair is already a lethal double ─────
                if any(frozenset(p) in d_ess
                       for p in [(g1, g2), (g1, g3), (g2, g3)]):
                    continue

                # ── Test triple knockout ──────────────────────────────────────
                with model:
                    for gid in [g1, g2, g3]:
                        try: model.genes.get_by_id(gid).knock_out()
                        except KeyError: pass
                    sol = model.optimize()
                    bm  = sol.objective_value if sol.status == "optimal" else 0.0

                if bm < threshold:
                    found += 1
                    ce_rows.append({
                        "substrate"  : cs_name,
                        "gene1"      : g1, "gene2": g2, "gene3": g3,
                        "gene1_name" : model.genes.get_by_id(g1).name,
                        "gene2_name" : model.genes.get_by_id(g2).name,
                        "gene3_name" : model.genes.get_by_id(g3).name,
                        "ko_type"    : "triple",
                        "wt_growth"  : round(wt, 4),
                    })

        dt = time.time() - t0
        print(f"  {cs_name:<35} WT={wt:.4f}  CE_triples={found}  ({dt:.0f}s)")

    print(f"\n  Total triple screening: {time.time()-t_total:.0f}s")

    # ── ResMut rescue for all triples ─────────────────────────────────────────
    print("\nRunning rescue for triple KO designs...")
    for ce in ce_rows:
        off   = off_rxns.get(ce["substrate"], set())
        ex_id = cs_map[ce["substrate"]]
        rescues = find_rescue_genes_from_library(
            model,
            ko_gene_ids=[ce["gene1"], ce["gene2"], ce["gene3"]],
            carbon_exchange_id=ex_id, off_rxn_ids=off,
            wt_growth=ce["wt_growth"],
            library_df=lib_df, multi_index=mi, external_model=external_model)

        for r in rescues:
            rescue_rows.append({
                "substrate": ce["substrate"],
                "ko_gene1" : ce["gene1"],
                "ko_gene2" : ce["gene2"],
                "ko_gene3" : ce["gene3"],
                "ko_type"  : "triple", **r})

    # ── Save ──────────────────────────────────────────────────────────────────
    ce_df     = pd.DataFrame(ce_rows)
    rescue_df = pd.DataFrame(rescue_rows)
    rescue_dedup = deduplicate_rescues(rescue_df) if len(rescue_df) else rescue_df

    ce_df.to_csv(OUT / "triple_ko_conditionally_essential.csv", index=False)
    rescue_df.to_csv(OUT / "triple_ko_rescue_raw.csv", index=False)
    rescue_dedup.to_csv(OUT / "triple_ko_rescue.csv", index=False)

    print(f"\n  ✓ triple_ko_conditionally_essential.csv  {len(ce_df)} CE triples")
    print(f"  ✓ triple_ko_rescue_raw.csv               {len(rescue_df)} rows")
    print(f"  ✓ triple_ko_rescue.csv                   {len(rescue_dedup)} unique rescue functions")


if __name__ == "__main__":
    multiprocessing.set_start_method("fork", force=True)
    main()