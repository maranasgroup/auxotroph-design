"""
double_knockout.py
==================
Double-gene knockout analysis — synthetic lethal pairs across 46
carbon sources.  Uses COBRApy's double_gene_deletion (vectorised +
parallel).  Only gene pairs where NEITHER gene is individually
essential on the target substrate are tested (true synthetic lethals).

Run order: single_knockout.py → double_knockout.py

Outputs
-------
  double_ko_conditionally_essential.csv -- CE pairs per substrate
  double_ko_rescue.csv                  -- rescue reactions per pair
"""
import sys, time, multiprocessing, pandas as pd

from knockout_utils import (load_model, load_off_reactions, load_carbon_sources,
                             deduplicate_rescues,
                             prepare_condition, get_essential_genes,
                             get_essential_pairs, load_reaction_library,
                             find_rescue_genes_from_library, OUT)

def main():
    print("Loading...")
    model      = load_model()
    off_rxns   = load_off_reactions()
    cs_map     = load_carbon_sources()
    lib_df, mi = load_reaction_library()
    external_model = None  # set to loaded model for Layer 2 rescue

    all_genes = [g.id for g in model.genes if g.reactions]

    # ── Glucose baseline ──────────────────────────────────────────────────────
    glc_off = off_rxns.get("Glucose", set())
    with model:
        prepare_condition(model, cs_map["Glucose"], glc_off)
        wt_glc  = model.slim_optimize()
        glc_ess = get_essential_genes(model, wt_glc, processes=4)
    screen = [g for g in all_genes[:20] if g not in glc_ess]
    print(f"Glucose WT={wt_glc:.4f}  non-essential candidates: {len(screen)}")

    # ── Per-substrate loop ────────────────────────────────────────────────────
    ce_rows, rescue_rows = [], []

    for cs_name, ex_id in cs_map.items():
        if cs_name == "Glucose":
            continue
        off = off_rxns.get(cs_name, set())
        t0  = time.time()

        with model:
            prepare_condition(model, ex_id, off)
            wt = model.slim_optimize()
            if wt < 1e-6:
                continue
            print("finding single essential")
            # Single essentials on this substrate → excluded from pairs
            s_ess      = get_essential_genes(model, wt, gene_list=screen, processes=4)
            pair_cands = [g for g in screen if g not in s_ess]
            print("Double pair deletion")
            # Double deletion via COBRApy built-in
            lethal_pairs = get_essential_pairs(model, wt,
                                               gene_list=pair_cands, processes=4)

        # CE2^l: pairs essential on substrate, neither gene essential on glucose
        ce_pairs = [p for p in lethal_pairs
                    if not any(g in glc_ess for g in p)]

        dt = time.time() - t0
        print(f"  {cs_name:<35} WT={wt:.4f}  CE_pairs={len(ce_pairs)}"
              f"  ({dt:.0f}s)")

        for pair in ce_pairs:
            g1, g2 = sorted(pair)
            ce_rows.append({
                "substrate"  : cs_name,
                "gene1"      : g1, "gene2": g2,
                "gene1_name" : model.genes.get_by_id(g1).name,
                "gene2_name" : model.genes.get_by_id(g2).name,
                "ko_type"    : "double",
                "wt_growth"  : round(wt, 4),
            })

            rescues = find_rescue_genes_from_library(
                model, ko_gene_ids=[g1, g2],
                carbon_exchange_id=ex_id, off_rxn_ids=off, wt_growth=wt,
                library_df=lib_df, multi_index=mi, external_model=external_model)

            for r in rescues:
                rescue_rows.append({"substrate": cs_name,
                                     "ko_gene1": g1, "ko_gene2": g2,
                                     "ko_type": "double", **r})

    # ── Save ──────────────────────────────────────────────────────────────────
    ce_df     = pd.DataFrame(ce_rows)
    rescue_df = pd.DataFrame(rescue_rows)

    ce_df.to_csv(OUT / "double_ko_conditionally_essential.csv", index=False)
    rescue_df.to_csv(OUT / "double_ko_rescue_raw.csv", index=False)
    rescue_dedup = deduplicate_rescues(rescue_df)
    rescue_dedup.to_csv(OUT / "double_ko_rescue_dedup.csv", index=False)

    print(f"\n  ✓ double_ko_conditionally_essential.csv  {len(ce_df)} CE pairs")
    print(f"  ✓ double_ko_rescue.csv                   {len(rescue_df)} rescues")


if __name__ == "__main__":
    multiprocessing.set_start_method("fork", force=True)
    main()