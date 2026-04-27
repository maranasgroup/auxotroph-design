"""
single_knockout.py
==================
Single-gene knockout analysis across 46 carbon sources.
Uses COBRApy's single_gene_deletion (vectorised + parallel) for
essentiality, and the multi-key reaction library for ResMut rescue.

Run order: build_reaction_library.py → single_knockout.py

Outputs
-------
  single_ko_es_matrix.csv               -- binary ES[substrate × gene]
  single_ko_conditionally_essential.csv -- CE^l per substrate
  single_ko_rescue.csv                  -- rescue reactions per CE gene
"""
import sys, time, multiprocessing, pandas as pd
from pathlib import Path
REPO_ROOT = Path(__file__).resolve().parent      # or .parent.parent if scripts live in src/
DATA      = REPO_ROOT / "data"
OUT       = REPO_ROOT / "outputs"
OUT.mkdir(exist_ok=True)
from src.knockout_utils import (load_model, load_off_reactions, load_carbon_sources,
                             prepare_condition, get_essential_genes,
                             load_reaction_library,
                             find_rescue_genes_from_library,
                             deduplicate_rescues, OUT)

def main():
    print("Loading model, data, and reaction library...")
    model      = load_model()
    off_rxns   = load_off_reactions()
    cs_map     = load_carbon_sources()
    lib_df, mi = load_reaction_library()

    external_model = None

    all_genes = [g.id for g in model.genes if g.reactions]
    print(f"  Library: {len(lib_df)} reactions | {len(mi)} index keys")

    # ── Glucose baseline ──────────────────────────────────────────────────────
    glc_off = off_rxns.get("Glucose", set())
    glc_ex  = cs_map["Glucose"]
    with model:
        prepare_condition(model, glc_ex, glc_off)
        wt_glc  = model.slim_optimize()
        glc_ess = get_essential_genes(model, wt_glc, processes=4)
    print(f"Glucose  WT={wt_glc:.4f}  essential={len(glc_ess)}")
    screen = [g for g in all_genes if g not in glc_ess]

    # ── Per-substrate loop ────────────────────────────────────────────────────
    es_rows, ce_rows, rescue_rows = [], [], []

    for cs_name, ex_id in cs_map.items():
        if cs_name == "Glucose":
            continue
        off = off_rxns.get(cs_name, set())
        t0  = time.time()

        with model:
            prepare_condition(model, ex_id, off)
            wt  = model.slim_optimize()
            if wt < 1e-6:
                print(f"  {cs_name}: no growth, skipping"); continue
            ess = get_essential_genes(model, wt, gene_list=screen, processes=4)

        ce = [g for g in ess if g not in glc_ess]
        print(f"  {cs_name:<35} WT={wt:.4f}  ess={len(ess)}  CE={len(ce)}"
              f"  ({time.time()-t0:.1f}s)")

        for gid in all_genes:
            es_rows.append({"substrate": cs_name, "gene": gid,
                            "essential": int(gid in ess)})

        for gid in ce:
            gname = model.genes.get_by_id(gid).name
            ce_rows.append({"substrate": cs_name, "gene_id": gid,
                            "gene_name": gname, "ko_type": "single",
                            "wt_growth": round(wt, 4)})

            rescues = find_rescue_genes_from_library(
                model, ko_gene_ids=[gid],
                carbon_exchange_id=ex_id, off_rxn_ids=off, wt_growth=wt,
                library_df=lib_df, multi_index=mi, external_model=external_model)

            for r in rescues:
                rescue_rows.append({"substrate": cs_name, "ko_gene": gid,
                                     "ko_gene_name": gname, "ko_type": "single",
                                     **r})

    # ── Save ──────────────────────────────────────────────────────────────────
    es_df    = pd.DataFrame(es_rows)
    es_pivot = (es_df.pivot_table(index="substrate", columns="gene",
                                   values="essential", fill_value=0)
                  .reset_index())
    es_pivot.columns.name = None
    ce_df     = pd.DataFrame(ce_rows)
    rescue_df = pd.DataFrame(rescue_rows)

    es_pivot.to_csv(OUT / "single_ko_es_matrix.csv", index=False)
    ce_df.to_csv(OUT / "single_ko_conditionally_essential.csv", index=False)
    rescue_df.to_csv(OUT / "single_ko_rescue_raw.csv", index=False)
    rescue_dedup = deduplicate_rescues(rescue_df)
    rescue_dedup.to_csv(OUT / "single_ko_rescue_dedup.csv", index = False)

    print(f"\n  ✓ single_ko_es_matrix.csv              {es_pivot.shape}")
    print(f"  ✓ single_ko_conditionally_essential.csv {len(ce_df)} CE genes")
    print(f"  ✓ single_ko_rescue.csv                  {len(rescue_df)} rescues")


if __name__ == "__main__":
    multiprocessing.set_start_method("fork", force=True)
    main()