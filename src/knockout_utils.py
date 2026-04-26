"""
knockout_utils.py
=================
Shared utilities for single / double / triple knockout analyses.
Uses COBRApy built-in deletion functions (vectorised + parallel) for
essentiality screening, and implements the ResMut rescue mechanism
using a multi-key reaction library (EC + BiGG + MetaNetX).
"""
import cobra
import pandas as pd
import warnings
from collections import defaultdict
from pathlib import Path

from cobra.flux_analysis import single_gene_deletion, double_gene_deletion

warnings.filterwarnings("ignore")

BASE = Path("../imL1515")
OUT = Path("")

UPTAKE_RATE = 10.0
O2_RATE = 20.0
ATPM_MAINT = 8.39
RESCUE_THRESH = 0.10
ESS_THRESH = 0.10

ALWAYS_OPEN = {
    "EX_pi_e", "EX_nh4_e", "EX_so4_e", "EX_mg2_e", "EX_fe2_e", "EX_fe3_e",
    "EX_h2o_e", "EX_h_e", "EX_co2_e", "EX_k_e", "EX_na1_e", "EX_cl_e",
    "EX_mobd_e", "EX_cobalt2_e", "EX_zn2_e", "EX_mn2_e", "EX_ni2_e",
    "EX_cu2_e", "EX_ca2_e", "EX_sel_e", "EX_slnt_e", "EX_tungs_e",
}


# ── Loaders ───────────────────────────────────────────────────────────────────
def load_model():
    return cobra.io.read_sbml_model(str(BASE / "iML1515.xml"))


def load_off_reactions():
    df = pd.read_csv(OUT / "multi_carbon_off_reactions.csv")
    return (df[df.flag == "applied"]
            .groupby("carbon_source")["reaction_id"].apply(set).to_dict())


def load_carbon_sources():
    df = pd.read_csv(OUT / "multi_carbon_growth_summary.csv")
    return dict(zip(df.carbon_source, df.exchange_id))


# ── Medium + regulatory setup ─────────────────────────────────────────────────
def setup_medium(model, carbon_exchange_id, uptake=UPTAKE_RATE, o2=O2_RATE):
    for r in model.reactions:
        if r.id.startswith("EX_") and r.lower_bound < 0 and r.id not in ALWAYS_OPEN:
            r.lower_bound = 0.0
    model.reactions.get_by_id("EX_o2_e").lower_bound = -o2
    model.reactions.get_by_id(carbon_exchange_id).lower_bound = -uptake
    model.reactions.get_by_id("ATPM").lower_bound = ATPM_MAINT


def apply_regulatory_off(model, off_rxn_ids):
    for rid in off_rxn_ids:
        try:
            r = model.reactions.get_by_id(rid)
            r.lower_bound = r.upper_bound = 0.0
        except KeyError:
            pass


def prepare_condition(model, carbon_exchange_id, off_rxn_ids):
    setup_medium(model, carbon_exchange_id)
    apply_regulatory_off(model, off_rxn_ids)


# ── COBRApy built-in essentiality ─────────────────────────────────────────────
def get_essential_genes(model, wt_growth, gene_list=None, processes=4):
    threshold = wt_growth * ESS_THRESH
    result = single_gene_deletion(model, gene_list=gene_list,
                                  method="fba", processes=processes)
    essential = result.loc[
        result["growth"].isna() | (result["growth"] < threshold), "ids"
    ]
    return {g for ids in essential for g in ids}


def get_essential_pairs(model, wt_growth, gene_list, processes=4):
    threshold = wt_growth * ESS_THRESH
    result = double_gene_deletion(model,
                                  gene_list1=gene_list, gene_list2=gene_list,
                                  method="fba", processes=processes)
    lethal = result.loc[
        result["growth"].isna() | (result["growth"] < threshold), "ids"
    ]
    return {frozenset(ids) for ids in lethal}


# ── Reaction library ──────────────────────────────────────────────────────────
def load_reaction_library(path=None):
    p = path or OUT / "reaction_library.csv"
    lib = pd.read_csv(p)
    idx = defaultdict(list)
    for _, row in lib.iterrows():
        for ec in str(row["ec_numbers"]).split(";"):
            ec = ec.strip()
            if ec and ec != "nan": idx[("ec", ec)].append(row["reaction_id"])
        bigg = str(row["bigg_id"]).strip()
        if bigg and bigg != "nan": idx[("bigg", bigg)].append(row["reaction_id"])
        mnx = str(row["metanetx_id"]).strip()
        if mnx and mnx != "nan":  idx[("mnx", mnx)].append(row["reaction_id"])
    return lib, dict(idx)


def _get_rxn_keys(rxn):
    """Extract (ecs, bigg, mnx) from a COBRApy reaction."""
    import re
    EC_PATTERN = re.compile(r'^\d+\.\d+\.\d+[\.\d\-]*$')
    ann = rxn.annotation
    notes = getattr(rxn, "notes", {}) or {}

    def _collect(aliases):
        vals = []
        for a in aliases:
            v = ann.get(a, [])
            if isinstance(v, str):
                v = [v]
            elif not isinstance(v, list):
                v = [str(v)]
            vals.extend(v)
        return [str(v).strip() for v in vals if str(v).strip() not in ("", "nan")]

    ecs = sorted({v for v in _collect({"ec-code", "ec", "EC"}) if EC_PATTERN.match(v)})
    bigg = (_collect({"bigg.reaction", "bigg", "BIGG"}) or [""])[0]
    if not bigg:
        orig = notes.get("original_bigg_ids", [])
        bigg = (orig[0] if isinstance(orig, list) and orig else str(orig or "")).strip()
    if not bigg and not rxn.id.startswith(("MNXR", "R0")):
        bigg = rxn.id
    mnx_raw = _collect({"metanetx.reaction", "metanetx", "mnx", "MNX"})
    mnx = mnx_raw[0] if mnx_raw else (rxn.id if rxn.id.startswith("MNXR") else "")
    return ecs, bigg, mnx


# ── ResMut rescue (CORRECTED) ─────────────────────────────────────────────────
def find_rescue_genes_from_library(
        model,
        ko_gene_ids: list,
        carbon_exchange_id: str,
        off_rxn_ids: set,
        wt_growth: float,
        library_df: pd.DataFrame,
        multi_index: dict,
        external_model=None,
) -> list[dict]:
    """
    Correct ResMut implementation matching the paper's formulation:

      For each r in D^k (reactions inactivated by KO):
        Force r open (set lb=0 / ub=1000) — this simulates adding a gene
        from the library that catalyzes reaction r.
        If biomass >= 10% WT → r rescues growth.
        Report all library entries sharing EC/BiGG/MNX with r as the
        set of rescue gene candidates (from any organism).

    The PREVIOUS (wrong) approach: unblock a DIFFERENT reaction sharing
    an EC with the blocked one. This was wrong because:
      - The model may already have that alternative reaction open (no effect)
      - For AND-complex KOs, only the exact same reaction can rescue
      - The paper's rescue gene adds the missing reaction, not an isozyme

    Layer 2 (external model): reactions not in iML1515 are temporarily
    added to the model to test if they rescue growth.
    """
    rescue_threshold = wt_growth * RESCUE_THRESH
    host_ids = {r.id for r in model.reactions}

    # ── Step 1: Identify D^k ─────────────────────────────────────────────────
    with model:
        prepare_condition(model, carbon_exchange_id, off_rxn_ids)
        bounds_pre = {r.id: (r.lower_bound, r.upper_bound) for r in model.reactions}
        for gid in ko_gene_ids:
            try:
                model.genes.get_by_id(gid).knock_out()
            except KeyError:
                pass
        dk = [
            r for r in model.reactions
            if r.gene_reaction_rule
               and (bounds_pre[r.id][0] != 0 or bounds_pre[r.id][1] != 0)
               and r.lower_bound == 0 and r.upper_bound == 0
        ]

    if not dk:
        return []

    results = []

    with model:
        prepare_condition(model, carbon_exchange_id, off_rxn_ids)
        for gid in ko_gene_ids:
            try:
                model.genes.get_by_id(gid).knock_out()
            except KeyError:
                pass

        for ko_rxn in dk:
            orig_lb, orig_ub = bounds_pre[ko_rxn.id]

            # ── Step 2: Force D^k reaction open ──────────────────────────────
            with model:
                # Set ub first to avoid lb > ub when both are 0
                ko_rxn.upper_bound = orig_ub if orig_ub > 0 else 1000.0
                ko_rxn.lower_bound = orig_lb if orig_lb < 0 else 0.0

                sol = model.optimize()
                bm = sol.objective_value if sol.status == "optimal" else 0.0

                if bm < rescue_threshold:
                    continue  # reaction doesn't rescue even if restored

                # ── Step 3: Find library analogues ───────────────────────────
                ecs, bigg, mnx = _get_rxn_keys(ko_rxn)

                lib_candidates = []  # (library_reaction_id, match_key, match_type)
                seen = set()

                # Layer 1: internal iML1515 reactions sharing a key
                for ec in ecs:
                    for lid in multi_index.get(("ec", ec), []):
                        if lid != ko_rxn.id and lid not in seen:
                            lib_candidates.append((lid, f"ec:{ec}", "ec"))
                            seen.add(lid)
                if bigg:
                    for lid in multi_index.get(("bigg", bigg), []):
                        if lid != ko_rxn.id and lid not in seen:
                            lib_candidates.append((lid, f"bigg:{bigg}", "bigg"))
                            seen.add(lid)
                if mnx:
                    for lid in multi_index.get(("mnx", mnx), []):
                        if lid != ko_rxn.id and lid not in seen:
                            lib_candidates.append((lid, f"mnx:{mnx}", "mnx"))
                            seen.add(lid)

                # Layer 2: external model reactions sharing a key
                if external_model is not None or (library_df.layer == 2).any():
                    layer2 = library_df[library_df.layer == 2]
                    for _, lr in layer2.iterrows():
                        lr_ecs = [e.strip() for e in str(lr.ec_numbers).split(";") if e.strip() and e.strip() != "nan"]
                        lr_bigg = str(lr.bigg_id).strip()
                        lr_mnx = str(lr.metanetx_id).strip()
                        if any(ec in lr_ecs for ec in ecs) \
                                or (bigg and bigg == lr_bigg) \
                                or (mnx and mnx == lr_mnx and lr_mnx not in ("", "nan")):
                            if lr.reaction_id not in seen:
                                mtype = ("ec" if any(ec in lr_ecs for ec in ecs)
                                         else "bigg" if bigg == lr_bigg else "mnx")
                                mk = (f"ec:{ecs[0]}" if mtype == "ec"
                                      else f"bigg:{bigg}" if mtype == "bigg"
                                else f"mnx:{mnx}")
                                lib_candidates.append((lr.reaction_id, mk, mtype))
                                seen.add(lr.reaction_id)

                # If no library match, record as direct-restore
                if not lib_candidates:
                    results.append({
                        "ko_reaction": ko_rxn.id,
                        "rescue_reaction": ko_rxn.id,
                        "rescue_rxn_name": ko_rxn.name,
                        "match_key": "direct_restore",
                        "match_type": "direct_restore",
                        "rescue_source": "no_library_analogue",
                        "rescue_layer": 0,
                        "rescue_gene_rule": ko_rxn.gene_reaction_rule,
                        "rescued_biomass": round(bm, 5),
                        "wt_biomass": round(wt_growth, 5),
                        "pct_rescue": round(100 * bm / wt_growth, 1),
                    })
                    continue

                for (lib_rxn_id, match_key, match_type) in lib_candidates:
                    lib_rows = library_df[library_df.reaction_id == lib_rxn_id]
                    if len(lib_rows) == 0:
                        continue
                    lr = lib_rows.iloc[0]
                    results.append({
                        "ko_reaction": ko_rxn.id,
                        "rescue_reaction": lib_rxn_id,
                        "rescue_rxn_name": lr.reaction_name,
                        "match_key": match_key,
                        "match_type": match_type,
                        "rescue_source": lr.source,
                        "rescue_layer": int(lr.layer),
                        "rescue_gene_rule": lr.gene_rule,
                        "rescued_biomass": round(bm, 5),
                        "wt_biomass": round(wt_growth, 5),
                        "pct_rescue": round(100 * bm / wt_growth, 1),
                    })

    return results


# ── Rescue deduplication ──────────────────────────────────────────────────────
def deduplicate_rescues(rescue_df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse the raw rescue table to one row per (substrate, ko_gene, rescue_function).

    Raw rescue has one row per (ko_reaction × rescue_reaction) combination.
    For a gene controlling N reactions all with the same EC, and M library
    entries sharing that EC, this produces N×M rows for one biological fact.

    The paper's unit of interest is:
      'Gene k on substrate l requires function F to restore growth'

    Deduplication key: (substrate, ko_gene, match_key)
      match_key encodes the rescue FUNCTION (e.g. 'ec:2.4.1.25', 'bigg:ATPS4rpp')

    For each unique (substrate, ko_gene, match_key) group:
      - Keep the row with the highest pct_rescue
      - Collect all unique rescue_reactions and rescue_sources as lists
    """
    if len(rescue_df) == 0:
        return rescue_df

    ko_col = "ko_gene" if "ko_gene" in rescue_df.columns else \
        "ko_gene1" if "ko_gene1" in rescue_df.columns else None
    sub_col = "substrate"

    if ko_col is None:
        return rescue_df  # double/triple KO — different structure, skip

    # Group by the biological rescue unit
    group_cols = [sub_col, ko_col, "match_key", "match_type"]
    group_cols = [c for c in group_cols if c in rescue_df.columns]

    records = []
    for keys, grp in rescue_df.groupby(group_cols):
        best = grp.loc[grp["pct_rescue"].idxmax()]
        unique_rescues = "; ".join(sorted(grp["rescue_reaction"].unique()))
        unique_sources = "; ".join(sorted(grp["rescue_source"].unique()))
        unique_ko_rxns = "; ".join(sorted(grp["ko_reaction"].unique()))
        unique_names = "; ".join(dict.fromkeys(grp["rescue_rxn_name"]))

        row = best.to_dict()
        row["rescue_reaction"] = unique_rescues
        row["rescue_rxn_name"] = unique_names
        row["rescue_source"] = unique_sources
        row["ko_reaction"] = unique_ko_rxns
        row["n_rescue_variants"] = len(grp["rescue_reaction"].unique())
        row["n_ko_reactions"] = len(grp["ko_reaction"].unique())
        records.append(row)

    dedup = pd.DataFrame(records)
    # Sort by substrate, ko_gene, then pct_rescue descending
    sort_cols = [c for c in [sub_col, ko_col, "pct_rescue"] if c in dedup.columns]
    return dedup.sort_values(sort_cols, ascending=[True, True, False]).reset_index(drop=True)
