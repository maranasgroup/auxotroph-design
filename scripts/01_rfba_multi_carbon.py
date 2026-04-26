"""
rFBA Multi-Carbon-Source Analysis
==================================
Runs the full regulatory constraint analysis for all 45 aerobic
'growth' substrates identified from the Biolog phenotyping data.

For each carbon source:
  1. WT FBA  → get growth rate + flux solution
  2. Evaluate all Boolean rules using WT fluxes
  3. GPR → candidate OFF reactions
  4. Feasibility-preserving constraint application
  5. rFBA → constrained growth rate

Outputs
-------
  multi_carbon_gene_activity.csv      — gene ON/OFF per gene × substrate
  multi_carbon_off_reactions.csv      — all applied OFF reactions
  multi_carbon_growth_summary.csv     — WT vs rFBA growth per substrate
  multi_carbon_reaction_heatmap.csv   — reaction × substrate OFF matrix
"""

import re, warnings, pandas as pd
import cobra
from pathlib import Path

warnings.filterwarnings("ignore")
BASE = Path("/Users/rbm5893/Documents/1.1 Roghaye Obisdian/1.1Research/GAP/publication/revision1/imL1515")
OUT  = Path("/Users/rbm5893/Documents/1.1 Roghaye Obisdian/1.1Research/GAP/publication/revision1/analysis")

# ── Load model & rules ────────────────────────────────────────────────────────
print("Loading model and rules...")
model = cobra.io.read_sbml_model(str(BASE / "iML1515.xml"))
imc   = pd.read_csv(BASE / "iMC1010.csv")
imc_rules = {str(r["Alias_3"]): str(r["Rule"])
             for _, r in imc.iterrows()
             if str(r["Alias_3"]).startswith("b") and str(r["Rule"]) not in ("","nan")}
norm_v2 = pd.read_csv(OUT / "delta_gene_rules_normalized_v2.csv")
delta_rules = {str(r["bnumber"]): str(r["normalized_rule_v2"])
               for _, r in norm_v2.iterrows()
               if str(r["bnumber"]).startswith("b") and str(r["normalized_rule_v2"]) not in ("","nan")}
all_rules = {**imc_rules, **delta_rules}
model_gene_ids = {g.id for g in model.genes}
print(f"  {len(model.genes)} genes | {len(model.reactions)} reactions | "
      f"{len(all_rules)} regulatory rules")

# ── Carbon source definitions ─────────────────────────────────────────────────
# All 45 aerobic 'growth' substrates from the Biolog data
# (glucose already analysed separately; included here for completeness)
CARBON_SOURCES = {
    # Sugars & sugar derivatives
    "Glucose"               : "EX_glc__D_e",
    "D-Fructose"            : "EX_fru_e",
    "D-Mannose"             : "EX_man_e",
    "D-Galactose"           : "EX_gal_e",
    "D-Ribose"              : "EX_rib__D_e",
    "D-Xylose"              : "EX_xyl__D_e",
    "L-Arabinose"           : "EX_arab__L_e",
    "beta-D-Allose"         : "EX_all__D_e",
    "L-Fucose"              : "EX_fuc__L_e",
    "L-Rhamnose"            : "EX_rmn_e",
    "D-Sorbitol"            : "EX_sbt__D_e",
    "D-Mannitol"            : "EX_mnl_e",
    # Sugar phosphates
    "Glucose-1-Phosphate"   : "EX_g1p_e",
    "Glucose-6-Phosphate"   : "EX_g6p_e",
    "Fructose-6-Phosphate"  : "EX_f6p_e",
    # Disaccharides
    "Maltose"               : "EX_malt_e",
    "Maltotriose"           : "EX_malttr_e",
    "Lactose"               : "EX_lcts_e",
    "D-Melibiose"           : "EX_melib_e",
    "D-Trehalose"           : "EX_tre_e",
    # Organic acids (TCA-related)
    "Succinic Acid"         : "EX_succ_e",
    "Fumaric Acid"          : "EX_fum_e",
    "L-Malic Acid"          : "EX_mal__L_e",
    "D-Malic Acid"          : "EX_mal__D_e",
    "Pyruvic Acid"          : "EX_pyr_e",
    "alpha-Keto-Glutaric"   : "EX_akg_e",
    # Organic acids (other)
    "L-Lactic Acid"         : "EX_lac__L_e",
    "D-Gluconic Acid"       : "EX_glcn_e",
    "5-Keto-D-Gluconate"    : "EX_5dglcn_e",
    "D-Glucuronic Acid"     : "EX_glcur_e",
    "D-Saccharic Acid"      : "EX_glcr_e",
    "Mucic Acid"            : "EX_galct__D_e",
    "D-Galactonic Acid"     : "EX_galctn__D_e",
    "L-Galactonic Acid"     : "EX_galctn__L_e",
    # Amino acids
    "D-Alanine"             : "EX_ala__D_e",
    "L-Serine"              : "EX_ser__L_e",
    "D-Serine"              : "EX_ser__D_e",
    "L-Aspartic Acid"       : "EX_asp__L_e",
    # Nucleosides
    "Adenosine"             : "EX_adn_e",
    "Inosine"               : "EX_ins_e",
    "2-DeoxyAdenosine"      : "EX_dad_2_e",
    "Thymidine"             : "EX_thymd_e",
    "Uridine"               : "EX_uri_e",
    # Modified sugars
    "N-Acetyl-D-Glucosamine": "EX_acgam_e",
    "N-Acetyl-Neuraminic"   : "EX_acnam_e",
    # Polyols
    "Glycerol"              : "EX_glyc_e",
}

# ── Helper: minimal M9 medium ─────────────────────────────────────────────────
ALWAYS_OPEN = {
    "EX_pi_e","EX_nh4_e","EX_so4_e","EX_mg2_e","EX_fe2_e","EX_fe3_e",
    "EX_h2o_e","EX_h_e","EX_co2_e","EX_k_e","EX_na1_e","EX_cl_e",
    "EX_mobd_e","EX_cobalt2_e","EX_zn2_e","EX_mn2_e","EX_ni2_e",
    "EX_cu2_e","EX_ca2_e","EX_sel_e","EX_slnt_e","EX_tungs_e",
}

def set_medium(model_obj, carbon_exchange_id, uptake_rate=-10.0):
    for r in model_obj.reactions:
        if r.id.startswith("EX_") and r.lower_bound < 0 and r.id not in ALWAYS_OPEN:
            r.lower_bound = 0.0
    model_obj.reactions.get_by_id("EX_o2_e").lower_bound = -15.0
    model_obj.reactions.get_by_id(carbon_exchange_id).lower_bound = uptake_rate

# ── Signal dictionary builder ─────────────────────────────────────────────────
BASE_SIGNALS = {
    "o2_e": 15.0, "fe2_e": 1.0, "nh4_e": 10.0, "so4_e": 10.0,
    "pi_e": 10.0, "mg2_e": 1.0, "mobd_e": 1.0,
    "cu2_e": 0.001, "zn2_e": 0.001, "ni2_e": 0.001, "mn2_e": 0.001,
    "h2o2_e": 0.0, "no_e": 0.0, "no3_e": 0.0, "no2_e": 0.0,
    "ac_e": 0.0, "lac__D_e": 0.0, "lac__L_e": 0.0,
    "leu_DASH_L_e": 0.0, "lys_DASH_L_e": 0.0,
    "acgam_e": 0.0, "chb_e": 0.0, "phenac_e": 0.0,
    "glc__D_e": 0.0, "glc_DASH_D_e": 0.0,   # default: no glucose
    "BiomassEcoli": 1.0, "Stress": False, "Stringent": False,
    "Oxidative-Stress": False, "pH": 7.4, "Rich-Medium": False,
}

def make_signals(carbon_exchange_id, conc=10.0):
    sig = dict(BASE_SIGNALS)
    # Map exchange ID to signal key(s)
    met_key = carbon_exchange_id.replace("EX_","").replace("_e","_e")
    # Some metabolites have _DASH_ notation in rules
    sig[met_key] = conc
    # Also add common aliases
    aliases = {
        "EX_glc__D_e"  : ["glc__D_e","glc_DASH_D_e"],
        "EX_lac__L_e"  : ["lac__L_e","lac_DASH_L_e"],
        "EX_lac__D_e"  : ["lac__D_e","lac_DASH_D_e"],
        "EX_arab__L_e" : ["arab__L_e","larab_e"],
        "EX_xyl__D_e"  : ["xyl__D_e","xyl_DASH_D_e"],
        "EX_leu__L_e"  : ["leu__L_e","leu_DASH_L_e"],
        "EX_lys__L_e"  : ["lys__L_e","lys_DASH_L_e"],
        "EX_thymd_e"   : ["thymd_e","thy_e"],
        "EX_mal__L_e"  : ["mal__L_e","mal_DASH_L_e"],
        "EX_mal__D_e"  : ["mal__D_e","mal_DASH_D_e"],
        "EX_fuc__L_e"  : ["fuc__L_e","fuc_DASH_L_e"],
        "EX_ser__L_e"  : ["ser__L_e","ser_DASH_L_e"],
        "EX_ser__D_e"  : ["ser__D_e","ser_DASH_D_e"],
        "EX_asp__L_e"  : ["asp__L_e","asp_DASH_L_e"],
        "EX_ala__D_e"  : ["ala__D_e","ala_DASH_D_e"],
        "EX_lcts_e"    : ["lcts_e","lac__D_e","lac_DASH_D_e"],
        "EX_4abut_e"   : ["4abut_e","gam_e"],
        "EX_galctn__D_e": ["galctn__D_e","galctn_DASH_D_e"],
        "EX_all__D_e"  : ["all__D_e","all_DASH_D_e"],
    }
    for alt_key in aliases.get(carbon_exchange_id, []):
        sig[alt_key] = conc
    return sig

# ── Rule evaluator ────────────────────────────────────────────────────────────
def eval_rule(rule, cond, tf_act):
    if rule in ("TRUE","ON"):   return True
    if rule in ("FALSE","OFF"): return False
    w = rule
    def rmet(m):
        return "True" if float(cond.get(m.group(1),0)) > float(m.group(2)) else "False"
    w = re.sub(r'([\w_]+(?:_e|__[\w]+_e))\s*>\s*([\d.]+)', rmet, w)
    w = re.sub(r'BiomassEcoli\s*>\s*([\d.]+)',
        lambda m: "True" if float(cond.get("BiomassEcoli",0)) > float(m.group(1)) else "False", w)
    w = re.sub(r'pH\s*<\s*([\d.]+)',
        lambda m: "True" if cond.get("pH",7.4) < float(m.group(1)) else "False", w)
    def rflux(m):
        return "True" if float(cond.get(f"_flux_{m.group(1)}",0)) > float(m.group(2)) else "False"
    def rfluxlt(m):
        return "True" if float(cond.get(f"_flux_{m.group(1)}",0)) < float(m.group(2)) else "False"
    w = re.sub(r'\b([A-Z][A-Z0-9_]{1,})\s*>\s*([\d.]+)', rflux, w)
    w = re.sub(r'\b([A-Z][A-Z0-9_]{1,})\s*<\s*(-?[\d.]+)', rfluxlt, w)
    for bn in set(re.findall(r'\bb\d{4,5}\b', w)):
        w = re.sub(r'\b'+bn+r'\b', str(tf_act.get(bn,True)), w)
    for sig, val in cond.items():
        if sig.startswith("_flux_"): continue
        if isinstance(val, bool):
            w = re.sub(r'\b'+re.escape(sig)+r'\b', str(val), w)
        elif sig == "BiomassEcoli":
            w = re.sub(r'\bBiomassEcoli\b', str(val), w)
        elif sig == "pH":
            w = re.sub(r'\bpH\b', str(val), w)
    w = re.sub(r'\bNOT\b','not',w); w = re.sub(r'\bAND\b','and',w); w = re.sub(r'\bOR\b','or',w)
    w = w.replace("TRUE","True").replace("FALSE","False").replace("OFF","False")
    w = re.sub(r'\b(?!not\b|and\b|or\b|True\b|False\b)[A-Za-z][A-Za-z0-9_\-]*\b','True',w)
    try:    return bool(eval(w,{"__builtins__":{}}))
    except: return True

def precompute_tf_activity(cond, rules):
    act = {}
    for bn,rule in rules.items():
        if not re.search(r'\bb\d{4,5}\b', rule):
            act[bn] = eval_rule(rule, cond, act)
    for _ in range(5):
        changed = False
        for bn,rule in rules.items():
            if bn in act: continue
            deps = set(re.findall(r'\bb\d{4,5}\b', rule))
            if deps.issubset(act.keys()):
                act[bn] = eval_rule(rule, cond, act); changed = True
        if not changed: break
    for bn in rules: act.setdefault(bn, True)
    return act

def reaction_is_off(rxn, gene_act):
    gpr = rxn.gene_reaction_rule.strip()
    if not gpr: return False
    expr = gpr
    for gene in rxn.genes:
        expr = re.sub(r'\b'+re.escape(gene.id)+r'\b', str(gene_act.get(gene.id,True)), expr)
    expr = re.sub(r'\b(?!True\b|False\b)[A-Za-z][A-Za-z0-9_]*\b','True',expr)
    try:    return not bool(eval(expr,{"__builtins__":{}}))
    except: return False

def apply_rfba_safe(model_obj, candidate_ids, wt_fluxes):
    """Apply constraints in order of increasing |WT flux|; relax if infeasible."""
    sorted_ids = sorted(candidate_ids, key=lambda rid: abs(wt_fluxes.get(rid,0)))
    applied, relaxed = [], []
    for rid in sorted_ids:
        try: rxn = model_obj.reactions.get_by_id(rid)
        except: continue
        orig_lb, orig_ub = rxn.lower_bound, rxn.upper_bound
        rxn.lower_bound = rxn.upper_bound = 0.0
        sol = model_obj.optimize()
        if sol.status == "optimal" and sol.objective_value > 1e-6:
            applied.append(rid)
        else:
            rxn.lower_bound = orig_lb; rxn.upper_bound = orig_ub
            relaxed.append(rid)
    return applied, relaxed

# ── Subsystem assignment ──────────────────────────────────────────────────────
SUBSYS_MAP = [
    ('Transport',              ['tex','tpp','t2r','abc','antiport','symport','transport','permease','diffusion']),
    ('Iron acquisition',       ['iron','enterobactin','siderophore','fep','fhu']),
    ('Oxidative phosphorylation',['cytochrome','oxidase','atpase','atp synth','nadh dehydrog','electron']),
    ('Fatty acid metabolism',  ['fatty acid','acyl-coa','acoad','enoyl-coa','fadA','fadB']),
    ('Amino acid metabolism',  ['amino','aminotransfer','serine','glycine','threonine','lysine','methionine','cysteine','tryptophan','phenylalanine','tyrosine','proline','arginine','histidine','aspartase','tryptophanase']),
    ('Nucleotide metabolism',  ['nucleotide','purine','pyrimidine','adenylosuccin','phosphoribosyl','guanylate','uridine','thymidine','nucleoside','deoxyribose']),
    ('Cofactor biosynthesis',  ['cobalamin','thiamine','riboflavin','folate','biotin','heme','porphyrin','siroheme','ubiquinone','menaquinone','nad ','fad ','molybdo','coenzyme']),
    ('Sugar catabolism',       ['gluconate','galactonate','arabinose','xylose','fucose','rhamnose','glucuronate','galacturonate','fructose catab','sorbitol','mannitol']),
    ('Glyoxylate/Acetate',     ['glyoxylate','acetate','acetyl-coa synth','acs ','pta ','ackA']),
    ('TCA cycle',              ['tca','citrate synth','isocitrate','succinate','fumarase','malate dehydrogenase']),
    ('Cell envelope',          ['lps','lipopolysaccharide','lipid a','murein','peptidoglycan','outer membrane']),
    ('Inorganic ion',          ['sulfate','phosphate','nitrogen','molybdate','copper','zinc','nickel','manganese']),
    ('Oxidative stress',       ['glutathione','peroxidase','catalase','superoxide','thioredoxin']),
]
def get_subsystem(rid, rname):
    rname_l = rname.lower(); rid_l = rid.lower()
    for cat, kws in SUBSYS_MAP:
        if any(k in rname_l or k in rid_l for k in kws):
            return cat
    return 'Other'

rxn_name_map = {r.id: r.name for r in model.reactions}

# ── Main loop ─────────────────────────────────────────────────────────────────
print(f"\nRunning rFBA for {len(CARBON_SOURCES)} carbon sources...\n")
print(f"{'Carbon source':<28} {'WT':>7} {'rFBA':>7}  {'OFF':>4}  {'FlagB':>5}")
print("-" * 58)

growth_rows  = []
off_rxn_rows = []
gene_act_rows_all = {}  # {carbon_source: {gene: active}}

for cs_name, ex_id in CARBON_SOURCES.items():
    # Step 1: WT FBA
    with model:
        set_medium(model, ex_id)
        sol_wt = model.optimize()
    if sol_wt.status != "optimal" or sol_wt.objective_value < 1e-6:
        print(f"  {cs_name:<28} NO GROWTH")
        continue
    wt_growth  = sol_wt.objective_value
    wt_fluxes  = dict(sol_wt.fluxes)

    # Step 2: Build condition signal dict + add WT fluxes
    cond_sig = make_signals(ex_id, conc=10.0)
    for rxn in model.reactions:
        cond_sig[f"_flux_{rxn.id}"] = float(sol_wt.fluxes.get(rxn.id,0.0))
    cond_sig["_flux_GLCpts"] = abs(sol_wt.fluxes.get("EX_glc__D_e",0.0))

    # Step 3: Evaluate rules
    tf_act   = precompute_tf_activity(cond_sig, all_rules)
    gene_act = {gid: tf_act.get(gid,True) for gid in model_gene_ids}
    gene_act_rows_all[cs_name] = gene_act

    # Step 4: Candidate OFF reactions
    candidates = [rxn.id for rxn in model.reactions if reaction_is_off(rxn, gene_act)]

    # Step 5: Feasibility-preserving rFBA
    with model:
        set_medium(model, ex_id)
        applied, relaxed = apply_rfba_safe(model, candidates, wt_fluxes)
        sol_rfba = model.optimize()
    rfba_growth = sol_rfba.objective_value if sol_rfba.status=="optimal" else 0.0

    print(f"  {cs_name:<28} {wt_growth:>7.4f} {rfba_growth:>7.4f}  {len(applied):>4}  {len(relaxed):>5}")

    # Record OFF reactions
    for rid in applied + relaxed:
        rxn    = model.reactions.get_by_id(rid)
        rname  = rxn_name_map.get(rid,"")
        off_g  = [g.id for g in rxn.genes if not gene_act.get(g.id,True)]
        off_rxn_rows.append({
            "carbon_source" : cs_name,
            "exchange_id"   : ex_id,
            "reaction_id"   : rid,
            "reaction_name" : rname,
            "subsystem"     : get_subsystem(rid, rname),
            "gpr"           : rxn.gene_reaction_rule,
            "off_genes"     : "; ".join(off_g),
            "n_off_genes"   : len(off_g),
            "wt_flux"       : round(wt_fluxes.get(rid,0.0), 5),
            "flag"          : "applied" if rid in applied else "relaxed_FlagB",
        })

    growth_rows.append({
        "carbon_source"     : cs_name,
        "exchange_id"       : ex_id,
        "wt_growth"         : round(wt_growth, 4),
        "rfba_growth"       : round(rfba_growth, 4),
        "growth_impact"     : round(wt_growth - rfba_growth, 4),
        "pct_impact"        : round(100*(wt_growth-rfba_growth)/wt_growth, 2) if wt_growth>0 else 0,
        "n_applied"         : len(applied),
        "n_relaxed_FlagB"   : len(relaxed),
        "feasible"          : sol_rfba.status == "optimal",
    })

# ── Build output tables ───────────────────────────────────────────────────────
print("\n" + "="*65)
print("BUILDING SUMMARY TABLES")
print("="*65)

off_df    = pd.DataFrame(off_rxn_rows)
growth_df = pd.DataFrame(growth_rows).sort_values("wt_growth", ascending=False)

if len(off_df):
    applied_df = off_df[off_df.flag=="applied"]

    # Reaction × carbon-source heatmap (1=OFF, 0=ON)
    cs_list  = list(CARBON_SOURCES.keys())
    rxn_list = applied_df["reaction_id"].unique().tolist()

    heatmap = pd.DataFrame(0, index=rxn_list, columns=cs_list)
    for _, row in applied_df.iterrows():
        if row["reaction_id"] in heatmap.index and row["carbon_source"] in heatmap.columns:
            heatmap.loc[row["reaction_id"], row["carbon_source"]] = 1

    heatmap["n_cs_off"]      = heatmap.sum(axis=1)
    heatmap["reaction_name"] = heatmap.index.map(rxn_name_map)
    heatmap["subsystem"]     = [get_subsystem(rid, rxn_name_map.get(rid,""))
                                  for rid in heatmap.index]
    heatmap = heatmap.sort_values("n_cs_off", ascending=False)

    # Gene activity table (all carbon sources)
    gene_act_rows = []
    all_genes = sorted(model_gene_ids)
    for cs, gact in gene_act_rows_all.items():
        for gid, active in gact.items():
            gene_act_rows.append({"carbon_source":cs,"bnumber":gid,"active":int(active)})
    gene_df = pd.DataFrame(gene_act_rows)
    gene_pivot = gene_df.pivot_table(
        index="bnumber", columns="carbon_source", values="active"
    ).reset_index()
    gene_pivot.columns.name = None
    cs_cols = [c for c in cs_list if c in gene_pivot.columns]
    gene_pivot["n_conditions_on"]  = gene_pivot[cs_cols].sum(axis=1)
    gene_pivot["always_on"]        = (gene_pivot[cs_cols]==1).all(axis=1).astype(int)
    gene_pivot["always_off"]       = (gene_pivot[cs_cols]==0).all(axis=1).astype(int)
    gene_pivot["condition_variable"]= (~((gene_pivot[cs_cols]==1).all(axis=1)|
                                         (gene_pivot[cs_cols]==0).all(axis=1))).astype(int)

# ── Print key results ─────────────────────────────────────────────────────────
print("\n📊 Growth summary (top 20 by rFBA growth):")
top = growth_df.head(20)[["carbon_source","wt_growth","rfba_growth","pct_impact","n_applied"]]
print(top.to_string(index=False))

if len(off_df):
    print(f"\n📊 Reactions OFF across all carbon sources: {heatmap.shape[0]} unique")
    print(f"  OFF in all {len(growth_df)} conditions: {(heatmap['n_cs_off']==len(growth_df)).sum()}")
    print(f"  OFF in >half conditions:               {(heatmap['n_cs_off']>len(growth_df)/2).sum()}")
    print(f"  OFF in only 1 condition (unique):       {(heatmap['n_cs_off']==1).sum()}")

    print("\n📊 Top subsystems (reactions OFF aerobically on glucose):")
    glc = applied_df[applied_df.carbon_source=="Glucose"]
    for sub, n in glc.groupby("subsystem").size().sort_values(ascending=False).head(10).items():
        print(f"  {sub:<45}: {n}")

    print("\n📊 Gene regulatory states across all carbon sources:")
    print(f"  Always ON  (all conditions): {gene_pivot.always_on.sum()}")
    print(f"  Always OFF (all conditions): {gene_pivot.always_off.sum()}")
    print(f"  Condition-variable:          {gene_pivot.condition_variable.sum()}")

    print("\n📊 Most condition-variable OFF reactions (top 15):")
    cvar = heatmap[(heatmap.n_cs_off>0)&(heatmap.n_cs_off<len(growth_df)-2)]
    print(cvar[["reaction_name","subsystem","n_cs_off"]].head(15).to_string(index=False))

# ── Save ──────────────────────────────────────────────────────────────────────
print("\n" + "="*65)
growth_df.to_csv(OUT/"multi_carbon_growth_summary.csv", index=False)
off_df.to_csv(OUT/"multi_carbon_off_reactions.csv", index=False)
if len(off_df):
    heatmap.reset_index().rename(columns={"index":"reaction_id"}).to_csv(
        OUT/"multi_carbon_reaction_heatmap.csv", index=False)
    gene_pivot.to_csv(OUT/"multi_carbon_gene_activity.csv", index=False)

print(f"  ✓ multi_carbon_growth_summary.csv   ({len(growth_df)} substrates)")
print(f"  ✓ multi_carbon_off_reactions.csv    ({len(off_df)} rows)")
if len(off_df):
    print(f"  ✓ multi_carbon_reaction_heatmap.csv ({heatmap.shape[0]} reactions × {len(growth_df)} substrates)")
    print(f"  ✓ multi_carbon_gene_activity.csv    ({len(gene_pivot)} genes × {len(growth_df)} substrates)")
print("\nDone.")
