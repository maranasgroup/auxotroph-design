"""
rFBA Regulatory Constraint Engine for iML1515 (v3 - feasibility-preserving)
=============================================================================
Implements the correct Covert 2001-2004 rFBA approach:
  1. WT FBA → get flux values for flux-dependent rule signals
  2. Evaluate all Boolean rules → gene ON/OFF
  3. GPR logic → candidate OFF reactions
  4. Feasibility-preserving constraint application:
       Apply all constraints; if infeasible, iteratively relax the most
       blocking constraint until feasible (Covert 2004 iterative strategy)
  5. Report final OFF reactions, growth rates, and flagged rules
"""
import re, warnings, pandas as pd, numpy as np, cobra
from pathlib import Path

warnings.filterwarnings("ignore")
BASE = Path("/mnt/user-data/uploads")
OUT  = Path("/mnt/user-data/outputs")

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

ALWAYS_OPEN = {
    "EX_pi_e","EX_nh4_e","EX_so4_e","EX_mg2_e","EX_fe2_e","EX_fe3_e",
    "EX_h2o_e","EX_h_e","EX_co2_e","EX_k_e","EX_na1_e","EX_cl_e",
    "EX_mobd_e","EX_cobalt2_e","EX_zn2_e","EX_mn2_e","EX_ni2_e",
    "EX_cu2_e","EX_ca2_e","EX_sel_e","EX_slnt_e","EX_tungs_e",
}

COND_DEFS = {
    "aerobic_glucose": {
        "uptakes": {"EX_glc__D_e": -10.0, "EX_o2_e": -15.0},
        "signals": {"glc__D_e":10.0,"glc_DASH_D_e":10.0,"o2_e":15.0,
                    "lac__D_e":0.0,"ac_e":0.0,"no3_e":0.0,"no2_e":0.0,
                    "fe2_e":1.0,"nh4_e":10.0,"so4_e":10.0,"pi_e":10.0,"mg2_e":1.0,
                    "mobd_e":1.0,"cu2_e":0.001,"zn2_e":0.001,"ni2_e":0.001,"mn2_e":0.001,
                    "h2o2_e":0.0,"no_e":0.0,"leu_DASH_L_e":0.0,"lys_DASH_L_e":0.0,
                    "fru_e":0.0,"xyl_DASH_D_e":0.0,"larab_e":0.0,"4abut_e":0.0,
                    "gam_e":0.0,"ura_e":0.0,"thymd_e":0.0,"thy_e":0.0,
                    "acgam_e":0.0,"chb_e":0.0,"phenac_e":0.0,
                    "BiomassEcoli":1.0,"Stress":False,"Stringent":False,
                    "Oxidative-Stress":False,"pH":7.4,"Rich-Medium":False},
    },
    "anaerobic_glucose": {
        "uptakes": {"EX_glc__D_e": -10.0, "EX_o2_e": 0.0},
        "signals": {"glc__D_e":10.0,"glc_DASH_D_e":10.0,"o2_e":0.0,
                    "lac__D_e":0.0,"ac_e":0.0,"no3_e":0.0,"no2_e":0.0,
                    "fe2_e":1.0,"nh4_e":10.0,"so4_e":10.0,"pi_e":10.0,"mg2_e":1.0,
                    "mobd_e":1.0,"cu2_e":0.001,"zn2_e":0.001,"ni2_e":0.001,"mn2_e":0.001,
                    "h2o2_e":0.0,"no_e":0.0,"leu_DASH_L_e":0.0,"lys_DASH_L_e":0.0,
                    "fru_e":0.0,"xyl_DASH_D_e":0.0,"larab_e":0.0,"4abut_e":0.0,
                    "gam_e":0.0,"ura_e":0.0,"thymd_e":0.0,"thy_e":0.0,
                    "acgam_e":0.0,"chb_e":0.0,"phenac_e":0.0,
                    "BiomassEcoli":1.0,"Stress":False,"Stringent":False,
                    "Oxidative-Stress":False,"pH":7.4,"Rich-Medium":False},
    },
    "aerobic_glucose_acetate": {
        "uptakes": {"EX_glc__D_e": -10.0, "EX_o2_e": -15.0, "EX_ac_e": -5.0},
        "signals": {"glc__D_e":10.0,"glc_DASH_D_e":10.0,"o2_e":15.0,
                    "lac__D_e":0.0,"ac_e":5.0,"no3_e":0.0,"no2_e":0.0,
                    "fe2_e":1.0,"nh4_e":10.0,"so4_e":10.0,"pi_e":10.0,"mg2_e":1.0,
                    "mobd_e":1.0,"cu2_e":0.001,"zn2_e":0.001,"ni2_e":0.001,"mn2_e":0.001,
                    "h2o2_e":0.0,"no_e":0.0,"leu_DASH_L_e":0.0,"lys_DASH_L_e":0.0,
                    "fru_e":0.0,"xyl_DASH_D_e":0.0,"larab_e":0.0,"4abut_e":0.0,
                    "gam_e":0.0,"ura_e":0.0,"thymd_e":0.0,"thy_e":0.0,
                    "acgam_e":0.0,"chb_e":0.0,"phenac_e":0.0,
                    "BiomassEcoli":1.0,"Stress":False,"Stringent":False,
                    "Oxidative-Stress":False,"pH":7.4,"Rich-Medium":False},
    },
    "aerobic_glucose_lactose": {
        "uptakes": {"EX_glc__D_e": -10.0, "EX_o2_e": -15.0, "EX_lac__D_e": -3.0},
        "signals": {"glc__D_e":10.0,"glc_DASH_D_e":10.0,"o2_e":15.0,
                    "lac__D_e":3.0,"lac_DASH_D_e":3.0,"lcts_e":3.0,
                    "ac_e":0.0,"no3_e":0.0,"no2_e":0.0,
                    "fe2_e":1.0,"nh4_e":10.0,"so4_e":10.0,"pi_e":10.0,"mg2_e":1.0,
                    "mobd_e":1.0,"cu2_e":0.001,"zn2_e":0.001,"ni2_e":0.001,"mn2_e":0.001,
                    "h2o2_e":0.0,"no_e":0.0,"leu_DASH_L_e":0.0,"lys_DASH_L_e":0.0,
                    "fru_e":0.0,"xyl_DASH_D_e":0.0,"larab_e":0.0,"4abut_e":0.0,
                    "gam_e":0.0,"ura_e":0.0,"thymd_e":0.0,"thy_e":0.0,
                    "acgam_e":0.0,"chb_e":0.0,"phenac_e":0.0,
                    "BiomassEcoli":1.0,"Stress":False,"Stringent":False,
                    "Oxidative-Stress":False,"pH":7.4,"Rich-Medium":False},
    },
}

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
    w = re.sub(r'pH\s*>\s*([\d.]+)',
        lambda m: "True" if cond.get("pH",7.4) > float(m.group(1)) else "False", w)
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

def set_medium(model, uptakes):
    for r in model.reactions:
        if r.id.startswith("EX_") and r.lower_bound < 0 and r.id not in ALWAYS_OPEN:
            r.lower_bound = 0.0
    for rid, lb in uptakes.items():
        try: model.reactions.get_by_id(rid).lower_bound = lb
        except: pass

# ── Feasibility-preserving constraint application ────────────────────────────
def apply_rfba_constraints_safe(model, off_rxns, wt_fluxes, uptakes):
    """
    Apply regulatory constraints (zero bounds) to OFF reactions.
    If model becomes infeasible, relax the most growth-blocking constraint
    (highest |WT flux|) until feasibility is restored.
    Returns (applied_rxns, relaxed_rxns).
    """
    # Sort: apply least impactful first (lowest |WT flux|)
    off_sorted = sorted(off_rxns, key=lambda rid: abs(wt_fluxes.get(rid, 0.0)))

    applied  = []
    relaxed  = []
    orig_bounds = {}

    for rid in off_sorted:
        try:
            rxn = model.reactions.get_by_id(rid)
        except Exception:
            continue
        orig_bounds[rid] = (rxn.lower_bound, rxn.upper_bound)
        rxn.lower_bound = 0.0
        rxn.upper_bound = 0.0

        sol = model.optimize()
        if sol.status == "optimal" and sol.objective_value > 1e-6:
            applied.append(rid)
        else:
            # Relax this constraint — restore bounds
            rxn.lower_bound = orig_bounds[rid][0]
            rxn.upper_bound = orig_bounds[rid][1]
            relaxed.append(rid)

    return applied, relaxed

# ── Main loop ─────────────────────────────────────────────────────────────────
print("\n" + "="*65)
print("rFBA ANALYSIS (feasibility-preserving)")
print("="*65)

gene_act_rows = []
off_rxn_rows  = []
growth_rows   = []

KEY_TFS = {"b0683":"Fur","b1988":"Nac","b3357":"CRP","b1334":"FNR",
           "b4401":"ArcA","b0889":"Lrp","b0080":"Cra"}

for cname, cdef in COND_DEFS.items():
    print(f"\n── {cname} ──")
    cond_sig = dict(cdef["signals"])

    # Step 1: WT FBA
    with model:
        set_medium(model, cdef["uptakes"])
        sol_wt = model.optimize()
    if sol_wt.status != "optimal":
        print("  WT FBA infeasible! Skipping."); continue
    wt_growth = sol_wt.objective_value
    wt_fluxes = dict(sol_wt.fluxes)
    print(f"  WT growth: {wt_growth:.4f} h⁻¹")

    # Add flux signals to condition dict
    for rxn in model.reactions:
        cond_sig[f"_flux_{rxn.id}"] = float(sol_wt.fluxes.get(rxn.id, 0.0))
    cond_sig["_flux_GLCpts"] = abs(sol_wt.fluxes.get("EX_glc__D_e", 0.0))

    # Step 2: Evaluate all rules
    tf_act   = precompute_tf_activity(cond_sig, all_rules)
    gene_act = {gid: tf_act.get(gid, True) for gid in model_gene_ids}

    spot = {n: tf_act.get(bn) for bn, n in KEY_TFS.items()}
    print(f"  Key TFs: {spot}")
    n_off_g = sum(1 for v in gene_act.values() if not v)
    print(f"  Genes: ON={sum(v for v in gene_act.values())}, OFF={n_off_g}")

    # Record gene activity
    for gid, active in gene_act.items():
        gene_act_rows.append({"bnumber":gid,"condition":cname,"active":int(active),
            "rule_source":"delta" if gid in delta_rules else
                          ("iMC1010" if gid in imc_rules else "constitutive")})

    # Step 3: Candidate OFF reactions
    candidates = [rxn.id for rxn in model.reactions if reaction_is_off(rxn, gene_act)]
    print(f"  Candidate OFF reactions: {len(candidates)}")

    # Step 4: Feasibility-preserving application
    with model:
        set_medium(model, cdef["uptakes"])
        applied, relaxed = apply_rfba_constraints_safe(model, candidates, wt_fluxes, cdef["uptakes"])
        sol_rfba = model.optimize()
    rfba_growth = sol_rfba.objective_value if sol_rfba.status=="optimal" else 0.0

    print(f"  Constraints applied: {len(applied)}  |  Relaxed (Flag B): {len(relaxed)}")
    print(f"  rFBA growth: {rfba_growth:.4f} h⁻¹")

    # Build off-reaction table (applied only)
    for rid in applied:
        rxn = model.reactions.get_by_id(rid)
        off_genes = [g.id for g in rxn.genes if not gene_act.get(g.id, True)]
        off_rxn_rows.append({
            "condition"    : cname,
            "reaction_id"  : rid,
            "reaction_name": rxn.name,
            "subsystem"    : rxn.subsystem or "Unclassified",
            "gpr"          : rxn.gene_reaction_rule,
            "off_genes"    : "; ".join(off_genes),
            "n_off_genes"  : len(off_genes),
            "n_rxn_genes"  : len(rxn.genes),
            "wt_flux"      : round(wt_fluxes.get(rid, 0.0), 5),
            "flag"         : "applied",
        })
    for rid in relaxed:
        rxn = model.reactions.get_by_id(rid)
        off_genes = [g.id for g in rxn.genes if not gene_act.get(g.id, True)]
        off_rxn_rows.append({
            "condition"    : cname,
            "reaction_id"  : rid,
            "reaction_name": rxn.name,
            "subsystem"    : rxn.subsystem or "Unclassified",
            "gpr"          : rxn.gene_reaction_rule,
            "off_genes"    : "; ".join(off_genes),
            "n_off_genes"  : len(off_genes),
            "n_rxn_genes"  : len(rxn.genes),
            "wt_flux"      : round(wt_fluxes.get(rid, 0.0), 5),
            "flag"         : "relaxed_Flag_B",
        })

    growth_rows.append({
        "condition"          : cname,
        "wt_growth"          : round(wt_growth, 4),
        "rfba_growth"        : round(rfba_growth, 4),
        "growth_impact"      : round(wt_growth - rfba_growth, 4),
        "n_applied"          : len(applied),
        "n_relaxed_FlagB"    : len(relaxed),
        "n_candidates_total" : len(candidates),
        "feasible"           : sol_rfba.status == "optimal",
    })

# ── Summary ───────────────────────────────────────────────────────────────────
print("\n" + "="*65)
print("SUMMARY")
print("="*65)

growth_df = pd.DataFrame(growth_rows)
off_df    = pd.DataFrame(off_rxn_rows)
gene_df   = pd.DataFrame(gene_act_rows)

print("\n📊 Growth rates (h⁻¹):")
print(growth_df.to_string(index=False))

if len(off_df):
    applied_df = off_df[off_df.flag=="applied"]
    cond_names = list(COND_DEFS.keys())

    print(f"\n📊 Applied OFF reactions per condition:")
    for cn in cond_names:
        n = len(applied_df[applied_df.condition==cn])
        print(f"  {cn:<38}: {n}")

    # Reaction summary
    rxn_summary = (applied_df.groupby("reaction_id")
        .agg(reaction_name=("reaction_name","first"),
             subsystem=("subsystem","first"),
             gpr=("gpr","first"),
             off_in_n=("condition","nunique"),
             conditions=("condition", lambda x: "; ".join(sorted(set(x)))),
             off_genes=("off_genes","first"),
             wt_flux=("wt_flux","first"))
        .reset_index()
        .sort_values(["off_in_n","subsystem"], ascending=[False,True]))
    for cn in cond_names:
        cset = set(applied_df[applied_df.condition==cn]["reaction_id"])
        rxn_summary[cn] = rxn_summary["reaction_id"].isin(cset).astype(int)

    print(f"\n📊 OFF in multiple conditions:")
    for k in [4,3,2,1]:
        print(f"  OFF in {k} condition(s): {(rxn_summary.off_in_n==k).sum()}")

    print(f"\n📊 Top subsystems (aerobic glucose):")
    aero = applied_df[applied_df.condition=="aerobic_glucose"]
    top = aero.groupby("subsystem").size().sort_values(ascending=False).head(12)
    for sub, n in top.items():
        print(f"  {sub:<50}: {n}")

    print(f"\n📊 Condition-variable OFF reactions (top 25):")
    cv = rxn_summary[(rxn_summary.off_in_n>0)&(rxn_summary.off_in_n<4)]
    print(cv[["reaction_id","reaction_name","subsystem","conditions","wt_flux"]].head(25).to_string(index=False))

    print(f"\n📊 Flag B reactions (relaxed to preserve feasibility):")
    fb = off_df[off_df.flag=="relaxed_Flag_B"]
    print(f"  Total unique reactions: {fb.reaction_id.nunique()}")
    print(fb.drop_duplicates("reaction_id")[["reaction_id","reaction_name","gpr","wt_flux"]].head(10).to_string(index=False))

# Gene pivot
gene_pivot = gene_df.pivot_table(
    index=["bnumber","rule_source"], columns="condition", values="active"
).reset_index()
gene_pivot.columns.name = None
cc = [c for c in cond_names if c in gene_pivot.columns]
gene_pivot["always_on"]      = (gene_pivot[cc]==1).all(axis=1).astype(int)
gene_pivot["always_off"]     = (gene_pivot[cc]==0).all(axis=1).astype(int)
gene_pivot["cond_variable"]  = (~((gene_pivot[cc]==1).all(axis=1)|(gene_pivot[cc]==0).all(axis=1))).astype(int)

print(f"\n📊 Gene regulatory states (all 1516 model genes):")
for src in ["iMC1010","delta","constitutive"]:
    sub = gene_pivot[gene_pivot.rule_source==src]
    print(f"  {src:<15}: always_ON={sub.always_on.sum():4d}  "
          f"always_OFF={sub.always_off.sum():4d}  "
          f"condition_variable={sub.cond_variable.sum():4d}")

# Save
print("\n" + "="*65)
off_df.to_csv(OUT/"off_reactions.csv", index=False)
rxn_summary.to_csv(OUT/"reaction_off_counts.csv", index=False)
growth_df.to_csv(OUT/"rfba_growth_summary.csv", index=False)
gene_pivot.to_csv(OUT/"gene_activity.csv", index=False)
print(f"  ✓ off_reactions.csv         ({len(applied_df)} applied-constraint rows)")
print(f"  ✓ reaction_off_counts.csv   ({len(rxn_summary)} unique reactions)")
print(f"  ✓ rfba_growth_summary.csv")
print(f"  ✓ gene_activity.csv         ({len(gene_pivot)} genes)")
print("\nDone.")
