"""
build_reaction_library.py
=========================
Builds the digital reaction library for the ResMut rescue mechanism.

What "22,689 keys found but no iML1515 match" meant and why it was wrong
-------------------------------------------------------------------------
The old code required external reactions to share a key (EC / BiGG / MNX)
with an iML1515 reaction before including them.  That was an unnecessary
pre-filter that caused false negatives.

At rescue time, the lookup is:
  "find library reactions sharing a key with the KNOCKED-OUT reaction"
The knocked-out reaction IS an iML1515 reaction, so the rescue lookup
already enforces specificity.  A library reaction only needs to have ANY
key at all — the FBA test itself determines if it actually rescues growth.

The fix: include ALL external reactions that have at least one extractable
key (EC, BiGG, or MNX).  Reactions with no keys at all cannot be looked up
and are still skipped.

Metabolite namespace compatibility
-----------------------------------
The BiGG universal model uses the same BiGG metabolite IDs as iML1515
(atp_c, fad_c, coa_c, …).  Reactions with metabolites not present in
iML1515 (e.g. mitochondrial compartment atp_m) will fail gracefully when
added to the model — COBRApy raises a KeyError which is caught.

Supported external model formats
---------------------------------
  BiGG universal model     notes.original_bigg_ids, empty gene_rule
  MetaNetX universal model EC / BIGG / KEGG uppercase keys, MNXR* ids
  Standard organism GEMs   bigg.reaction, ec-code, metanetx.reaction
"""

import json, warnings, io, re, pandas as pd, cobra
from pathlib import Path
from collections import defaultdict

warnings.filterwarnings("ignore")
from pathlib import Path
REPO_ROOT = Path(__file__).resolve().parent      # or .parent.parent if scripts live in src/
DATA      = REPO_ROOT / "data"
OUT       = REPO_ROOT / "outputs"
OUT.mkdir(exist_ok=True)
BASE = DATA


# ── Annotation key aliases ────────────────────────────────────────────────────
EC_ALIASES   = {"ec-code", "ec", "EC", "brenda", "BRENDA"}
BIGG_ALIASES = {"bigg.reaction", "bigg", "BIGG", "biggid", "BiGG"}
MNX_ALIASES  = {"metanetx.reaction", "metanetx", "mnx", "MNX", "MNXID"}
RHEA_ALIASES = {"rhea", "RHEA"}
SKIP_PREFIXES = ("EX_", "DM_", "SK_", "sink_", "BIOMASS", "biomass")

EC_PATTERN = re.compile(r'^\d+\.\d+\.\d+[\.\d\-]*$')

def _is_ec(val: str) -> bool:
    return bool(EC_PATTERN.match(val.strip()))

def _norm_annotation(ann) -> dict:
    if isinstance(ann, dict):   return ann
    if isinstance(ann, list):
        result = {}
        for item in ann:
            if isinstance(item, (list, tuple)) and len(item) == 2:
                k, v = item
                if isinstance(k, str) and k:
                    result[k] = v if isinstance(v, list) else str(v)
        return result
    return {}

def _normalise_model_dict(d: dict) -> dict:
    for coll in ("reactions", "metabolites", "genes"):
        for obj in d.get(coll, []):
            if "annotation" in obj:
                obj["annotation"] = _norm_annotation(obj["annotation"])
    return d

def load_model_safe(path: Path) -> cobra.Model | None:
    try:
        if path.suffix == ".json":
            with open(path, "r", encoding="utf-8") as fh:
                raw = json.load(fh)
            _normalise_model_dict(raw)
            buf = io.StringIO(json.dumps(raw))
            return cobra.io.load_json_model(buf)
        else:
            return cobra.io.read_sbml_model(str(path))
    except Exception as e:
        print(f"  Failed to load {path.name}: {e}")
        return None

def get_keys(rxn: cobra.Reaction) -> dict:
    """
    Extract EC, BiGG, MetaNetX keys from any model format.
    Priority order for BiGG: annotation → notes.original_bigg_ids → rxn.id
    """
    ann   = rxn.annotation
    notes = getattr(rxn, "notes", {}) or {}

    def _collect(aliases):
        vals = []
        for alias in aliases:
            v = ann.get(alias, [])
            if isinstance(v, str):            v = [v]
            elif not isinstance(v, list):     v = [str(v)]
            vals.extend(v)
        return [str(v).strip() for v in vals if str(v).strip() not in ("", "nan")]

    # EC — validate X.X.X.X pattern to exclude KEGG R-numbers
    ecs = sorted({v for v in _collect(EC_ALIASES) if _is_ec(v)})

    # BiGG — annotation → notes.original_bigg_ids → reaction id
    bigg_raw = _collect(BIGG_ALIASES)
    bigg = bigg_raw[0] if bigg_raw else ""
    if not bigg:
        orig = notes.get("original_bigg_ids", [])
        if isinstance(orig, list) and orig:   bigg = str(orig[0]).strip()
        elif isinstance(orig, str) and orig:  bigg = orig.strip()
    if not bigg and not rxn.id.startswith(("MNXR", "R0")):
        bigg = rxn.id   # reaction id IS the BiGG id in universal models

    # MetaNetX — annotation → MNXR* reaction id fallback
    mnx_raw = _collect(MNX_ALIASES)
    mnx = mnx_raw[0] if mnx_raw else ""
    if not mnx and rxn.id.startswith("MNXR"):
        mnx = rxn.id

    rhea = sorted(set(_collect(RHEA_ALIASES)))
    return {"ec": ecs, "bigg": bigg, "mnx": mnx, "rhea": rhea}


# ── Load host model ───────────────────────────────────────────────────────────
print("Loading iML1515...")
host     = cobra.io.read_sbml_model(str(BASE / "iML1515.xml"))
host_ids = {r.id for r in host.reactions}
print(f"  {len(host.reactions)} reactions | {len(host.genes)} genes")

# ── Build multi-key index for iML1515 ────────────────────────────────────────
host_index = defaultdict(list)
for r in host.reactions:
    if not r.gene_reaction_rule:
        continue
    k = get_keys(r)
    for ec in k["ec"]:
        host_index[("ec",   ec)].append(r.id)
    if k["bigg"]:
        host_index[("bigg", k["bigg"])].append(r.id)
    if k["mnx"]:
        host_index[("mnx",  k["mnx"])].append(r.id)

print(f"  Host index: {len(host_index)} (key_type, value) pairs")


# ── Layer 1: iML1515 internal ─────────────────────────────────────────────────
rows = []
for r in host.reactions:
    if not r.gene_reaction_rule:
        continue
    k = get_keys(r)
    if not (k["ec"] or k["bigg"] or k["mnx"]):
        continue

    alts = set()
    for ec in k["ec"]:
        alts |= set(host_index[("ec",   ec)])
    if k["bigg"]:
        alts |= set(host_index[("bigg", k["bigg"])])
    if k["mnx"]:
        alts |= set(host_index[("mnx",  k["mnx"])])
    alts.discard(r.id)

    rows.append({
        "layer": 1, "source": "iML1515",
        "organism": "E. coli K-12 MG1655",
        "reaction_id": r.id, "reaction_name": r.name,
        "ec_numbers": "; ".join(k["ec"]),
        "bigg_id": k["bigg"], "metanetx_id": k["mnx"],
        "rhea_ids": "; ".join(k["rhea"]),
        "gene_rule": r.gene_reaction_rule,
        "lower_bound": r.lower_bound, "upper_bound": r.upper_bound,
        "n_internal_alts": len(alts),
        "internal_alts": "; ".join(sorted(alts)),
    })

print(f"\nLayer 1: {len(rows)} iML1515 gene reactions indexed")


# ── Layer 2: external model(s) ───────────────────────────────────────────────
EXCLUDED  = {"iML1515.xml"}
ext_total = 0

for f in sorted(BASE.iterdir()):
    if f.suffix not in {".xml", ".json", ".sbml"} or f.name in EXCLUDED:
        continue

    print(f"\nLoading external model: {f.name}...")
    ext = load_model_safe(f)
    if ext is None:
        continue
    print(f"  {len(ext.reactions)} reactions loaded")

    # Diagnose annotation schema
    ann_keys = set()
    has_notes_bigg = False
    for r in list(ext.reactions)[:200]:
        ann_keys |= set(r.annotation.keys())
        if r.notes and r.notes.get("original_bigg_ids"):
            has_notes_bigg = True
    print(f"  Annotation keys (first 200): {sorted(ann_keys) or 'none'}")
    print(f"  Uses notes.original_bigg_ids: {has_notes_bigg}")

    added = already_in = no_keys = 0

    for r in ext.reactions:
        if any(r.id.startswith(p) for p in SKIP_PREFIXES):
            continue
        if r.id in host_ids:
            already_in += 1
            continue

        k = get_keys(r)

        # ── THE FIX ──────────────────────────────────────────────────────────
        # OLD (wrong): skip if no key matches iML1515 host_index
        # NEW (correct): include if reaction has ANY identifiable key
        # Rationale: the rescue lookup at runtime already enforces key-match
        # specificity against the knocked-out reaction.  Pre-filtering against
        # iML1515 host reactions caused 22,689 valid rescues to be excluded.
        # ─────────────────────────────────────────────────────────────────────
        if not (k["ec"] or k["bigg"] or k["mnx"]):
            no_keys += 1
            continue   # truly unidentifiable — cannot participate in rescue lookup

        rows.append({
            "layer": 2, "source": f.stem, "organism": f.stem,
            "reaction_id": r.id, "reaction_name": r.name,
            "ec_numbers": "; ".join(k["ec"]),
            "bigg_id": k["bigg"], "metanetx_id": k["mnx"],
            "rhea_ids": "; ".join(k["rhea"]),
            "gene_rule": r.gene_reaction_rule or "",
            "lower_bound": r.lower_bound, "upper_bound": r.upper_bound,
            "n_internal_alts": 0, "internal_alts": "",
        })
        added += 1
        ext_total += 1

    print(f"  Added {added} reactions to library")
    print(f"  Skipped: {already_in} already in iML1515 | "
          f"{no_keys} had no extractable keys (truly unidentifiable)")


# ── Save ──────────────────────────────────────────────────────────────────────
lib_df = pd.DataFrame(rows)
lib_df.to_csv(OUT / "reaction_library.csv", index=False)

print(f"\n{'='*55}")
print(f"Library summary:")
print(f"  Layer 1 (iML1515 internal): {(lib_df.layer==1).sum()}")
print(f"  Layer 2 (external models):  {(lib_df.layer==2).sum()}")
print(f"  Total:                      {len(lib_df)}")
if (lib_df.layer==2).sum() > 0:
    for src, grp in lib_df[lib_df.layer==2].groupby("source"):
        print(f"    {src}: {len(grp)} reactions")
print(f"\n  ✓ reaction_library.csv saved")