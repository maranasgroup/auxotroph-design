# Auxotrophic Strain Design via Conditional Essentiality and Rescue Analysis

Pipeline for designing auxotrophic *E. coli* strains: regulatory FBA across 45 carbon sources, single/double/triple knockout screens, and ResMut rescue via a multi-key (EC / BiGG / MetaNetX) reaction library.


## Layout

```
data/        iML1515.xml, iMC1010.csv, delta_gene_rules_normalized_v2.csv, carbon_sources.csv
scripts/     01_rfba_multi_carbon.py → 02_build_reaction_library.py → 03/04/05_*_knockout.py
src/         knockout_utils.py, visualization.py
outputs/     generated CSVs (gitignored)
```

## Install

```bash
git clone https://github.com/<user>/auxotroph-design.git && cd auxotroph-design
pip install -r requirements.txt
```

Python 3.10, COBRApy 0.29, GLPK. Download `iML1515.xml` from [BiGG](http://bigg.ucsd.edu/models) into `data/`.

## Run

```bash
python scripts/01_rfba_multi_carbon.py        
python scripts/02_build_reaction_library.py   
python scripts/03_single_knockout.py
python scripts/04_double_knockout.py
python scripts/05_triple_knockout.py
```

Each stage reads the previous stage's CSVs from `outputs/`. Total runtime ~6–10 h on 4 cores (stage 04 dominates).

`scripts/rfba_constraints.py` runs the rFBA validation on 4 hand-picked conditions; independent of the main pipeline.

## Key outputs

| File | Contents |
|------|----------|
| `multi_carbon_off_reactions.csv` | regulatory OFF reactions per substrate |
| `{single,double,triple}_ko_conditionally_essential.csv` | CE designs per substrate |
| `{single,double,triple}_ko_rescue_dedup.csv` | rescue functions per design |


## Parameters

Hardcoded at the top of `src/knockout_utils.py` — these are the values that produced the published figures.

