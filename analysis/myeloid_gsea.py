"""Preranked enrichment using public Enrichr snapshots of Hallmark and Reactome."""
import os
import tempfile
import sys
os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "myeloid-mpl"))
from pathlib import Path
import pandas as pd
import gseapy as gp

out = Path("Myeloid/disease_pathways")
sets = {}
for collection in ["MSigDB_Hallmark_2020", "Reactome_2022"]:
    for line in (out / (collection + ".gmt")).read_text().splitlines():
        fields = line.split("\t")
        sets[collection + " | " + fields[0]] = list(dict.fromkeys(fields[2:]))
jobs = []
within = pd.read_csv(out / "within_state_allgenes.csv")
jobs.extend(("COPD_" + key, frame) for key, frame in within.groupby("State"))
transition = pd.read_csv(out / "state_transition_ranks.csv")
jobs.extend((key, frame) for key, frame in transition.groupby("contrast"))
if "--base" not in sys.argv:
    if "--positions" in sys.argv:
        jobs = []
    positions = pd.read_csv(out / "disease_pseudotime_position_tests.csv")
    jobs.extend(("COPD_" + key, frame) for key, frame in positions.groupby("position"))
results = []
if "--positions" in sys.argv and (out / "gsea_all_results.csv").exists():
    results.append(pd.read_csv(out / "gsea_all_results.csv"))
for name, frame in jobs:
    print("GSEA", name, flush=True)
    rank = frame[["feature", "t"]].dropna().drop_duplicates("feature")
    rank = rank.sort_values(["t", "feature"], ascending=[False, True])
    result = gp.prerank(rnk=rank, gene_sets=sets, min_size=15, max_size=500,
                        permutation_num=1000, threads=4, seed=42, outdir=None,
                        no_plot=True, verbose=False)
    table = result.res2d.copy()
    table["contrast"] = name
    results.append(table)
    pd.concat(results, ignore_index=True).to_csv(out / "gsea_all_results.csv", index=False)
all_results = pd.concat(results, ignore_index=True)
antigen = all_results.Term.str.contains(
    "MHC|Antigen Presentation|Antigen [Pp]rocessing|Antigen processing|Interferon Gamma|Interferon-gamma",
    case=False, regex=True)
all_results.loc[antigen].to_csv(out / "gsea_antigen_results.csv", index=False)
print("Completed", len(jobs), "contrasts", flush=True)
