"""Disease-group trajectory, antigen pathways, and residual coexpression panels."""
import os
import tempfile
os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "myeloid-mpl"))
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False})
out = Path("Myeloid/disease_pathways")
colors = {"Healthy": "#4575A5", "COPD": "#B65B4C"}
rng = np.random.default_rng(42)
def save(fig, name):
    fig.savefig(out / (name + ".png"), dpi=180, bbox_inches="tight")
    fig.savefig(out / (name + ".pdf"), bbox_inches="tight")
    plt.close(fig)

meta = pd.read_csv("Myeloid/trajectory/cell_pseudotime_state.csv").set_index("cell")
a = meta.loc[meta.trajectory_scope.eq("monocyte_macrophage")].copy()
xy = pd.read_csv("Myeloid/results_cleaned/umap_coordinates.csv", index_col=0).loc[a.index].to_numpy()
status = a.groupby("myeloid_donor").status.first()
counts = pd.crosstab(a.myeloid_donor, a.State)
proportions = counts.div(counts.sum(axis=1), axis=0)
proportions.index.name = "donor"
proportions.to_csv(out / "donor_state_proportions.csv")
rows = []
for state in proportions.columns:
    case = proportions.loc[status.eq("COPD"), state]
    control = proportions.loc[status.eq("Healthy"), state]
    rows.append({"comparison": state + " proportion", "mean_COPD": case.mean(),
                 "mean_Healthy": control.mean(), "p": mannwhitneyu(case, control).pvalue})
medians = a.groupby("myeloid_donor").pseudotime.median()
case, control = medians[status.eq("COPD")], medians[status.eq("Healthy")]
rows.append({"comparison": "donor median pseudotime", "mean_COPD": case.mean(),
             "mean_Healthy": control.mean(), "p": mannwhitneyu(case, control).pvalue})
tests = pd.DataFrame(rows)
tests["FDR"] = multipletests(tests.p, method="fdr_bh")[1]
tests["n_COPD"], tests["n_Healthy"] = int(status.eq("COPD").sum()), int(status.eq("Healthy").sum())
tests.to_csv(out / "trajectory_disease_distribution_tests.csv", index=False)

fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
for ax, group in zip(axes[0], ["Healthy", "COPD"]):
    mask = a.status.eq(group).to_numpy()
    p = ax.scatter(xy[mask, 0], xy[mask, 1], c=a.pseudotime.to_numpy()[mask],
                   s=.4, cmap="viridis", vmin=0, vmax=1, linewidths=0, rasterized=True)
    ax.set(xlabel="UMAP 1", ylabel="UMAP 2", title=group + " (cells = " + str(mask.sum()) + ")",
           xlim=(xy[:, 0].min()-1, xy[:, 0].max()+1),
           ylim=(xy[:, 1].min()-1, xy[:, 1].max()+1), aspect="equal")
fig.colorbar(p, ax=list(axes[0]), fraction=.025, label="Pseudotime")
grid = np.linspace(0, 1, 200)
for group, color in colors.items():
    curves = []
    for donor in status.index[status.eq(group)]:
        values = np.sort(a.loc[a.myeloid_donor.eq(donor), "pseudotime"].to_numpy())
        curves.append(np.searchsorted(values, grid, side="right") / len(values))
    axes[1, 0].plot(grid, np.mean(curves, axis=0), color=color, linewidth=1.5,
                    label=group + " (" + str(len(curves)) + " donors)")
axes[1, 0].set(xlabel="Pseudotime", ylabel="Mean donor cumulative fraction",
               title="Equal weight per donor")
axes[1, 0].legend(frameon=False)
for j, state in enumerate(proportions.columns):
    for group, offset in [("Healthy", -.12), ("COPD", .12)]:
        values = proportions.loc[status.eq(group), state].to_numpy() * 100
        axes[1, 1].scatter(j + offset + rng.uniform(-.03, .03, len(values)), values,
                           color=colors[group], s=6, alpha=.6)
        axes[1, 1].errorbar(j+offset, values.mean(), yerr=values.std(ddof=1)/np.sqrt(len(values)),
                            fmt="_", color="black", markersize=8, capsize=2, linewidth=.8)
axes[1, 1].set_xticks(range(len(proportions.columns)), proportions.columns, rotation=35)
axes[1, 1].set(ylabel="Within-donor cells (%)", title="State proportions (mean +/- SE)")
save(fig, "trajectory_disease_comparison")

# GSEA and independent candidate-excluded pathway scores.
gsea = pd.read_csv(out / "gsea_all_results.csv")
definitions = {
    "MHC-I peptide loading": "Antigen Presentation: Folding, Assembly, Peptide Loading Of Class I MHC R-HSA-983170",
    "MHC-II presentation": "MHC Class II Antigen Presentation R-HSA-2132295",
    "Cross-presentation": "Antigen processing-Cross Presentation R-HSA-1236975",
    "IFN-gamma signaling": "Interferon Gamma Signaling R-HSA-877300",
    "IFN-gamma response": "Interferon Gamma Response"}
selected = []
for label, term in definitions.items():
    d = gsea[gsea.Term.str.split(" | ", regex=False).str[-1].eq(term) &
             gsea.contrast.str.match(r"COPD_State[1-6]$")].copy()
    d["pathway"] = label
    d["State"] = d.contrast.str.replace("COPD_", "", regex=False)
    selected.append(d)
selected = pd.concat(selected)
matrix = selected.pivot(index="pathway", columns="State", values="NES").reindex(definitions)
marks = selected.assign(mark=np.where(selected["FDR q-val"] < .05, "*", "")).pivot(
    index="pathway", columns="State", values="mark").reindex_like(matrix)
score = pd.read_csv(out / "antigen_score_disease_tests.csv")
score = score.loc[score.State.eq("State5")].copy()
score["path"] = score.pathway.map(lambda name: next(
    (label for label, term in definitions.items() if term == name), name))
score = score.set_index("path").reindex(list(definitions)[:-1])
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1.3, 1]},
                         constrained_layout=True)
sns.heatmap(matrix, cmap="RdBu_r", center=0, annot=marks, fmt="", ax=axes[0],
             cbar_kws={"label": "NES"}, annot_kws={"size": 10})
axes[0].set(title="COPD vs Healthy GSEA (* FDR < 0.05)", xlabel=None, ylabel=None)
for i, (name, row) in enumerate(score.iterrows()):
    axes[1].errorbar(row.effect_COPD, i, xerr=[[row.effect_COPD-row.ci_low], [row.ci_high-row.effect_COPD]],
                    fmt="o", color="black", markersize=4, capsize=2)
    axes[1].text(1.03, i, "q=" + format(row.FDR, ".3f"), transform=axes[1].get_yaxis_transform(),
                  ha="left", va="center", fontsize=8)
axes[1].axvline(0, color="gray", linewidth=.8)
axes[1].set_yticks(range(len(score)), score.index)
axes[1].invert_yaxis()
axes[1].set(xlabel="COPD - Healthy score (95% CI)",
            title="State5 score, candidates excluded")
save(fig, "antigen_pathway_evidence")

# Actual residual coexpression, accounting for donor and graph State.
cor = pd.read_csv(out / "antigen_residual_coexpression.csv")
candidate_order = ["HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DRB1", "B2M", "ERAP2"]
partner_order = ["CIITA", "CD74", "HLA-DMA", "HLA-DMB", "CTSS", "IFI30", "TAP1", "TAP2",
                 "TAPBP", "HLA-A", "HLA-B", "HLA-C", "NLRC5", "PSMB8", "PSMB9", "STAT1", "IRF1"]
matrix = cor.pivot(index="candidate", columns="partner", values="partial_correlation").reindex(
    index=candidate_order, columns=partner_order)
marks = cor.assign(mark=np.where(cor.FDR < .05, "*", "")).pivot(
    index="candidate", columns="partner", values="mark").reindex_like(matrix)
module = pd.read_csv(out / "candidate_module_reactome_enrichment.csv")
module_paths = {k: v for k, v in definitions.items() if k in ["MHC-I peptide loading", "MHC-II presentation", "IFN-gamma signaling"]}
module_rows = []
for label, name in module_paths.items():
    d = module[module.pathway.eq(name) & module.candidate.isin(candidate_order)].copy()
    d["path"] = label
    module_rows.append(d)
module = pd.concat(module_rows)
module["value"] = -np.log10(np.maximum(module.FDR, 1e-15))
m = module.pivot(index="candidate", columns="path", values="value").reindex(index=candidate_order, columns=module_paths)
sig = module.assign(mark=np.where(module.FDR < .05, "*", "")).pivot(
    index="candidate", columns="path", values="mark").reindex_like(m)
fig, axes = plt.subplots(1, 2, figsize=(12, 4.1), gridspec_kw={"width_ratios": [2.2, 1]}, constrained_layout=True)
sns.heatmap(matrix, cmap="RdBu_r", vmin=-1, vmax=1, center=0, annot=marks, fmt="", ax=axes[0],
            cbar_kws={"label": "Partial r"}, annot_kws={"size": 8})
axes[0].set(title="Coexpression, donor + State adjusted", xlabel=None, ylabel=None)
axes[0].set_yticklabels(axes[0].get_yticklabels(), rotation=0)
axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=60, ha="right")
sns.heatmap(m, cmap="Blues", annot=sig, fmt="", ax=axes[1],
            cbar_kws={"label": "-log10 FDR"}, annot_kws={"size": 8})
axes[1].set(title="Top-100 module enrichment", xlabel=None, ylabel=None)
axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=35, ha="right")
save(fig, "antigen_coexpression_evidence")
print("Completed disease comparison, pathway and residual coexpression panels")
