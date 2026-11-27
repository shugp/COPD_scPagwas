"""Plain scientific figures for current and previous Myeloid analyses."""
import os
import tempfile
os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "myeloid-mpl"))
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.dpi": 200})
root = Path("Myeloid")
case_colors = {"Healthy": "#4575A5", "COPD": "#B65B4C"}
rng = np.random.default_rng(42)

def save(fig, stem, png=True):
    fig.savefig(str(stem) + ".pdf", bbox_inches="tight")
    if png:
        fig.savefig(str(stem) + ".png", bbox_inches="tight")
    plt.close(fig)

def scatter_umap(xy, labels, title="", annotate=True):
    fig, ax = plt.subplots(figsize=(6, 5))
    categories = sorted(pd.Series(labels).dropna().unique(), key=str)
    palette = plt.get_cmap("tab10" if len(categories) <= 10 else "tab20")
    for i, name in enumerate(categories):
        mask = np.asarray(labels == name)
        ax.scatter(xy[mask, 0], xy[mask, 1], s=0.35, color=palette(i % palette.N),
                   rasterized=True, label=str(name), linewidths=0)
        if annotate:
            center = np.median(xy[mask], axis=0)
            ax.text(*center, str(name), ha="center", fontsize=8, color="black")
    ax.set(xlabel="UMAP 1", ylabel="UMAP 2", title=title, aspect="equal")
    if not annotate:
        ax.legend(frameon=False, markerscale=3, loc="center left", bbox_to_anchor=(1, .5))
    return fig

def dotplot(table, gene, group, value, percentage, stem, width=10):
    genes = list(dict.fromkeys(table[gene]))
    groups = list(dict.fromkeys(table[group]))[::-1]
    fig, ax = plt.subplots(figsize=(width, max(2.6, len(groups) * .4 + 1.3)))
    x = table[gene].map({g: i for i, g in enumerate(genes)})
    y = table[group].map({g: i for i, g in enumerate(groups)})
    points = ax.scatter(x, y, s=table[percentage] * .9, c=table[value],
                        cmap="Blues", edgecolors="none")
    ax.set_xticks(range(len(genes)), genes, rotation=60, ha="right")
    ax.set_yticks(range(len(groups)), groups)
    ax.set_xlim(-.7, len(genes) - .3)
    fig.colorbar(points, ax=ax, fraction=.025, pad=.02, label="Mean log1p expression")
    for size in [25, 50, 100]:
        ax.scatter([], [], s=size * .9, c="gray", label=str(size) + "%")
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.14, 1), title="Expressing")
    save(fig, stem)

def facet_expression(table, group, stem, y="log2CPM", color_by_group=False):
    genes = list(dict.fromkeys(table.feature))
    groups = list(dict.fromkeys(table[group]))
    cols = 3
    fig, axes = plt.subplots(int(np.ceil(len(genes)/cols)), cols,
                             figsize=(10, int(np.ceil(len(genes)/cols))*2.2), squeeze=False)
    for ax, gene in zip(axes.flat, genes):
        d = table[table.feature.eq(gene)]
        values = [d.loc[d[group].eq(g), y].dropna().to_numpy() for g in groups]
        ax.boxplot(values, widths=.5, showfliers=False,
                   medianprops={"color": "black"}, boxprops={"color": "gray"},
                   whiskerprops={"color": "gray"}, capprops={"color": "gray"})
        for i, (g, v) in enumerate(zip(groups, values)):
            color = case_colors.get(g, "gray") if color_by_group else "gray"
            ax.scatter(i + 1 + rng.uniform(-.12, .12, len(v)), v, s=4, color=color, alpha=.6)
        ax.set_xticks(range(1, len(groups)+1), groups, rotation=45 if len(groups)>2 else 0)
        ax.set_title(gene)
        ax.set_ylabel("log2 CPM")
    for ax in list(axes.flat)[len(genes):]:
        ax.set_visible(False)
    fig.tight_layout()
    save(fig, stem, png=False)

clean = root / "results_cleaned"
trajectory = root / "trajectory"
disease = root / "disease_pathways"
meta = pd.read_csv(trajectory / "cell_pseudotime_state.csv")
xy = pd.read_csv(clean / "umap_coordinates.csv", index_col=0).loc[meta.cell].to_numpy()
fig = scatter_umap(xy, meta.myeloid_subtype, annotate=False)
save(fig, clean / "myeloid_cleaned_umap")
dotplot(pd.read_csv(clean / "annotation_marker_expression.csv"),
        "gene", "subtype", "mean_log1p_expression", "percent_expressing",
        clean / "annotation_markers_dotplot", width=12)
facet_expression(pd.read_csv(clean / "coloc_overall_donor_expression.csv"), "status",
                 clean / "coloc_overall_expression", color_by_group=True)

active = meta.trajectory_scope.eq("monocyte_macrophage")
a = meta.loc[active].copy()
points = xy[active]
fig, axes = plt.subplots(1, 2, figsize=(10, 4.6))
p = axes[0].scatter(points[:, 0], points[:, 1], c=a.pseudotime, cmap="viridis",
                    s=.4, rasterized=True, linewidths=0)
fig.colorbar(p, ax=axes[0], label="Pseudotime", fraction=.04)
states = sorted(a.State.unique())
palette = plt.get_cmap("tab10")
for i, state in enumerate(states):
    mask = a.State.eq(state).to_numpy()
    axes[1].scatter(points[mask, 0], points[mask, 1], s=.4, color=palette(i),
                    rasterized=True, linewidths=0)
    center = np.median(points[mask], axis=0)
    axes[1].text(*center, state, ha="center", fontsize=8, color="black")
for ax, title in zip(axes, ["Diffusion pseudotime", "States"]):
    ax.set(xlabel="UMAP 1", ylabel="UMAP 2", title=title, aspect="equal")
fig.tight_layout()
save(fig, trajectory / "trajectory_pseudotime_states")
heat = pd.read_csv(trajectory / "candidate_state_contrasts.csv")
heat["State"] = heat.comparison.str.split(" vs ").str[0]
matrix = heat.pivot(index="feature", columns="State", values="logFC")
stars = heat.assign(mark=np.where(heat.FDR_all_candidate_contrasts < .05, "*", "")).pivot(
    index="feature", columns="State", values="mark").loc[matrix.index, matrix.columns]
fig, ax = plt.subplots(figsize=(5.8, 6))
sns.heatmap(matrix, cmap="RdBu_r", center=0, ax=ax, annot=stars, fmt="",
            cbar_kws={"label": "log2FC vs State1"}, annot_kws={"size": 9})
ax.set(xlabel=None, ylabel=None)
save(fig, trajectory / "candidate_state_heatmap")
facet_expression(pd.read_csv(trajectory / "candidate_state_donor_expression.csv"),
                 "State", trajectory / "candidate_state_expression")
dotplot(pd.read_csv(trajectory / "state_marker_expression.csv"),
        "gene", "State", "mean_log1p", "percent_expressing",
        trajectory / "state_markers", width=8)

# Descriptive mean profiles of donor bins; formal inference is in the model tables.
curve = pd.read_csv(trajectory / "candidate_pseudotime_donor_expression.csv")
fig, axes = plt.subplots(5, 3, figsize=(10, 11))
for ax, (gene, d) in zip(axes.flat, curve.groupby("feature", sort=True)):
    ax.scatter(d.pseudotime, d.log2CPM, s=3, color="gray", alpha=.4)
    bins = pd.qcut(d.pseudotime, 10, duplicates="drop")
    means = d.groupby(bins, observed=True)[["pseudotime", "log2CPM"]].mean()
    ax.plot(means.pseudotime, means.log2CPM, color="black", linewidth=1)
    ax.set(title=gene, xlabel="Pseudotime", ylabel="log2 CPM")
fig.tight_layout()
save(fig, trajectory / "candidate_pseudotime_curves", png=False)

# New disease and antigen-presentation figures.
within = pd.read_csv(disease / "within_state_candidates.csv")
matrix = within.pivot(index="feature", columns="State", values="logFC")
stars = within.assign(mark=np.where(within.FDR_candidates_all_states < .05, "*", "")).pivot(
    index="feature", columns="State", values="mark").loc[matrix.index, matrix.columns]
fig, ax = plt.subplots(figsize=(6, 6))
sns.heatmap(matrix, cmap="RdBu_r", center=0, annot=stars, fmt="", ax=ax,
            cbar_kws={"label": "COPD - Healthy log2FC"})
ax.set(xlabel=None, ylabel=None)
save(fig, disease / "within_state_candidate_effects")
scores = pd.read_csv(disease / "antigen_donor_state_scores.csv")
paths = scores.pathway.unique()
fig, axes = plt.subplots(len(paths), 1, figsize=(8, 2.3*len(paths)), squeeze=False)
short_names = {
    "R-HSA-2132295": "MHC-II presentation", "R-HSA-983170": "MHC-I peptide loading",
    "R-HSA-1236975": "Cross-presentation", "R-HSA-877300": "IFN-gamma signaling"}
for ax, path in zip(axes.flat, paths):
    d = scores[scores.pathway.eq(path)]
    sns.boxplot(data=d, x="State", y="score", hue="status", hue_order=["Healthy", "COPD"],
                palette=case_colors, ax=ax, fliersize=0, linewidth=.7)
    sns.stripplot(data=d, x="State", y="score", hue="status", hue_order=["Healthy", "COPD"],
                  dodge=True, jitter=.08, palette=case_colors, size=2, alpha=.7,
                  ax=ax, legend=False)
    title = next((name for key, name in short_names.items() if key in path), path)
    ax.set(title=title, xlabel=None, ylabel="Relative score")
    ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
save(fig, disease / "antigen_state_disease_scores")
association = pd.read_csv(disease / "candidate_antigen_associations.csv")
association["path"] = association.pathway.map(lambda path: next(
    (name for key, name in short_names.items() if key in path), path))
matrix = association.pivot(index="feature", columns="path", values="slope")
stars = association.assign(mark=np.where(association.FDR < .05, "*", "")).pivot(
    index="feature", columns="path", values="mark").loc[matrix.index, matrix.columns]
fig, ax = plt.subplots(figsize=(7, 6))
sns.heatmap(matrix, cmap="RdBu_r", center=0, annot=stars, fmt="", ax=ax,
            cbar_kws={"label": "Adjusted score / log2 CPM"})
ax.set(xlabel=None, ylabel=None)
ax.set_xticklabels(ax.get_xticklabels(), rotation=35, ha="right")
save(fig, disease / "candidate_antigen_associations")

observed = pd.read_csv(disease / "disease_pseudotime_observations.csv")
focus = ["HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DRB1", "B2M", "ERAP2",
         "CD74", "CIITA", "CTSS", "TAP1", "TAP2", "NLRC5"]
fig, axes = plt.subplots(4, 3, figsize=(10, 9))
for ax, gene in zip(axes.flat, focus):
    d = observed[observed.feature.eq(gene)]
    for status, color in case_colors.items():
        z = d[d.status.eq(status)]
        ax.scatter(z.pseudotime, z.log2CPM, s=3, color=color, alpha=.35)
        bins = pd.cut(z.pseudotime, 8)
        means = z.groupby(bins, observed=True)[["pseudotime", "log2CPM"]].mean()
        ax.plot(means.pseudotime, means.log2CPM, color=color, linewidth=1, label=status)
    ax.set(title=gene, xlabel="Pseudotime", ylabel="log2 CPM")
axes[0, 0].legend(frameon=False)
fig.tight_layout()
save(fig, disease / "antigen_disease_pseudotime")

gsea_path = disease / "gsea_antigen_results.csv"
if gsea_path.exists():
    gsea = pd.read_csv(gsea_path)
    gsea = gsea[~gsea.Term.str.contains("B Cell|Nef", case=False)]
    gsea["short"] = gsea.Term.str.split(" | ", regex=False).str[-1].str.replace(
        r" R-HSA-\d+$", "", regex=True)
    replacements = {
        "Antigen Presentation: Folding, Assembly, Peptide Loading Of Class I MHC": "MHC-I peptide loading",
        "MHC Class II Antigen Presentation": "MHC-II presentation",
        "Class I MHC Mediated Antigen Processing And Presentation": "MHC-I processing / presentation",
        "Antigen processing-Cross Presentation": "Cross-presentation",
        "Antigen Processing: Ubiquitination And Proteasome Degradation": "Antigen proteasomal processing",
        "Interferon Gamma Signaling": "IFN-gamma signaling",
        "Interferon Gamma Response": "IFN-gamma response (Hallmark)"}
    gsea["short"] = gsea["short"].replace(replacements)
    matrix = gsea.pivot_table(index="short", columns="contrast", values="NES")
    sig = gsea.assign(mark=np.where(gsea["FDR q-val"] < .05, "*", "")).pivot_table(
        index="short", columns="contrast", values="mark", aggfunc="first").reindex_like(matrix)
    fig, ax = plt.subplots(figsize=(11, max(3, .35*len(matrix))))
    sns.heatmap(matrix, cmap="RdBu_r", center=0, annot=sig, fmt="", ax=ax,
                cbar_kws={"label": "GSEA NES"})
    ax.set(xlabel=None, ylabel=None)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=55, ha="right")
    save(fig, disease / "antigen_gsea")
print("Completed plain scientific plots")
