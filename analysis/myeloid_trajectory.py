"""Infer a monocyte/macrophage transcriptional continuum using Harmony and DPT."""
import os
import tempfile
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(tempfile.gettempdir(), "myeloid-numba"))
os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "myeloid-mpl"))
import sys
import json
import random
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr
import scanpy as sc
import anndata as ad
import igraph as ig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42})

prefix, output = sys.argv[1:3]
out = Path(output)
out.mkdir(exist_ok=True, parents=True)
meta = pd.read_csv(prefix + "_meta.csv").set_index("cell")
h = pd.read_csv(prefix + "_harmony.csv", index_col=0).loc[meta.index]
u = pd.read_csv(prefix + "_umap.csv", index_col=0).loc[meta.index]
selected = meta.myeloid_subtype.isin(["Monocyte", "Macrophage"])
m = meta.loc[selected].copy()
a = ad.AnnData(sp.csr_matrix((len(m), 1)), obs=m)
a.obsm["X_harmony"] = h.loc[m.index].iloc[:, :30].to_numpy(dtype=np.float32)
a.obsm["X_umap"] = u.loc[m.index].to_numpy(dtype=np.float32)
print("Building 30-neighbor graph for", a.n_obs, "cells", flush=True)
sc.pp.neighbors(a, n_neighbors=30, use_rep="X_harmony", random_state=42,
                transformer="pynndescent")
n_components, components = connected_components(a.obsp["connectivities"])
root_cluster = a.obs.myeloid_cluster.astype(str).eq("3").to_numpy()
root_candidates = np.flatnonzero(root_cluster)
centroid = np.median(a.obsm["X_harmony"][root_candidates], axis=0)
distance = np.sum((a.obsm["X_harmony"][root_candidates] - centroid) ** 2, axis=1)
root = root_candidates[np.argmin(distance)]
root_component = components[root]
in_component = components == root_component
excluded_components = a.obs.index[~in_component]
a = a[in_component].copy()
root_candidates = np.flatnonzero(a.obs.myeloid_cluster.astype(str).eq("3"))
distance = np.sum((a.obsm["X_harmony"][root_candidates] - centroid) ** 2, axis=1)
ranked = root_candidates[np.argsort(distance)]
root = int(ranked[0])
print("Diffusion map for", a.n_obs, "cells", flush=True)
sc.tl.diffmap(a, n_comps=15, random_state=42)
a.uns["iroot"] = root
sc.tl.dpt(a, n_dcs=10, n_branchings=0)
pseudotime = a.obs.dpt_pseudotime.to_numpy().copy()
# Test a nearby root from a different donor while holding the graph fixed.
root_donor = a.obs.iloc[root].myeloid_donor
alternative = next(int(i) for i in ranked if a.obs.iloc[i].myeloid_donor != root_donor)
a.uns["iroot"] = alternative
sc.tl.dpt(a, n_dcs=10, n_branchings=0)
rho = spearmanr(pseudotime, a.obs.dpt_pseudotime).statistic
a.obs["dpt_pseudotime"] = pseudotime
a.uns["iroot"] = root

random.seed(42)
ig.set_random_number_generator(random.Random(42))
g = ig.Graph.Weighted_Adjacency(a.obsp["connectivities"], mode="undirected", loops=False)
partition = g.community_leiden(objective_function="modularity", weights="weight",
                               resolution=0.15, n_iterations=4)
clusters = np.asarray(partition.membership)
median_time = pd.Series(pseudotime).groupby(clusters).median().sort_values()
mapping = {int(k): "State" + str(i + 1) for i, k in enumerate(median_time.index)}
states = [mapping[int(k)] for k in clusters]
a.obs["State"] = pd.Categorical(states, categories=list(mapping.values()), ordered=True)
sc.tl.paga(a, groups="State")
state_edges = pd.DataFrame(a.uns["paga"]["connectivities"].toarray(),
                           index=a.obs.State.cat.categories, columns=a.obs.State.cat.categories)
state_edges.to_csv(out / "state_graph_connectivity.csv")

result = meta.copy()
result["trajectory_scope"] = np.where(selected, "disconnected_from_root", "outside_monocyte_macrophage")
result["State"] = None
result["pseudotime"] = np.nan
result.loc[a.obs_names, "trajectory_scope"] = "monocyte_macrophage"
result.loc[a.obs_names, "State"] = a.obs.State.astype(str).to_numpy()
result.loc[a.obs_names, "pseudotime"] = pseudotime
result.to_csv(out / "cell_pseudotime_state.csv")
pd.crosstab(a.obs.State, a.obs.myeloid_subtype).to_csv(out / "state_celltypes.csv")
pd.crosstab(a.obs.State, a.obs.myeloid_donor).to_csv(out / "state_donors.csv")
pd.crosstab(a.obs.State, a.obs.status).to_csv(out / "state_disease_composition.csv")
a.write_h5ad(out / "trajectory_graph.h5ad", compression="gzip")
with open(out / "trajectory_parameters.json", "w") as f:
    json.dump({"method": "Diffusion pseudotime + Leiden graph states + PAGA",
               "representation": "existing Harmony dimensions 1-30", "neighbors": 30,
               "state_resolution": 0.15, "root_cell": a.obs_names[root],
               "root_basis": "cluster 3 FCN1+ monocyte nearest Harmony centroid",
               "alternative_root_cell": a.obs_names[alternative],
               "alternative_root_spearman": float(rho), "cells": a.n_obs,
               "graph_components": int(n_components),
               "cells_outside_root_component": len(excluded_components),
               "scanpy_version": sc.__version__}, f, indent=2)

fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)
xy = a.obsm["X_umap"]
p = axes[0].scatter(xy[:, 0], xy[:, 1], c=pseudotime, cmap="viridis", s=0.3, rasterized=True)
fig.colorbar(p, ax=axes[0], label="Diffusion pseudotime")
categories = a.obs.State.cat.categories
colors = plt.get_cmap("tab20")(np.linspace(0, 1, len(categories)))
centers = []
for state, color in zip(categories, colors):
    mask = a.obs.State.eq(state).to_numpy()
    axes[1].scatter(xy[mask, 0], xy[mask, 1], color=color, s=0.3, rasterized=True)
    center = np.median(xy[mask], axis=0)
    centers.append(center)
    axes[1].text(*center, state, fontsize=8, color="black")
for i in range(len(categories)):
    for j in range(i + 1, len(categories)):
        weight = state_edges.iloc[i, j]
        if weight >= 0.01:
            c = np.asarray([centers[i], centers[j]])
            axes[1].plot(c[:, 0], c[:, 1], color="black", linewidth=0.4 + 4 * weight, alpha=0.6)
for ax, title in zip(axes, ["Monocyte–macrophage pseudotime", "Graph States (PAGA edges >= 0.01)"]):
    ax.scatter(*xy[root], marker="x", s=20, c="black", zorder=5)
    ax.set(xlabel="Original UMAP 1", ylabel="Original UMAP 2", title=title)
    ax.set_aspect("equal")
for extension in ["pdf", "png"]:
    fig.savefig(out / ("trajectory_pseudotime_states." + extension), dpi=220)
print("Completed", len(categories), "states; root sensitivity rho =", rho, flush=True)
