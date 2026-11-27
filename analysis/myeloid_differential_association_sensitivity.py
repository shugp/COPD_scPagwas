"""Small-sample checks for differential candidate-pathway associations."""
import os
import tempfile
os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "myeloid-mpl"))
from pathlib import Path
import itertools
import math
import json
import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out = Path("Myeloid/disease_pathways")
d = pd.read_csv(out / "differential_association_input.csv")
scores = pd.read_csv(out / "antigen_donor_state_scores.csv")
tests = pd.read_csv(out / "candidate_pathway_differential_association.csv")
candidates = list(d.columns[5:])
pathways = list(scores.pathway.unique())
rng = np.random.default_rng(42)

def partial_r(x, y, source):
    design = pd.get_dummies(pd.Series(source), dtype=float).to_numpy()
    rank = np.linalg.matrix_rank(design)
    xx = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    yy = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    divisor = np.sqrt(np.dot(xx, xx) * np.dot(yy, yy))
    r = np.clip(np.dot(xx, yy) / divisor, -.999999, .999999) if divisor > 1e-10 else np.nan
    return r, len(x) - rank - 2

def partial_corr_matrix(w, x, y, source_indices):
    covariance = (w @ (x[:, :, None] * y[:, None, :]).reshape(len(x), -1)).reshape(
        len(w), x.shape[1], y.shape[1])
    vx, vy = w @ (x*x), w @ (y*y)
    for idx in source_indices:
        n = w[:, idx].sum(axis=1)
        sx, sy = w[:, idx] @ x[idx], w[:, idx] @ y[idx]
        denominator = np.where(n > 0, n, 1)
        covariance -= sx[:, :, None] * sy[:, None, :] / denominator[:, None, None]
        vx -= sx*sx / denominator[:, None]
        vy -= sy*sy / denominator[:, None]
    denominator = np.sqrt(np.maximum(vx[:, :, None], 0) * np.maximum(vy[:, None, :], 0))
    r = np.divide(covariance, denominator, out=np.full_like(covariance, np.nan), where=denominator>1e-10)
    return np.clip(r, -.999999, .999999)

permutation_rows = []
for state in sorted(d.State.unique()):
    data = d.loc[d.State.eq(state)].reset_index(drop=True)
    pivot = scores.pivot(index="sample", columns="pathway", values="score")
    x = data[candidates].to_numpy()
    y = pivot.loc[data["sample"], pathways].to_numpy()
    observed = data.status.eq("COPD").to_numpy()
    strata = [np.flatnonzero(data.data_source.eq(source)) for source in data.data_source.unique()]
    n_cases = [int(observed[idx].sum()) for idx in strata]
    choices = [math.comb(len(idx), n) for idx, n in zip(strata, n_cases)]
    total = math.prod(choices)
    exact = total <= 10000
    if exact:
        possibilities = [list(itertools.combinations(idx, n)) for idx, n in zip(strata, n_cases)]
        masks = np.zeros((total, len(data)), dtype=float)
        for row, assignment in enumerate(itertools.product(*possibilities)):
            masks[row, list(itertools.chain.from_iterable(assignment))] = 1
    else:
        masks = np.zeros((10000, len(data)), dtype=float)
        for row in range(len(masks)):
            for idx, n in zip(strata, n_cases):
                masks[row, rng.choice(idx, n, replace=False)] = 1
    case = partial_corr_matrix(masks, x, y, strata)
    control = partial_corr_matrix(1-masks, x, y, strata)
    statistics = np.abs(np.arctanh(case)-np.arctanh(control))
    case0 = partial_corr_matrix(observed[None, :].astype(float), x, y, strata)[0]
    control0 = partial_corr_matrix((~observed)[None, :].astype(float), x, y, strata)[0]
    statistic0 = np.abs(np.arctanh(case0)-np.arctanh(control0))
    valid = np.isfinite(statistics)
    exceed = (statistics >= statistic0[None, :, :] - 1e-10) & valid
    numerator, denominator = exceed.sum(axis=0), valid.sum(axis=0)
    if not exact:
        numerator, denominator = numerator+1, denominator+1
    p = np.divide(numerator, denominator, out=np.full_like(numerator, np.nan, dtype=float),
                  where=denominator>0)
    for i, gene in enumerate(candidates):
        for j, path in enumerate(pathways):
            permutation_rows.append({"State": state, "feature": gene, "pathway": path,
              "conditional_permutation_p": p[i, j], "permutation_assignments": len(masks),
              "permutation_exact": exact, "permutations_valid": int(valid[:, i, j].sum())})
permutations = pd.DataFrame(permutation_rows)
valid = permutations.conditional_permutation_p.notna()
permutations["FDR_conditional_permutation"] = np.nan
permutations.loc[valid, "FDR_conditional_permutation"] = multipletests(
    permutations.loc[valid, "conditional_permutation_p"], method="fdr_bh")[1]
permutations.to_csv(out / "differential_association_conditional_permutations.csv", index=False)
keys = ["State", "feature", "pathway"]
tests = tests.drop(columns=[c for c in permutations.columns if c not in keys and c in tests.columns])
tests = tests.merge(permutations, on=keys, validate="one_to_one")
tests.to_csv(out / "candidate_pathway_differential_association.csv", index=False)

# Check the two Fisher-screen discoveries with donor removal and stratified bootstrap.
loo_rows, bootstrap_rows = [], []
for _, row in tests.loc[tests.FDR_partial.lt(.05)].iterrows():
    data = d.loc[d.State.eq(row.State)].merge(
        scores.loc[scores.pathway.eq(row.pathway), ["sample", "score"]], on="sample")
    def evaluate(frame):
        rs, dfs = {}, {}
        for status in ["COPD", "Healthy"]:
            g = frame.loc[frame.status.eq(status)]
            rs[status], dfs[status] = partial_r(g[row.feature].to_numpy(), g.score.to_numpy(), g.data_source.to_numpy())
        delta = rs["COPD"]-rs["Healthy"]
        p = np.nan
        if min(dfs.values()) > 0 and all(np.isfinite(list(rs.values()))):
            z = (np.arctanh(rs["COPD"])-np.arctanh(rs["Healthy"])) / np.sqrt(1/dfs["COPD"]+1/dfs["Healthy"])
            p = 2*norm.sf(abs(z))
        return delta, p, rs
    for donor in data.donor:
        delta, p, rs = evaluate(data.loc[data.donor.ne(donor)])
        loo_rows.append({"State": row.State, "feature": row.feature, "pathway": row.pathway,
                         "excluded_donor": donor, "delta_r": delta, "p_Fisher": p,
                         "r_COPD": rs["COPD"], "r_Healthy": rs["Healthy"]})
    groups = [g for _, g in data.groupby(["status", "data_source"], sort=False)]
    deltas = []
    for repeat in range(5000):
        frame = pd.concat([g.iloc[rng.integers(0, len(g), len(g))] for g in groups], ignore_index=True)
        delta, _, _ = evaluate(frame)
        if np.isfinite(delta):
            deltas.append(delta)
    low, high = np.quantile(deltas, [.025, .975])
    bootstrap_rows.append({"State": row.State, "feature": row.feature, "pathway": row.pathway,
                           "delta_r": row.partial_delta_r, "bootstrap_low": low, "bootstrap_high": high,
                           "valid_bootstraps": len(deltas), "bootstraps": 5000})
pd.DataFrame(loo_rows).to_csv(out / "differential_association_leave_one_donor.csv", index=False)
pd.DataFrame(bootstrap_rows).to_csv(out / "differential_association_bootstrap.csv", index=False)

# Directly contrast both groups and distinguish screening from robustness evidence.
show = tests.sort_values("partial_p").head(6).copy()
labels = show.feature + " / " + show.State + " / " + show.pathway.map(
    lambda p: "MHC-II" if "2132295" in p else ("Cross-presentation" if "1236975" in p else
              ("IFN-gamma" if "877300" in p else "MHC-I")))
fig, ax = plt.subplots(figsize=(10, 4.8))
for j, (_, row) in enumerate(show.iterrows()):
    for status, offset, color in [("Healthy", -.12, "#4575A5"), ("COPD", .12, "#B65B4C")]:
        r, low, high = (row["partial_"+status+"_"+field] for field in ["r","low","high"])
        ax.errorbar(r, j+offset, xerr=[[r-low], [high-r]], fmt="o", color=color,
                    markersize=4, capsize=2, linewidth=1, label=status if j == 0 else None)
    text = "Fisher q=" + format(row.FDR_partial, ".3g") + "; perm q=" + format(
        row.FDR_conditional_permutation, ".3g")
    ax.text(1.04, j, text, transform=ax.get_yaxis_transform(), va="center", fontsize=8)
ax.set_yticks(range(len(labels)), labels)
ax.invert_yaxis()
ax.axvline(0, color="gray", linewidth=.8)
ax.set(xlim=(-1.05,1.05), xlabel="Source-adjusted partial correlation (95% Fisher CI)")
ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0, -.12), ncol=2)
fig.tight_layout()
for extension in ["png","pdf"]:
    fig.savefig(out/("differential_association_group_correlations."+extension), dpi=180, bbox_inches="tight")
plt.close(fig)
print(tests.sort_values("partial_p")[["State","feature","FDR_partial","conditional_permutation_p",
      "FDR_conditional_permutation","FDR_partial_spearman","FDR_HC3_interaction"]].head(8).to_string(index=False))
print(pd.DataFrame(bootstrap_rows).to_string(index=False))
print("Significant conditional permutation tests:", int(tests.FDR_conditional_permutation.lt(.05).sum()))

parameters_path = out / "differential_association_parameters.json"
parameters = json.loads(parameters_path.read_text())
parameters["conditional_permutation"] = {
    "method": "Status labels permuted within data source, preserving source-specific group counts",
    "exact_if_assignments_at_most": 10000, "otherwise_draws": 10000,
    "interpretation": "Sensitivity test under within-source label exchangeability",
    "significant_after_BH": int(tests.FDR_conditional_permutation.lt(.05).sum())}
parameters["bootstrap"] = "5000 donor resamples within status and source for Fisher-screen hits; exploratory unadjusted intervals"
parameters_path.write_text(json.dumps(parameters, indent=2))
