"""Direct COPD-versus-Healthy comparison of candidate-pathway coupling within State."""
import os
import tempfile
os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "myeloid-mpl"))
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import norm, pearsonr, rankdata
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
import patsy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

out = Path("Myeloid/disease_pathways")
observed = pd.read_csv(out / "differential_association_input.csv")
scores = pd.read_csv(out / "antigen_donor_state_scores.csv")
candidate = list(observed.columns[5:])
pathways = list(scores.pathway.unique())
states = sorted(observed.State.unique())
colors = {"Healthy": "#4575A5", "COPD": "#B65B4C"}
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
                     "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42})
short = {
    "R-HSA-983170": "MHC-I peptide loading", "R-HSA-1236975": "Cross-presentation",
    "R-HSA-877300": "IFN-gamma signaling", "R-HSA-2132295": "MHC-II presentation"}
def path_label(path):
    return next(value for key, value in short.items() if key in path)
def correlation(d, gene, ranked=False, adjusted=True):
    a = d[gene].to_numpy(dtype=float)
    b = d.score.to_numpy(dtype=float)
    if ranked:
        a, b = rankdata(a), rankdata(b)
    design = np.ones((len(d), 1))
    if adjusted:
        design = pd.get_dummies(d.data_source, dtype=float).to_numpy()
    rank = np.linalg.matrix_rank(design)
    a = a - design @ np.linalg.lstsq(design, a, rcond=None)[0]
    b = b - design @ np.linalg.lstsq(design, b, rcond=None)[0]
    variance_a, variance_b = np.dot(a, a), np.dot(b, b)
    denominator_df = len(d) - rank - 2
    valid = denominator_df > 0 and variance_a > 1e-10 and variance_b > 1e-10
    if not valid:
        return {"r": np.nan, "low": np.nan, "high": np.nan, "df": denominator_df,
                "n": len(d), "source_rank": rank}, a, b
    r = float(np.clip(np.dot(a, b) / np.sqrt(variance_a * variance_b), -.999999, .999999))
    z, se = np.arctanh(r), 1 / np.sqrt(denominator_df)
    return {"r": r, "low": np.tanh(z - 1.96 * se), "high": np.tanh(z + 1.96 * se),
            "df": denominator_df, "n": len(d), "source_rank": rank}, a, b

rows, scatter = [], []
for state in states:
    for path in pathways:
        d = observed[observed.State.eq(state)].merge(
            scores.loc[scores.pathway.eq(path), ["sample", "score"]], on="sample", validate="one_to_one")
        for gene in candidate:
            row = {"State": state, "pathway": path, "feature": gene}
            for label, ranked, adjusted in [("partial", False, True), ("raw", False, False),
                                             ("partial_spearman", True, True)]:
                results = {}
                for status in ["Healthy", "COPD"]:
                    group = d[d.status.eq(status)]
                    result, a, b = correlation(group, gene, ranked=ranked, adjusted=adjusted)
                    results[status] = result
                    row.update({label + "_" + status + "_" + key: value for key, value in result.items()})
                    if label == "partial":
                        scatter.extend({"State": state, "pathway": path, "feature": gene,
                                        "donor": donor, "status": status,
                                        "residual_expression": aa, "residual_score": bb}
                                       for donor, aa, bb in zip(group.donor, a, b))
                case, control = results["COPD"], results["Healthy"]
                row[label + "_delta_r"] = case["r"] - control["r"]
                z = np.nan
                if np.isfinite(case["r"]) and np.isfinite(control["r"]):
                    z = (np.arctanh(case["r"]) - np.arctanh(control["r"])) / np.sqrt(
                        1/case["df"] + 1/control["df"])
                row[label + "_z"] = z
                row[label + "_p"] = 2 * norm.sf(abs(z))
            # Independent donors within State. HC3 protects against heteroskedasticity.
            # A source represented by one donor has leverage one and contributes no
            # within-source slope information. Remove it before HC3 estimation.
            source_n = d.data_source.map(d.data_source.value_counts())
            regression_data = d.loc[source_n >= 2].copy()
            regression_data["expression"] = regression_data[gene]
            y, X = patsy.dmatrices("score ~ C(data_source) + C(status, Treatment('Healthy')) * expression",
                                   regression_data, return_type="dataframe")
            term = "C(status, Treatment('Healthy'))[T.COPD]:expression"
            row["HC3_slope_difference"] = row["HC3_interaction_p"] = np.nan
            if np.linalg.matrix_rank(X) == X.shape[1] and len(X) - X.shape[1] >= 3:
                fit = sm.OLS(y, X).fit(cov_type="HC3", use_t=True)
                row["HC3_slope_difference"] = float(fit.params[term])
                row["HC3_interaction_p"] = float(fit.pvalues[term])
                row["HC3_n_donors"] = len(X)
            rows.append(row)
tests = pd.DataFrame(rows)
for p_column, q_column in [("partial_p", "FDR_partial"), ("raw_p", "FDR_raw"),
                           ("partial_spearman_p", "FDR_partial_spearman"),
                           ("HC3_interaction_p", "FDR_HC3_interaction")]:
    tests[q_column] = np.nan
    valid = tests[p_column].notna() & np.isfinite(tests[p_column])
    tests.loc[valid, q_column] = multipletests(tests.loc[valid, p_column], method="fdr_bh")[1]
tests["test_status"] = np.where(tests.partial_p.notna(), "tested", "insufficient_residual_variation_or_df")
tests = tests.sort_values(["FDR_partial", "partial_p"], na_position="last")
tests.to_csv(out / "candidate_pathway_differential_association.csv", index=False)
scatter = pd.DataFrame(scatter)
scatter.to_csv(out / "differential_association_residuals.csv", index=False)

fig, axes = plt.subplots(2, 2, figsize=(10, 11), constrained_layout=True)
for ax, path in zip(axes.flat, pathways):
    d = tests[tests.pathway.eq(path)]
    matrix = d.pivot(index="feature", columns="State", values="partial_delta_r").reindex(
        index=candidate, columns=states)
    marks = d.assign(mark=np.where(d.FDR_partial < .05, "*", "")).pivot(
        index="feature", columns="State", values="mark").reindex_like(matrix)
    sns.heatmap(matrix, cmap="RdBu_r", vmin=-2, vmax=2, center=0,
                annot=marks, fmt="", ax=ax, cbar_kws={"label": "r(COPD) - r(Healthy)"})
    ax.set(title=path_label(path) + " (* Fisher FDR < 0.05)", xlabel=None, ylabel=None)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)
fig.savefig(out / "differential_association_heatmap.png", dpi=180, bbox_inches="tight")
fig.savefig(out / "differential_association_heatmap.pdf", bbox_inches="tight")
plt.close(fig)

# Show all primary discoveries, or the strongest exploratory contrasts when none pass FDR.
chosen = tests.loc[tests.partial_p.notna()].head(6)
fig, axes = plt.subplots(2, 3, figsize=(11, 7), constrained_layout=True)
for ax, (_, row) in zip(axes.flat, chosen.iterrows()):
    d = scatter[scatter.State.eq(row.State) & scatter.pathway.eq(row.pathway) & scatter.feature.eq(row.feature)]
    for status, color in colors.items():
        group = d[d.status.eq(status)]
        ax.scatter(group.residual_expression, group.residual_score, color=color, s=14, alpha=.8,
                   label=status + " r=" + format(row["partial_" + status + "_r"], ".2f"))
        if group.residual_expression.std() > 1e-8:
            slope, intercept = np.polyfit(group.residual_expression, group.residual_score, 1)
            grid = np.linspace(group.residual_expression.min(), group.residual_expression.max(), 40)
            ax.plot(grid, intercept+slope*grid, color=color, linewidth=1)
    ax.set(title=row.feature + " / " + row.State + "\n" + path_label(row.pathway) +
                 " (q=" + format(row.FDR_partial, ".3g") + ")",
           xlabel="Expression, source-adjusted", ylabel="Pathway score, source-adjusted")
    ax.legend(frameon=False, fontsize=7)
fig.savefig(out / "differential_association_scatter.png", dpi=180, bbox_inches="tight")
fig.savefig(out / "differential_association_scatter.pdf", bbox_inches="tight")
plt.close(fig)
summary = {"primary_test": "Within-State Fisher-z difference of partial Pearson correlations",
           "nuisance": "Data-source intercepts fitted separately within COPD and Healthy",
           "unit": "Independent biological donor within each State",
           "candidate_self_inclusion": "All 15 coloc candidates excluded from pathway scores",
           "multiple_testing": "BH across all estimable candidate x pathway x State tests",
           "secondary": ["unadjusted Pearson difference", "source-adjusted partial Spearman difference",
                         "source-adjusted disease x expression slope, HC3 robust SE"],
           "tests_total": len(tests), "tests_estimable": int(tests.partial_p.notna().sum()),
           "significant_primary": int(tests.FDR_partial.lt(.05).sum()),
           "significant_slope": int(tests.FDR_HC3_interaction.lt(.05).sum()),
           "small_samples": "Fisher-z intervals and tests are approximate, especially with few residual degrees of freedom"}
(out / "differential_association_parameters.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
print(tests[["State", "feature", "pathway", "partial_Healthy_r", "partial_COPD_r",
             "partial_delta_r", "partial_p", "FDR_partial", "FDR_HC3_interaction"]].head(12).to_string(index=False))
