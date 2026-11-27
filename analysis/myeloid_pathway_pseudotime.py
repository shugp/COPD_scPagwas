"""Plot group-specific pathway trajectories and check interaction with donor fixed effects."""
import os
import tempfile
os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "myeloid-mpl"))
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out = Path("Myeloid/disease_pathways/pathway_pseudotime")
design = pd.read_csv(out / "donor_bin_design.csv")
observations = pd.read_csv(out / "antigen_pathway_donor_bin_observations.csv")
curves = pd.read_csv(out / "antigen_pathway_fitted_curves.csv")
differences = pd.read_csv(out / "antigen_pathway_curve_differences.csv")
tests = pd.read_csv(out / "antigen_pathway_disease_pseudotime_tests.csv")
short = {"R-HSA-983170": "MHC-I peptide loading", "R-HSA-2132295": "MHC-II presentation",
         "R-HSA-1236975": "Cross-presentation", "R-HSA-877300": "IFN-gamma signaling"}
def label(path):
    return next(value for key, value in short.items() if key in path)
rows = []
for path in tests.pathway:
    d = observations.loc[observations.pathway.eq(path)].merge(
        design[["sample", "b1", "b2", "b3"]], on="sample", validate="one_to_one")
    donor = pd.get_dummies(d.donor, dtype=float).to_numpy()
    basis = d[["b1", "b2", "b3"]].to_numpy()
    interaction = basis * d.status.eq("COPD").to_numpy()[:, None]
    for source_slopes in [False, True]:
        columns = [donor, basis]
        if source_slopes:
            source = pd.get_dummies(d.source, drop_first=True, dtype=float).to_numpy()
            columns.append((source[:, :, None] * basis[:, None, :]).reshape(len(d), -1))
        X = np.column_stack(columns + [interaction])
        assert np.linalg.matrix_rank(X) == X.shape[1]
        fit = sm.OLS(d.score.to_numpy(), X).fit(
            cov_type="cluster", cov_kwds={"groups": d.donor, "use_correction": True},
            use_t=True)
        L = np.zeros((3, X.shape[1]))
        L[:, -3:] = np.eye(3)
        test = fit.f_test(L)
        rows.append({"pathway": path, "sensitivity": "source_specific_curves" if source_slopes else "donor_fixed_effects",
                     "F": float(np.asarray(test.fvalue).squeeze()),
                     "p": float(np.asarray(test.pvalue).squeeze()), "numDF": test.df_num,
                     "denDF": test.df_denom, "n_donors": d.donor.nunique()})
sensitivity = pd.DataFrame(rows)
sensitivity["FDR"] = sensitivity.groupby("sensitivity").p.transform(
    lambda p: multipletests(p, method="fdr_bh")[1])
sensitivity.to_csv(out / "antigen_pathway_shape_sensitivity.csv", index=False)
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
                     "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42})
colors = {"Healthy": "#4575A5", "COPD": "#B65B4C"}
source_checks = pd.read_csv(out / "antigen_source_curve_sensitivity.csv") if (out / "antigen_source_curve_sensitivity.csv").exists() else None
fig, axes = plt.subplots(2, 2, figsize=(10, 7.3), constrained_layout=True)
for ax, (_, row) in zip(axes.flat, tests.iterrows()):
    for status, color in colors.items():
        d = curves[curves.pathway.eq(row.pathway) & curves.status.eq(status)]
        ax.plot(d.pseudotime, d.fitted_score, color=color, linewidth=1.3, label=status)
        ax.fill_between(d.pseudotime, d.ci_low, d.ci_high, color=color, alpha=.10, linewidth=0)
    obs = observations[observations.pathway.eq(row.pathway)]
    ax.plot(obs.pseudotime.to_numpy(), np.full(len(obs), ax.get_ylim()[0]),
            "|", color="gray", alpha=.3, markersize=3)
    ax.set(title=label(row.pathway) + "\nShape q=" + format(row.FDR_shape, ".3g") +
                 "; source-curve q=" + (format(source_checks.set_index("pathway").loc[row.pathway, "FDR_shape"], ".3g")
                                         if source_checks is not None else "NA"),
           xlabel="Pseudotime", ylabel="Relative pathway score")
axes[0, 0].legend(frameon=False)
for ext in ["png", "pdf"]:
    fig.savefig(out / ("antigen_pathway_disease_curves." + ext), dpi=180, bbox_inches="tight")
plt.close(fig)
fig, axes = plt.subplots(2, 2, figsize=(10, 7.3), constrained_layout=True)
for ax, (_, row) in zip(axes.flat, tests.iterrows()):
    d = differences[differences.pathway.eq(row.pathway)]
    ax.axhline(0, color="gray", linestyle="--", linewidth=.8)
    ax.plot(d.pseudotime, d.difference_COPD, color="black", linewidth=1.2)
    ax.fill_between(d.pseudotime, d.ci_low, d.ci_high, color="gray", alpha=.2, linewidth=0)
    ax.set(title=label(row.pathway), xlabel="Pseudotime", ylabel="COPD - Healthy score (95% CI)")
for ext in ["png", "pdf"]:
    fig.savefig(out / ("antigen_pathway_curve_differences." + ext), dpi=180, bbox_inches="tight")
plt.close(fig)
print(tests[["pathway", "FDR_shape", "FDR_any_disease"]].to_string(index=False))
print(sensitivity.to_string(index=False))

screen = pd.read_csv(out / "all_pathways_individual_mixed_models.csv").set_index("pathway")
source_screen = pd.read_csv(out / "source_curve_sensitivity.csv").set_index("pathway")
explore = pd.read_csv(out / "exploratory_pathway_fitted_curves.csv")
fig, axes = plt.subplots(2, 2, figsize=(10, 7.3), constrained_layout=True)
for ax, path in zip(axes.flat, explore.pathway.unique()):
    for status, color in colors.items():
        d = explore[explore.pathway.eq(path) & explore.status.eq(status)]
        ax.plot(d.pseudotime, d.fitted_score, color=color, linewidth=1.3, label=status)
        ax.fill_between(d.pseudotime, d.ci_low, d.ci_high, color=color, alpha=.1, linewidth=0)
    name = path.split(" R-HSA-")[0]
    if len(name) > 50:
        name = name[:47] + "..."
    ax.set(title=name + "\nShape q=" + format(screen.loc[path,"FDR_shape"],".3g") +
           "; source-curve q=" + format(source_screen.loc[path,"FDR_shape_screen_sensitivity"],".3g"),
           xlabel="Pseudotime",ylabel="Relative pathway score")
axes[0,0].legend(frameon=False)
for ext in ["png","pdf"]:
    fig.savefig(out/("exploratory_pathway_disease_curves."+ext),dpi=180,bbox_inches="tight")
plt.close(fig)
