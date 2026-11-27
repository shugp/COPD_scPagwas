#!/usr/bin/env python3
"""Run the Myeloid analysis or regenerate the four manuscript deliverables."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

REPO = Path(__file__).resolve().parent
STAGES = {
    "scpagwas": ["inte_all.rds", "copdgwas.txt"],
    "annotation": ["Myeloid/Myeloid.Rds", "coloc.CSV", "pagwas/inte_all_gene_PCC.csv"],
    "trajectory": ["Myeloid/results_cleaned/Myeloid.cleaned.rds",
                   "Myeloid/results_cleaned/coloc_overall_gene_tests.csv"],
    "disease": ["Myeloid/results_cleaned/Myeloid.cleaned.rds",
                "Myeloid/trajectory/cell_pseudotime_state.csv",
                "Myeloid/results_cleaned/coloc_overall_gene_tests.csv",
                "Myeloid/disease_pathways/MSigDB_Hallmark_2020.gmt",
                "Myeloid/disease_pathways/Reactome_2022.gmt"],
    "documents": ["coloc.CSV"],
    "audit": ["Myeloid/Myeloid.Rds"],
}
STAGES["figures"] = [
    "gPAS/scPagwas.inte_all.bar_positie_nagtive1.pdf",
    "Myeloid/trajectory/cell_pseudotime_state.csv",
    "Myeloid/results_cleaned/umap_coordinates.csv",
    "Myeloid/results_cleaned/annotation_marker_expression.csv",
    "Myeloid/results_cleaned/coloc_overall_gene_tests.csv",
    "Myeloid/disease_pathways/gsea_all_results.csv",
    "Myeloid/disease_pathways/antigen_residual_coexpression.csv",
    "Myeloid/disease_pathways/candidate_module_reactome_enrichment.csv",
    "Myeloid/disease_pathways/Reactome_2022.gmt",
    "Myeloid/disease_pathways/pathway_pseudotime/antigen_pathway_fitted_curves.csv",
    "Myeloid/disease_pathways/pathway_pseudotime/antigen_pathway_disease_pseudotime_tests.csv",
    "Myeloid/disease_pathways/pathway_pseudotime/antigen_source_curve_sensitivity.csv",
]
STAGES["render"] = list(dict.fromkeys(STAGES["figures"] + STAGES["documents"]))
STAGES["analysis"] = list(dict.fromkeys(STAGES["annotation"] + [
    "Myeloid/disease_pathways/MSigDB_Hallmark_2020.gmt",
    "Myeloid/disease_pathways/Reactome_2022.gmt",
]))
STAGES["diagnostics"] = ["Myeloid/Myeloid.Rds"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True,
                        help="Directory containing Myeloid/, gPAS/, coloc.CSV, etc.")
    parser.add_argument("--stage", choices=STAGES, default="render")
    parser.add_argument("--output-dir", type=Path,
                        help="Destination for render/documents/figures; default DATA_ROOT/Myeloid")
    parser.add_argument("--check", action="store_true", help="Check inputs without running the stage")
    parser.add_argument("--repository-url", help="Real public repository URL for the manuscript declaration")
    args = parser.parse_args()
    data_root = args.data_root.expanduser().resolve()
    required = STAGES[args.stage]
    missing = [name for name in required if not (data_root / name).is_file()]
    if missing:
        parser.exit(2, "Missing inputs for " + args.stage + ":\n" +
                    "\n".join("  " + str(data_root / name) for name in missing) +
                    "\nSee docs/INPUTS.md. No analysis has been run.\n")
    if args.check:
        print(json.dumps({"stage": args.stage, "data_root": str(data_root),
                          "required_inputs": required, "status": "ready"}, indent=2))
        return
    env = os.environ.copy()
    env["MYELOID_CODE_DIR"] = str(REPO / "analysis")
    # Nested R system2('python', ...) uses the interpreter that launched this runner.
    env["PATH"] = str(Path(sys.executable).parent) + os.pathsep + env.get("PATH", "")
    out = (args.output_dir or data_root / "Myeloid").expanduser().resolve()

    def run(command):
        print("Running:", " ".join(map(str, command)), flush=True)
        subprocess.run(list(map(str, command)), cwd=data_root, env=env, check=True)

    if args.stage in ("render", "figures"):
        if not shutil.which("pdftocairo"):
            parser.exit(2, "pdftocairo is required; install Poppler (see README.md).\n")
        run([sys.executable, REPO / "tools/render_figure.py", "--data-root", data_root,
             "--output-dir", out])
    if args.stage in ("render", "documents"):
        command = [sys.executable, REPO / "tools/write_documents.py", "--data-root", data_root,
                   "--output-dir", out]
        if args.repository_url:
            if not args.repository_url.startswith("https://"):
                parser.exit(2, "--repository-url must be an actual public HTTPS URL.\n")
            command.extend(["--repository-url", args.repository_url])
        run(command)
    if args.stage == "audit":
        if not shutil.which("Rscript"):
            parser.exit(2, "The donor audit reads the original Seurat object and requires Rscript.\n")
        run(["Rscript", REPO / "tools/audit_donors.R"])
    scripts = {
        "scpagwas": ["01_scPagwas_fast.R"],
        "annotation": ["02_Myeloid.R"],
        "trajectory": ["03_Myeloid_trajectory.R"],
        "disease": ["04_Myeloid_disease_pathways.R"],
        "analysis": ["02_Myeloid.R", "03_Myeloid_trajectory.R", "04_Myeloid_disease_pathways.R"],
        "diagnostics": ["00_resolution_diagnostics.R"],
    }
    if args.stage in scripts:
        if args.output_dir:
            parser.exit(2, "Analysis stages write results under --data-root; --output-dir is for rendering only.\n")
        if not shutil.which("Rscript"):
            parser.exit(2, "Rscript was not found. Install the R dependencies described in README.md.\n")
        for script in scripts[args.stage]:
            run(["Rscript", REPO / "analysis" / script])


if __name__ == "__main__":
    main()
