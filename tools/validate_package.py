"""Check syntax, resource consistency, and portable entry points without analysis data."""
import ast
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

REPO = Path(__file__).resolve().parents[1]


def main():
    python_files = sorted(REPO.rglob("*.py"))
    for path in python_files:
        ast.parse(path.read_text(), filename=str(path))
    for path in REPO.rglob("*.json"):
        json.loads(path.read_text())
    manuscript = json.loads((REPO / "resources/manuscript_cn.json").read_text())
    references = json.loads((REPO / "resources/references.json").read_text())
    assert len(manuscript["body"]) == 5, "The Research Letter must have five body paragraphs"
    assert 5 <= len(references) <= 10
    assert "{repository_url}" in manuscript["availability"]
    assert len({r["doi"] for r in references}) == len(references)
    assert references[0]["doi"] == "10.1016/j.xgen.2023.100383"
    source_dois = {"10.1038/s41586-022-04552-0", "10.1186/s12931-022-02293-2",
                   "10.1165/rcmb.2021-0555OC", "10.1016/j.immuni.2023.01.032"}
    assert source_dois.issubset({r["doi"] for r in references})
    cites = set()
    for paragraph in manuscript["body"]:
        for group in re.findall(r"\[([0-9,–—-]+)\]", paragraph):
            for item in group.split(","):
                ends = re.split("[–—-]", item)
                cites.update(range(int(ends[0]), int(ends[-1]) + 1))
    assert cites == set(range(1, len(references) + 1)), "Missing or invalid manuscript citation"
    assert "PP.H4 > 0.6" in " ".join(manuscript["body"])
    for path in (REPO / "analysis").glob("*"):
        if path.suffix in (".R", ".py"):
            text = path.read_text()
            assert "/data/sgp/" not in text
            assert not re.search(r"/tmp/.*\.(rds|csv|pkl)", text, re.I), path
    manifest = json.loads((REPO / "resources/source_code_manifest.json").read_text())
    for entry in manifest:
        path = REPO / entry["packaged"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == entry["packaged_sha256"], path
    rscript = shutil.which("Rscript")
    r_files = sorted(REPO.rglob("*.R"))
    if rscript:
        # Parse only; never source analysis scripts or create analysis outputs.
        expr = 'args <- commandArgs(TRUE); for (f in args) invisible(parse(file=f)); cat(length(args), "R files parsed\\n")'
        subprocess.run([rscript, "-e", expr, *map(str, r_files)], check=True)
    print(f"PASS: {len(python_files)} Python files, {len(r_files)} R files, five paragraphs and {len(references)} references.")
    print("R syntax parsed." if rscript else "Rscript absent: R syntax check skipped (run locally with R installed).")


if __name__ == "__main__":
    main()
