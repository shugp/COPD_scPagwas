"""Verify sample labels and biological donor units from the analysis metadata."""
import argparse
import json
from pathlib import Path
import pandas as pd


def audit(metadata):
    required = ["data_source", "orig.ident", "status"]
    # Some exports keep donor IDs but omit the original sample labels.
    if "orig.ident" not in metadata and "myeloid_donor" in metadata:
        raise ValueError("orig.ident is absent; use the original/retained object metadata, not donor IDs alone")
    samples = metadata[required].drop_duplicates().copy()
    if samples.duplicated(["data_source", "orig.ident"]).any():
        raise ValueError("Conflicting disease labels for the same source/sample")
    samples["donor"] = samples["orig.ident"].astype(str)
    repeated = samples["data_source"].eq("GSE168191")
    samples.loc[repeated, "donor"] = samples.loc[repeated, "donor"].str.replace(
        r"_rep[0-9]+\.h5$", "", regex=True)
    samples["donor"] = samples["data_source"] + ":" + samples["donor"]
    donors = samples[["data_source", "donor", "status"]].drop_duplicates()
    if donors.duplicated("donor").any():
        raise ValueError("Conflicting disease labels for the same biological donor")
    summary = {"sample_labels": len(samples), "donor_units": len(donors),
               "samples_by_group": samples["status"].value_counts().to_dict(),
               "donors_by_group": donors["status"].value_counts().to_dict(),
               "donors_by_source": donors.groupby("data_source").size().to_dict(),
               "basis": "GSE168191 rep1/rep2 labels are merged within donor; other source-qualified orig.ident labels are donor units."}
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--metadata", type=Path, required=True)
    args = p.parse_args()
    print(json.dumps(audit(pd.read_csv(args.metadata)), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
