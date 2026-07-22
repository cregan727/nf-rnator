#!/usr/bin/env python3
"""
Generate a self-contained HTML QC report for one project, aggregating
across every plate/library that contributed samples to it.

Usage:
    project_report.py --h5ad ProjectX.h5ad --project ProjectX --output ProjectX_report.html
"""
import argparse
import base64
import io

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def fig_to_base64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=110)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("utf-8")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5ad", required=True, help="Merged per-project .h5ad from merge_project_h5ad.py")
    ap.add_argument("--project", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    adata = ad.read_h5ad(args.h5ad)
    adata.obs["n_counts"] = adata.X.sum(axis=1).A1 if hasattr(adata.X, "A1") else adata.X.sum(axis=1)
    adata.obs["n_genes"] = (adata.X > 0).sum(axis=1).A1 if hasattr(adata.X, "A1") else (adata.X > 0).sum(axis=1)

    per_sample = adata.obs[["barcode", "sample_name", "library", "n_counts", "n_genes"]] \
        .sort_values("n_counts", ascending=False).reset_index(drop=True)

    fig1, ax1 = plt.subplots(figsize=(10, 4))
    ax1.bar(range(len(per_sample)), per_sample["n_counts"], color="#4c72b0")
    ax1.set_xlabel("Sample (ranked by UMI count)")
    ax1.set_ylabel("UMI count")
    ax1.set_title(f"{args.project}: UMI counts per sample (all plates)")
    umi_plot = fig_to_base64(fig1)

    n_libraries = adata.obs["library"].nunique()
    per_library = adata.obs.groupby("library", dropna=False)["n_counts"].agg(["count", "sum", "mean"])
    library_table = per_library.to_html(classes="tbl", border=0)
    sample_table = per_sample.to_html(classes="tbl", border=0, index=False)

    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>{args.project} -- Project QC Report</title>
<style>
body {{ font-family: -apple-system, sans-serif; margin: 2em; color: #222; }}
h1, h2 {{ color: #2c3e50; }}
table.tbl {{ border-collapse: collapse; margin-bottom: 2em; }}
table.tbl td, table.tbl th {{ padding: 4px 10px; border-bottom: 1px solid #ddd; text-align: left; }}
img {{ max-width: 100%; margin-bottom: 2em; }}
</style></head><body>
<h1>Project QC Report: {args.project}</h1>
<p>{adata.n_obs} samples across {n_libraries} plate(s)/library(ies).</p>

<h2>UMI counts per sample</h2>
<img src="data:image/png;base64,{umi_plot}">

<h2>Per-plate summary</h2>
{library_table}

<h2>Per-sample detail</h2>
{sample_table}

</body></html>"""

    with open(args.output, "w") as fh:
        fh.write(html)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
