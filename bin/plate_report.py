#!/usr/bin/env python3
"""
Generate a self-contained HTML QC report for one BRB-seq plate/library:
per-well UMI/gene counts, project breakdown, and STARsolo mapping stats
(one summary table per genome, if the plate has more than one).

Usage:
    plate_report.py --h5ad Library.h5ad --log-final Library.mouse.Log.final.out Library.human.Log.final.out \
                     --library Library --output Library_report.html
"""
import argparse
import base64
import io
import re
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np
import pandas as pd

WELL_RE = re.compile(r'^\s*([A-Ha-h])(\d{1,2})\s*$')


def parse_well(well):
    """'A01' -> (row=0, col=0); 'H12' -> (row=7, col=11). None if unparseable."""
    m = WELL_RE.match(str(well))
    if not m:
        return None
    row = ord(m.group(1).upper()) - ord('A')
    col = int(m.group(2)) - 1
    if not (0 <= row < 8 and 0 <= col < 12):
        return None
    return row, col


def plate_grid(obs, value_col):
    """8x12 float grid of value_col, keyed by obs['well']. NaN = empty/unused well."""
    grid = np.full((8, 12), np.nan)
    for _, r in obs.iterrows():
        pos = parse_well(r.get("well"))
        if pos is not None:
            grid[pos] = r[value_col]
    return grid


def render_numeric_plate(obs, value_col, title, cmap="viridis", text_grid=None, fmt="{:.0f}"):
    grid = plate_grid(obs, value_col)
    fig, ax = plt.subplots(figsize=(11, 6.5))
    im = ax.imshow(grid, cmap=cmap, aspect="equal")
    _label_plate_axes(ax)
    ax.set_title(title)
    for r in range(8):
        for c in range(12):
            if not np.isnan(grid[r, c]):
                text = text_grid[r][c] if text_grid is not None else fmt.format(grid[r, c])
                ax.text(c, r, text, ha="center", va="center", fontsize=6, color="white")
    fig.colorbar(im, ax=ax, shrink=0.8)
    return fig


def render_categorical_plate(obs, value_col, title):
    categories = sorted(obs[value_col].dropna().unique())
    cat_to_idx = {c: i for i, c in enumerate(categories)}
    grid = np.full((8, 12), np.nan)
    labels = np.full((8, 12), "", dtype=object)
    for _, r in obs.iterrows():
        pos = parse_well(r.get("well"))
        if pos is not None and pd.notna(r[value_col]):
            grid[pos] = cat_to_idx[r[value_col]]
            labels[pos] = str(r.get("sample_name", ""))

    n = max(len(categories), 1)
    cmap = ListedColormap(plt.cm.tab20.colors[:n] if n <= 20 else plt.cm.nipy_spectral(np.linspace(0, 1, n)))
    norm = BoundaryNorm(np.arange(-0.5, n, 1), cmap.N)

    fig, ax = plt.subplots(figsize=(11, 6.5))
    ax.imshow(grid, cmap=cmap, norm=norm, aspect="equal")
    _label_plate_axes(ax)
    ax.set_title(title)
    for r in range(8):
        for c in range(12):
            if labels[r, c]:
                ax.text(c, r, labels[r, c], ha="center", va="center", fontsize=5, color="black")
    legend_handles = [Patch(facecolor=cmap(i), label=str(cat)) for cat, i in cat_to_idx.items()]
    ax.legend(handles=legend_handles, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8, title=value_col)
    return fig


def _label_plate_axes(ax):
    ax.set_xticks(range(12))
    ax.set_xticklabels([f"{i+1:02d}" for i in range(12)])
    ax.set_yticks(range(8))
    ax.set_yticklabels(list("ABCDEFGH"))
    ax.set_xticks(np.arange(-0.5, 12, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, 8, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="minor", length=0)


def fig_to_base64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=110)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("utf-8")


def parse_log_final(path):
    stats = {}
    with open(path) as fh:
        for line in fh:
            if "|" in line:
                key, val = line.split("|", 1)
                stats[key.strip()] = val.strip()
    return stats


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--log-final", required=True, nargs="+",
                     help="STARsolo *.Log.final.out, one or more (one per genome on this plate)")
    ap.add_argument("--well-totals", required=True, nargs="+",
                     help="well_totals.csv from mtx_to_h5ad.py, one per genome on this plate -- "
                          "carries every well's UMI count from every genome's alignment, needed to spot mixups")
    ap.add_argument("--library", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    adata = ad.read_h5ad(args.h5ad)
    adata.obs["n_counts"] = adata.X.sum(axis=1).A1 if hasattr(adata.X, "A1") else adata.X.sum(axis=1)
    adata.obs["n_genes"] = (adata.X > 0).sum(axis=1).A1 if hasattr(adata.X, "A1") else (adata.X > 0).sum(axis=1)

    well_cols = ["barcode", "well", "sample_name", "project", "genome", "n_counts", "n_genes"]
    well_cols = [c for c in well_cols if c in adata.obs.columns]
    per_well = adata.obs[well_cols].sort_values("n_counts", ascending=False).reset_index(drop=True)

    # --- cross-genome well totals: every well's UMI count from EVERY ---
    # --- genome's alignment, not just its assigned one -- this is what ---
    # --- catches mixups that the (already genome-filtered) h5ad above can't ---
    totals = pd.concat([pd.read_csv(f) for f in args.well_totals], ignore_index=True)
    genomes_present = sorted(totals["this_genome"].unique())
    multi_genome = len(genomes_present) > 1

    wide = totals.pivot_table(index=["barcode", "well", "sample_name", "project", "assigned_genome"],
                               columns="this_genome", values="umis_this_genome", fill_value=0).reset_index()
    genome_cols = [g for g in genomes_present if g in wide.columns]

    def assigned_umis(row):
        g = row.get("assigned_genome")
        return row[g] if g in genome_cols and pd.notna(g) else float("nan")

    def off_target_umis(row):
        g = row.get("assigned_genome")
        others = [row[gc] for gc in genome_cols if gc != g]
        return sum(others) if others else 0.0

    wide["assigned_umis"] = wide.apply(assigned_umis, axis=1)
    wide["off_target_umis"] = wide.apply(off_target_umis, axis=1)
    total_all = wide["assigned_umis"].fillna(0) + wide["off_target_umis"]
    wide["total_umis"] = total_all
    wide["off_target_frac"] = np.where(total_all > 0, wide["off_target_umis"] / total_all, 0.0)

    MIXUP_THRESHOLD = 0.10   # flag wells where >10% of UMIs came from an unassigned genome
    MIN_UMIS_TO_FLAG = 500   # ...but only if there's enough total signal for that ratio to mean anything --
                              # a near-empty well (e.g. 22 vs 33 UMIs) can show a wild ratio from pure noise
    eligible = wide["total_umis"] >= MIN_UMIS_TO_FLAG
    mixups = wide[eligible & (wide["off_target_frac"] > MIXUP_THRESHOLD)].sort_values("off_target_frac", ascending=False)
    low_depth_skipped = int((~eligible & (wide["off_target_frac"] > MIXUP_THRESHOLD)).sum())

    if "well" not in adata.obs.columns or adata.obs["well"].isna().all():
        plate_map_html = "<p><em>No 'well' column in the plate map -- can't render a spatial layout. " \
                          "Add well positions (A01-H12) to the plate map to enable this.</em></p>"
    else:
        sections = []

        if multi_genome:
            # per-well text showing every genome's UMI count, so both
            # numbers are readable directly off the main heatmap
            text_grid = [["" for _ in range(12)] for _ in range(8)]
            for _, r in wide.iterrows():
                pos = parse_well(r["well"])
                if pos is None:
                    continue
                lines = [f"{g[:3]}:{int(r[g])}" for g in genome_cols]
                text_grid[pos[0]][pos[1]] = "\n".join(lines)
            fig_umi_plate = render_numeric_plate(adata.obs, "n_counts",
                                                  f"{args.library}: UMI count by well position (assigned genome)",
                                                  text_grid=text_grid)
        else:
            fig_umi_plate = render_numeric_plate(adata.obs, "n_counts", f"{args.library}: UMI count by well position")
        sections.append(f'<h3>UMI count</h3><img src="data:image/png;base64,{fig_to_base64(fig_umi_plate)}">')

        fig_genes_plate = render_numeric_plate(adata.obs, "n_genes", f"{args.library}: genes detected by well position", cmap="plasma")
        sections.append(f'<h3>Genes detected</h3><img src="data:image/png;base64,{fig_to_base64(fig_genes_plate)}">')

        if "project" in adata.obs.columns and adata.obs["project"].nunique(dropna=True) > 0:
            fig_project_plate = render_categorical_plate(adata.obs, "project", f"{args.library}: project by well position")
            sections.append(f'<h3>Project layout</h3><img src="data:image/png;base64,{fig_to_base64(fig_project_plate)}">')

        if "genome" in adata.obs.columns and adata.obs["genome"].nunique(dropna=True) > 1:
            fig_genome_plate = render_categorical_plate(adata.obs, "genome", f"{args.library}: genome by well position")
            sections.append(f'<h3>Genome layout (mixed-genome plate)</h3><img src="data:image/png;base64,{fig_to_base64(fig_genome_plate)}">')

        if multi_genome:
            wide_indexed = wide.set_index("well")
            obs_for_mixup = adata.obs.copy()
            obs_for_mixup["off_target_frac"] = obs_for_mixup["well"].map(wide_indexed["off_target_frac"]) * 100
            fig_mixup_plate = render_numeric_plate(
                obs_for_mixup, "off_target_frac",
                f"{args.library}: % UMIs from an UNASSIGNED genome (mixup indicator)",
                cmap="Reds", fmt="{:.1f}%",
            )
            sections.append(f'<h3>Possible mixups (&gt;{int(MIXUP_THRESHOLD*100)}% off-target flagged below)</h3>'
                             f'<img src="data:image/png;base64,{fig_to_base64(fig_mixup_plate)}">')

        plate_map_html = "".join(sections)

    if multi_genome and len(mixups):
        mixup_cols = ["well", "sample_name", "project", "assigned_genome"] + genome_cols + ["off_target_frac"]
        mixup_display = mixups[mixup_cols].copy()
        mixup_display["off_target_frac"] = (mixup_display["off_target_frac"] * 100).round(1).astype(str) + "%"
        mixup_html = mixup_display.to_html(classes="tbl", border=0, index=False)
        if low_depth_skipped:
            mixup_html += (f"<p><em>{low_depth_skipped} well(s) also showed &gt;{int(MIXUP_THRESHOLD*100)}% "
                            f"off-target signal but were skipped: fewer than {MIN_UMIS_TO_FLAG} total UMIs, "
                            f"where the ratio isn't statistically meaningful (these are more likely just "
                            f"low-depth/failed wells -- check the UMI count plate map for them).</em></p>")
    elif multi_genome:
        mixup_html = f"<p>No wells exceeded the {int(MIXUP_THRESHOLD*100)}% off-target-genome threshold " \
                      f"(among wells with at least {MIN_UMIS_TO_FLAG} total UMIs).</p>"
        if low_depth_skipped:
            mixup_html += (f"<p><em>{low_depth_skipped} low-depth well(s) showed an elevated ratio but were "
                            f"skipped as statistically unreliable below {MIN_UMIS_TO_FLAG} total UMIs.</em></p>")
    else:
        mixup_html = "<p><em>Only one genome on this plate -- nothing to cross-check.</em></p>"

    # UMI count per well, bar chart (ranked -- complements the spatial view above)
    fig1, ax1 = plt.subplots(figsize=(10, 4))
    ax1.bar(range(len(per_well)), per_well["n_counts"], color="#4c72b0")
    ax1.set_xlabel("Well (ranked by UMI count)")
    ax1.set_ylabel("UMI count")
    ax1.set_title(f"{args.library}: UMI counts per well")
    umi_plot = fig_to_base64(fig1)

    # genes detected per well
    fig2, ax2 = plt.subplots(figsize=(10, 4))
    ax2.bar(range(len(per_well)), per_well["n_genes"], color="#55a868")
    ax2.set_xlabel("Well (ranked by UMI count)")
    ax2.set_ylabel("Genes detected")
    ax2.set_title(f"{args.library}: genes detected per well")
    genes_plot = fig_to_base64(fig2)

    star_sections = []
    for log_path in args.log_final:
        # filenames are {library}.{genome}.Log.final.out
        base = Path(log_path).name
        genome_label = base[len(args.library) + 1:].removesuffix(".Log.final.out") if base.startswith(args.library + ".") else base
        stats = parse_log_final(log_path)
        rows = "".join(f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in stats.items())
        star_sections.append(f"<h3>Genome: {genome_label}</h3><table class=\"tbl\">{rows}</table>")
    star_rows = "".join(star_sections)

    project_counts = per_well.groupby("project", dropna=False)["n_counts"].agg(["count", "sum", "mean"])
    project_table = project_counts.to_html(classes="tbl", border=0)

    well_table = per_well.to_html(classes="tbl", border=0, index=False)

    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>{args.library} -- Plate QC Report</title>
<style>
body {{ font-family: -apple-system, sans-serif; margin: 2em; color: #222; }}
h1, h2 {{ color: #2c3e50; }}
table.tbl {{ border-collapse: collapse; margin-bottom: 2em; }}
table.tbl td, table.tbl th {{ padding: 4px 10px; border-bottom: 1px solid #ddd; text-align: left; }}
img {{ max-width: 100%; margin-bottom: 2em; }}
</style></head><body>
<h1>Plate QC Report: {args.library}</h1>

<h2>Plate map</h2>
{plate_map_html}

<h2>Possible sample mixups</h2>
{mixup_html}

<h2>STARsolo mapping summary</h2>
{star_rows}

<h2>UMI counts per well</h2>
<img src="data:image/png;base64,{umi_plot}">

<h2>Genes detected per well</h2>
<img src="data:image/png;base64,{genes_plot}">

<h2>Per-project summary (wells on this plate)</h2>
{project_table}

<h2>Per-well detail</h2>
{well_table}

</body></html>"""

    with open(args.output, "w") as fh:
        fh.write(html)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
