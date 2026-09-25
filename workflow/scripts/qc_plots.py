"""QC figures for the new ATAC anchors and cCREs.

Each page is saved as SVG and PNG under --figures, and all pages go into one
PDF. Page 3 (ATAC vs DNase) is drawn only when the support table has DNase.
"""

import argparse
from pathlib import Path

import numpy as np
import polars as pl
from _style import apply_style, despine, save_figure

apply_style()

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

THRESHOLD = 1.64
ZERO_Z = -10.0  # zscore.py's value for anchors with no signal
TOP_SAMPLES = 40
WINDOW = 10_000


def build_parser():
    """Return the argument parser for qc_plots.py."""
    p = argparse.ArgumentParser(description="QC figures for the new ATAC anchors and cCREs.")
    p.add_argument("--support", required=True, help="ccre_support.py table.")
    p.add_argument("--atac-max-z", required=True, help="ATAC anchor/max_zscore TSV.")
    p.add_argument("--peaks", required=True, help="Pooled peaks BED.")
    p.add_argument("--rpeaks", required=True, help="rPeaks BED.")
    p.add_argument("--filtered", required=True, help="rPeaks with enough experiments.")
    p.add_argument("--no-rdhs", required=True, help="Filtered rPeaks off the rDHSs.")
    p.add_argument("--figures", required=True, help="Directory for the per-page SVG and PNG.")
    p.add_argument("-o", "--output", required=True, help="Multi-page PDF to write.")
    return p


def count_lines(path):
    """Number of lines in a file."""
    with open(path, "rb") as f:
        return sum(chunk.count(b"\n") for chunk in iter(lambda: f.read(1 << 24), b""))


def funnel(args, support):
    """Regions left after each filtering step, from pooled peaks to ATAC cCREs."""
    steps = {
        "Peaks": count_lines(args.peaks),
        "rPeaks": count_lines(args.rpeaks),
        "≥ 5\nexperiments": count_lines(args.filtered),
        "No rDHS\noverlap": count_lines(args.no_rdhs),
        "Not\nblacklisted": support.height,
        "ATAC\ncCREs": support["is_ccre"].sum(),
    }
    fig, ax = plt.subplots(figsize=(6, 3.5))
    bars = ax.bar(list(steps), list(steps.values()), color="C0")
    ax.bar_label(bars, labels=[f"{n:,}" for n in steps.values()], padding=2)
    ax.set_yscale("log")
    ax.set_ylabel("Regions")
    ax.set_title("Filtering funnel")
    despine(ax, categorical_x=True)
    return fig


def atac_max_z(args, support):
    """ATAC max z-score of registry and new anchors."""
    max_z = pl.read_csv(args.atac_max_z, separator="\t", schema_overrides={"anchor": pl.Utf8})
    is_new = max_z["anchor"].is_in(support["anchor"].implode()).to_numpy()
    z = max_z["max_zscore"].to_numpy()
    bins = np.linspace(np.nanmin(z[z > ZERO_Z]), np.nanmax(z), 100)
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.hist(z[~is_new], bins=bins, density=True, histtype="step", label=f"Registry ({(~is_new).sum():,})")
    ax.hist(z[is_new], bins=bins, density=True, histtype="step", label=f"New ATAC ({is_new.sum():,})")
    ax.axvline(THRESHOLD, color="black", linestyle="--", linewidth=0.8)
    ax.set_xlabel("ATAC max z-score")
    ax.set_ylabel("Density")
    ax.set_title("ATAC max z-score by anchor origin")
    ax.legend(frameon=False)
    despine(ax)
    return fig


def atac_vs_dnase(support):
    """ATAC vs DNase max z-score of the new anchors, with quadrant counts."""
    atac = support["atac_max_z"].to_numpy()
    dnase = support["dnase_max_z"].to_numpy()
    shown = dnase > ZERO_Z
    fig, ax = plt.subplots(figsize=(5, 4))
    hb = ax.hexbin(atac[shown], dnase[shown], gridsize=60, bins="log", mincnt=1, cmap="viridis")
    fig.colorbar(hb, ax=ax, label="Anchors")
    ax.axvline(THRESHOLD, color="black", linestyle="--", linewidth=0.8)
    ax.axhline(THRESHOLD, color="black", linestyle="--", linewidth=0.8)
    for x_high, y_high, x, y in [(0, 1, 0.02, 0.98), (1, 1, 0.98, 0.98), (0, 0, 0.02, 0.02), (1, 0, 0.98, 0.02)]:
        n = ((atac > THRESHOLD) == x_high) & ((dnase > THRESHOLD) == y_high)
        ax.text(x, y, f"{n.sum():,}", transform=ax.transAxes, ha="right" if x_high else "left",
                va="top" if y_high else "bottom")
    ax.set_xlabel("ATAC max z-score")
    ax.set_ylabel("DNase max z-score")
    ax.set_title(f"New anchors; {(~shown).sum():,} with no DNase signal not shown")
    despine(ax)
    return fig


def support_page(support):
    """Per-sample z-score support and peak support of the new anchors."""
    ccre = support["is_ccre"].to_numpy()
    passing = support["samples_passing"].to_numpy()
    experiments = support["peak_experiments"].to_numpy()
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(11, 3.5))

    bins = np.arange(0, np.percentile(passing, 99) + 2)
    ax1.hist([passing[ccre], passing[~ccre]], bins=bins, stacked=True, label=["ATAC cCRE", "Not called"])
    ax1.set_xlabel(f"Samples with z ≥ {THRESHOLD}")
    ax1.set_ylabel("New anchors")
    ax1.legend(frameon=False)

    bins = np.arange(experiments.min(), np.percentile(experiments, 99) + 2)
    ax2.hist([experiments[ccre], experiments[~ccre]], bins=bins, stacked=True)
    ax2.set_xlabel("Experiments supporting the rPeak")
    ax2.set_ylabel("New anchors")

    hb = ax3.hexbin(experiments, passing, gridsize=50, bins="log", mincnt=1, cmap="viridis")
    fig.colorbar(hb, ax=ax3, label="Anchors")
    ax3.set_xlabel("Experiments supporting the rPeak")
    ax3.set_ylabel(f"Samples with z ≥ {THRESHOLD}")
    for ax in (ax1, ax2, ax3):
        despine(ax)
    fig.suptitle("Support per new anchor (histograms clipped at the 99th percentile)")
    return fig


def nearest_registry(support):
    """Distance from each ATAC cCRE midpoint to the nearest registry anchor midpoint."""
    ccres = support.filter("is_ccre")
    left = -ccres["left_registry_bp"].drop_nulls().to_numpy()
    right = ccres["right_registry_bp"].drop_nulls().to_numpy()
    bins = np.arange(-WINDOW, WINDOW + 100, 100)
    fig, ax = plt.subplots(figsize=(8, 3.5))
    ax.hist(left[left >= -WINDOW], bins=bins, alpha=0.7, label="Left")
    ax.hist(right[right <= WINDOW], bins=bins, alpha=0.7, label="Right")
    for x in (-175, 175):
        ax.axvline(x, color="black", linestyle=":", linewidth=0.6)
    for x in (-75, 75):
        ax.axvline(x, color="black", linestyle="--", linewidth=0.6)
    ax.set_xlabel("Distance to nearest registry anchor (bp)")
    ax.set_ylabel("ATAC cCREs")
    ax.set_title("Nearest registry anchor around ATAC cCREs")
    ax.legend(frameon=False)
    despine(ax)
    return fig


def sample_contribution(support):
    """ATAC cCREs whose rPeak came from each sample, top samples only."""
    counts = (
        support.filter("is_ccre")
        .group_by("rpeak_sample")
        .len()
        .sort(["len", "rpeak_sample"], descending=[True, False])
    )
    top = counts.head(TOP_SAMPLES)
    colors = ["C0" if s.startswith("ENCSR") else "C1" for s in top["rpeak_sample"]]
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.bar(top["rpeak_sample"], top["len"], color=colors)
    ax.tick_params(axis="x", labelrotation=90)
    ax.set_ylabel("ATAC cCREs")
    ax.set_title(f"rPeak source, top {len(top)} of {counts.height} samples")
    ax.legend(handles=[Patch(color="C0", label="ENCODE"), Patch(color="C1", label="Custom")], frameon=False)
    despine(ax, categorical_x=True)
    return fig


def main(argv=None):
    """Draw every page, save each as SVG and PNG, and collect them into one PDF."""
    args = build_parser().parse_args(argv)
    support = pl.read_csv(args.support, separator="\t")

    pages = {
        "1-funnel": lambda: funnel(args, support),
        "2-atac-max-z": lambda: atac_max_z(args, support),
        "3-atac-vs-dnase": lambda: atac_vs_dnase(support),
        "4-support": lambda: support_page(support),
        "5-nearest-registry": lambda: nearest_registry(support),
        "6-sample-contribution": lambda: sample_contribution(support),
    }
    if "dnase_max_z" not in support.columns:
        del pages["3-atac-vs-dnase"]

    with PdfPages(args.output) as pdf:
        for name, draw in pages.items():
            fig = draw()
            save_figure(fig, Path(args.figures) / name)
            pdf.savefig(fig)
            plt.close(fig)


if __name__ == "__main__":
    main()
