"""Keep rPeaks whose summit overlaps come from at least N experiments.

Input is `bedtools intersect -wo -a rpeaks -b summits`: 14 rPeak columns, 14
summit columns, overlap. Per rPeak (in input order) prints chrom, start, end,
id, then the sorted ';'-joined targets, biosamples and experiments of its
overlapping summits, and the experiment count.
"""

import argparse
import sys

import polars as pl

# 0-based columns of the intersection.
COLUMNS = {0: "chrom", 1: "start", 2: "end", 13: "id", 24: "exp", 25: "target", 26: "biosample"}


def build_parser():
    """Return the argument parser for filter_rpeaks.py."""
    p = argparse.ArgumentParser(description="Filter rPeaks by supporting experiments.")
    p.add_argument("intersection", help="bedtools intersect -wo of rPeaks and summits.")
    p.add_argument("--min-experiments", type=int, default=5, help="Minimum experiments.")
    return p


def joined(col):
    """Sorted, deduplicated items of a column, ';'-joined."""
    return pl.col(col).unique().sort().str.join(";")


def main(argv=None):
    """Aggregate summits per rPeak and print the supported ones."""
    args = build_parser().parse_args(argv)
    overlaps = pl.read_csv(
        args.intersection,
        separator="\t",
        has_header=False,
        columns=list(COLUMNS),
        new_columns=list(COLUMNS.values()),
        infer_schema_length=0,
        quote_char=None,
    )
    (
        overlaps.group_by("id", maintain_order=True)
        .agg(
            pl.col("chrom", "start", "end").first(),
            joined("target"),
            joined("biosample"),
            joined("exp"),
            n=pl.col("exp").n_unique(),
        )
        .filter(pl.col("n") >= args.min_experiments)
        .select("chrom", "start", "end", "id", "target", "biosample", "exp", "n")
        .write_csv(sys.stdout, separator="\t", include_header=False, quote_style="never")
    )


if __name__ == "__main__":
    main()
