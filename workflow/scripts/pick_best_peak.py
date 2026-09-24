"""Pick the strongest peak of each merged cluster.

Input is `bedtools merge -c 14,7 -o collapse,collapse` output: comma-joined peak
ids and signal values per cluster. Prints the id with the highest signal (the
first one on ties), one per cluster, in input order.
"""

import argparse
import sys

import polars as pl


def build_parser():
    """Return the argument parser for pick_best_peak.py."""
    p = argparse.ArgumentParser(description="Strongest peak id per merged cluster.")
    p.add_argument("merged", help="bedtools merge output with collapsed ids and signals.")
    return p


def main(argv=None):
    """Print the best peak id of every cluster."""
    args = build_parser().parse_args(argv)
    clusters = pl.read_csv(
        args.merged,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "start", "end", "ids", "signals"],
        schema_overrides={"ids": pl.Utf8, "signals": pl.Utf8},
    )
    best = clusters.select(
        pl.col("ids")
        .str.split(",")
        .list.get(pl.col("signals").str.split(",").cast(pl.List(pl.Float64)).list.arg_max())
    )
    best.write_csv(sys.stdout, include_header=False)


if __name__ == "__main__":
    main()
