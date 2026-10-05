"""Pick the strongest peak of each merged cluster.

Input is `bedtools merge -c 14,7 -o collapse,collapse` output: comma-joined peak
ids and signal values per cluster. Prints the id with the highest signal (the
first one on ties), one per cluster, in input order.
"""

import argparse
import sys

import polars as pl


def build_parser():
    """Build the pick_best_peak.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the pick_best_peak.py flags.
    """
    p = argparse.ArgumentParser(description="Strongest peak id per merged cluster.")
    p.add_argument(
        "merged", help="bedtools merge output with collapsed ids and signals."
    )
    return p


def main(argv=None):
    """Print the best peak id of every cluster.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
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
        .list.get(
            pl.col("signals").str.split(",").cast(pl.List(pl.Float64)).list.arg_max()
        )
    )
    best.write_csv(sys.stdout, include_header=False)


if __name__ == "__main__":
    main()
