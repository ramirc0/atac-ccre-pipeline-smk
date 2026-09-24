"""Per-anchor maximum z-score across z-score files.

Every file must list the same anchors in the same order. Missing (NaN)
z-scores never win; an anchor missing in every file gets an empty max.
"""

import argparse

import numpy as np
import polars as pl


def build_parser():
    """Return the argument parser for max_z.py."""
    p = argparse.ArgumentParser(description="Per-anchor max z-score across files.")
    p.add_argument("-i", "--inputs", nargs="+", required=True, help="zscore.py parquet files.")
    p.add_argument("-o", "--output", required=True, help="anchor/max_zscore TSV to write.")
    return p


def read_zscores(path):
    """Read the anchor and zscore columns of a zscore.py parquet file."""
    return pl.read_parquet(path, columns=["anchor", "zscore"])


def main(argv=None):
    """Fold every file into a running max and write it."""
    args = build_parser().parse_args(argv)

    first = read_zscores(args.inputs[0])
    anchor = first["anchor"]
    max_z = first["zscore"].to_numpy(writable=True)

    for i, path in enumerate(args.inputs[1:], start=2):
        print(f"Processing {i}/{len(args.inputs)}: {path}")
        df = read_zscores(path)
        if len(df) != len(anchor):
            raise ValueError(f"Row count mismatch in {path}: expected {len(anchor)}, got {len(df)}")
        if not df["anchor"].equals(anchor):
            raise ValueError(f"Anchor order mismatch in {path}")

        z = df["zscore"].to_numpy()
        better = np.nan_to_num(z, nan=-np.inf) > np.nan_to_num(max_z, nan=-np.inf)
        max_z[better] = z[better]

    pl.DataFrame({"anchor": anchor, "max_zscore": max_z}).with_columns(
        pl.col("max_zscore").fill_nan(None)
    ).write_csv(args.output, separator="\t")


if __name__ == "__main__":
    main()
