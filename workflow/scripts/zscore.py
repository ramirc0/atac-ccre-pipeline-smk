"""Z-score of log10 mean bigWig signal per anchor.

Reads bigWigAverageOverBed output (name, size, covered, sum, mean0, mean). The
mean and standard deviation come from anchors with non-zero signal only;
anchors with zero signal or zero coverage get z = -10. Writes anchor, zscore,
signal and the descending min-rank of zscore, sorted by anchor.
"""

import argparse
import math

import numpy as np
import polars as pl

ZERO_Z = -10.0


def build_parser():
    """Return the argument parser for zscore.py."""
    p = argparse.ArgumentParser(description="Z-score log10 bigWig signal per anchor.")
    p.add_argument("-i", "--input", required=True, help="bigWigAverageOverBed tab output.")
    p.add_argument("-o", "--output", required=True, help="Z-score parquet to write.")
    return p


def zscores(signal, nonzero):
    """Z-score of log10(signal) over `nonzero` anchors; ZERO_Z elsewhere."""
    # math.log(x, 10), not np.log10: matches the original values bit for bit.
    log_signal = np.fromiter((math.log(x, 10) for x in signal[nonzero]), float)
    z = np.full(len(signal), ZERO_Z)
    z[nonzero] = (log_signal - np.mean(log_signal)) / np.std(log_signal)
    return z


def main(argv=None):
    """Score one bigWig's anchor signal and write it sorted by anchor."""
    args = build_parser().parse_args(argv)
    df = pl.read_csv(
        args.input,
        separator="\t",
        has_header=False,
        new_columns=["anchor", "size", "covered", "sum", "mean0", "mean"],
        schema_overrides={"anchor": pl.Utf8, "covered": pl.Float64, "mean0": pl.Float64},
    )
    signal = df["mean0"].to_numpy()
    nonzero = (signal != 0) & (df["covered"].to_numpy() != 0)
    z = zscores(signal, nonzero)

    (
        pl.DataFrame({"anchor": df["anchor"], "zscore": z, "signal": np.where(nonzero, signal, 0.0)})
        .with_columns(rank=pl.col("zscore").rank("min", descending=True))
        .sort("anchor")
        .write_parquet(args.output)
    )


if __name__ == "__main__":
    main()
