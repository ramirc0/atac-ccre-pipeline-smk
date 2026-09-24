"""Standardize narrowPeak files into one unsorted peak BED.

Keeps chr1-22, X, Y. Peaks narrower than 150 bp or wider than 350 bp are
re-centered on the summit (+/-75 or +/-175). Appends exp, target, biosample and
a per-sample unique id `<sample>-<row>`, row numbered before filtering.
"""

import argparse

import polars as pl

SCHEMA = {
    "chrom": pl.Utf8,
    "start": pl.Int64,
    "end": pl.Int64,
    "name": pl.Utf8,
    "score": pl.Float64,
    "strand": pl.Utf8,
    "signalValue": pl.Float64,
    "pValue": pl.Float64,
    "qValue": pl.Float64,
    "peak": pl.Int64,
}


def build_parser():
    """Return the argument parser for prepare_peaks.py."""
    p = argparse.ArgumentParser(description="Standardize narrowPeak files.")
    p.add_argument("--peaks", nargs="+", required=True, help="narrowPeak files.")
    p.add_argument("--samples", nargs="+", required=True, help="Sample id per file.")
    p.add_argument("--biosamples", nargs="+", help="Biosample per file (default: sample id).")
    p.add_argument("--target", default="ATAC", help="Target column value.")
    p.add_argument("-o", "--output", required=True, help="Peak BED to write.")
    return p


def standardize(path, sample, biosample, target):
    """Read one narrowPeak file and return its standardized peaks."""
    width = pl.col("end") - pl.col("start")
    summit = pl.col("start") + pl.col("peak")
    return (
        pl.read_csv(
            path,
            separator="\t",
            has_header=False,
            new_columns=list(SCHEMA),
            schema_overrides=SCHEMA,
        )
        .with_row_index("row_nr", offset=1)
        .filter(pl.col("chrom").str.contains(r"^chr([1-9]|1[0-9]|2[0-2]|X|Y)$"))
        .with_columns(
            pl.when(width < 150)
            .then(summit - 75)
            .when(width > 350)
            .then(summit - 175)
            .otherwise(pl.col("start"))
            .alias("start"),
            pl.when(width < 150)
            .then(summit + 75)
            .when(width > 350)
            .then(summit + 175)
            .otherwise(pl.col("end"))
            .alias("end"),
            pl.lit(sample).alias("exp"),
            pl.lit(target).alias("target"),
            pl.lit(biosample).alias("biosample"),
            (pl.lit(f"{sample}-") + pl.col("row_nr").cast(pl.Utf8)).alias("unique_id"),
        )
        .drop("row_nr")
        .filter(pl.col("start") < pl.col("end"))
    )


def main(argv=None):
    """Standardize every sample's peaks and write them in input order."""
    args = build_parser().parse_args(argv)
    biosamples = args.biosamples or args.samples
    if not len(args.peaks) == len(args.samples) == len(biosamples):
        raise SystemExit("--peaks, --samples and --biosamples must have the same length.")

    with open(args.output, "w") as out:
        for path, sample, biosample in zip(args.peaks, args.samples, biosamples):
            print(f"{sample}: {path}")
            standardize(path, sample, biosample, args.target).write_csv(
                out, separator="\t", include_header=False
            )


if __name__ == "__main__":
    main()
