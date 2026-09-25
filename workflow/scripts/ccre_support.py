"""Support table for the new ATAC anchors.

One row per new anchor (EH38A... for hg38): coordinates, whether it became an
ATAC cCRE, ATAC and DNase max z-scores, the number of ATAC samples with z-score
>= threshold, the number of experiments whose summits support its rPeak, the
sample whose peak won the rPeak, and the midpoint distance to the nearest
registry anchor on the left and right.
"""

import argparse

import numpy as np
import polars as pl
from accession_regions import GENOME_IDS, MODE_IDS
from call_ccres import read_bed

BED = ["chrom", "start", "end", "anchor"]


def build_parser():
    """Return the argument parser for ccre_support.py."""
    p = argparse.ArgumentParser(description="Support table for the new ATAC anchors.")
    p.add_argument("--anchors", required=True, help="Anchors BED: chrom, start, end, anchor.")
    p.add_argument("--summary", required=True, help="ATAC-Summary.txt from accession_rpeaks.")
    p.add_argument("--ccres", required=True, help="cCREs.bed from call_ccres.")
    p.add_argument("--atac-max-z", required=True, help="ATAC anchor/max_zscore TSV.")
    p.add_argument("--dnase-max-z", nargs="*", default=[], help="DNase anchor/max_zscore TSV.")
    p.add_argument("--zscores", nargs="+", required=True, help="ATAC zscore.py parquet files.")
    p.add_argument("--genome", required=True, choices=GENOME_IDS, help="Genome build.")
    p.add_argument("--threshold", type=float, default=1.64, help="Per-sample z-score cutoff.")
    p.add_argument("-o", "--output", required=True, help="Support TSV to write.")
    return p


def read_max_z(path, name):
    """Read an anchor/max_zscore TSV, renaming the score to `name`."""
    return pl.read_csv(path, separator="\t", schema_overrides={"anchor": pl.Utf8}).rename(
        {"max_zscore": name}
    )


def count_passing(paths, prefix, threshold):
    """Per new anchor, the number of files whose z-score is >= threshold."""
    anchor = pl.read_parquet(paths[0], columns=["anchor"])["anchor"]
    is_new = anchor.str.starts_with(prefix).to_numpy()
    counts = np.zeros(is_new.sum(), dtype=np.int32)
    for path in paths:
        df = pl.read_parquet(path, columns=["anchor", "zscore"])
        if not df["anchor"].equals(anchor):
            raise ValueError(f"Anchor order mismatch in {path}")
        counts += df["zscore"].to_numpy()[is_new] >= threshold
    return pl.DataFrame({"anchor": anchor.filter(is_new), "samples_passing": counts})


def nearest_registry(anchors, is_new):
    """Midpoint distance from each new anchor to the nearest registry anchor per side."""
    mid = ((pl.col("start") + pl.col("end")) // 2).alias("mid")
    registry = anchors.filter(~is_new).select("chrom", mid.alias("registry_mid")).sort("registry_mid")
    new = anchors.filter(is_new).select("anchor", "chrom", mid).sort("mid")

    def side(strategy):
        return new.join_asof(
            registry,
            left_on="mid",
            right_on="registry_mid",
            by="chrom",
            strategy=strategy,
            check_sortedness=False,
        )

    left = side("backward").select("anchor", left_registry_bp=pl.col("mid") - pl.col("registry_mid"))
    right = side("forward").select("anchor", right_registry_bp=pl.col("registry_mid") - pl.col("mid"))
    return left.join(right, on="anchor")


def main(argv=None):
    """Join every per-anchor measure onto the new anchors and write the table."""
    args = build_parser().parse_args(argv)
    prefix = GENOME_IDS[args.genome] + MODE_IDS["ATAC"]

    anchors = read_bed(args.anchors, BED)
    is_new = pl.col("anchor").str.starts_with(prefix)
    ccres = pl.read_csv(args.ccres, separator="\t", columns=["anchor", "type"]).filter(
        pl.col("type") == "ATAC-cCRE"
    )
    summary = pl.read_csv(
        args.summary,
        separator="\t",
        has_header=False,
        columns=[3, 7, 8],
        new_columns=["rpeak", "peak_experiments", "anchor"],
        quote_char=None,
    ).select(
        "anchor",
        "peak_experiments",
        rpeak_sample=pl.col("rpeak").str.replace(r"-\d+$", ""),
    )

    table = (
        anchors.filter(is_new)
        .with_columns(is_ccre=pl.col("anchor").is_in(ccres["anchor"].implode()))
        .join(read_max_z(args.atac_max_z, "atac_max_z"), on="anchor", how="left")
    )
    for path in args.dnase_max_z:
        table = table.join(read_max_z(path, "dnase_max_z"), on="anchor", how="left")
    table = (
        table.join(count_passing(args.zscores, prefix, args.threshold), on="anchor", how="left")
        .join(summary, on="anchor", how="left")
        .join(nearest_registry(anchors, is_new), on="anchor", how="left")
    )
    table.write_csv(args.output, separator="\t")


if __name__ == "__main__":
    main()
