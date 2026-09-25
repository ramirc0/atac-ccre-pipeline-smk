"""Select ENCODE ATAC-seq experiments from an encode-metadata snapshot.

Experiments: released, unperturbed (unless --include-perturbed) Homo sapiens
ATAC-seq. Per experiment:
  peak    the released GRCh38 preferred-default narrowPeak
  bigWig  the released GRCh38 fold-change bigWig from the experiment's default
          analysis that covers the most biological replicates (the pooled track)
Filters on the peak file's FRiP and, optionally, read depth. Fails if a chosen
file is missing from the local ENCODE mirror.

Writes a sample sheet (header; sample_id, narrowpeak, bigwig, biosample) and,
optionally, the headerless list the original pipeline reads (experiment, peak,
bigWig, biosample).
"""

import argparse
import sys
from pathlib import Path

import polars as pl

PEAK_DIR = "/data/projects/encode/data"
BIGWIG_DIR = "/zata/data/zlab/projects/encode/data"


def build_parser():
    """Return the argument parser for encode_atac_samples.py."""
    p = argparse.ArgumentParser(description="ENCODE ATAC samples from a metadata snapshot.")
    p.add_argument("snapshot", type=Path, help="Snapshot parquets dir (results/<run_id>/parquets).")
    p.add_argument("-o", "--output", required=True, help="Sample sheet TSV to write.")
    p.add_argument("--legacy-output", help="Headerless list for the original pipeline.")
    p.add_argument("--min-frip", type=float, default=0.2, help="Minimum FRiP of the peak file.")
    p.add_argument("--min-reads", type=float, default=None, help="Minimum usable fragments.")
    p.add_argument("--include-perturbed", action="store_true", help="Keep perturbed experiments.")
    p.add_argument("--peak-dir", default=PEAK_DIR, help="Local ENCODE peak mirror.")
    p.add_argument("--bigwig-dir", default=BIGWIG_DIR, help="Local ENCODE bigWig mirror.")
    return p


def select(snapshot, min_frip, min_reads, include_perturbed=False):
    """One row per experiment: experiment, peak, bigwig, biosample, frip, reads."""
    experiments = pl.read_parquet(
        snapshot / "experiment.parquet",
        columns=["accession", "assay_title", "status", "organism", "perturbed",
                 "biosample_summary", "default_analysis"],
    ).filter(
        (pl.col("assay_title") == "ATAC-seq")
        & (pl.col("status") == "released")
        & (pl.col("organism") == "Homo sapiens")
        & (pl.lit(include_perturbed) | ~pl.col("perturbed"))
    )
    files = (
        pl.read_parquet(
            snapshot / "file.parquet",
            columns=["id", "accession", "dataset", "file_type", "output_type", "status",
                     "genome_assembly", "preferred_default", "biological_replicates", "analyses"],
        )
        .filter((pl.col("status") == "released") & (pl.col("genome_assembly") == "GRCh38"))
        .with_columns(experiment=pl.col("dataset").str.extract(r"^/experiments/(\w+)/$"))
        .join(
            experiments.select(experiment="accession", default_analysis="default_analysis"),
            on="experiment",
        )
    )

    peaks = files.filter(
        (pl.col("file_type") == "bed narrowPeak") & pl.col("preferred_default").fill_null(False)
    ).select("experiment", peak="accession", peak_id="id")

    bigwigs = (
        files.filter(
            (pl.col("output_type") == "fold change over control")
            & pl.col("analyses").str.contains(pl.col("default_analysis"), literal=True)
        )
        .with_columns(n_reps=pl.col("biological_replicates").str.count_matches(",") + 1)
        .sort("experiment", "n_reps", "accession", descending=[False, True, False])
        .group_by("experiment", maintain_order=True)
        .first()
        .select("experiment", bigwig="accession")
    )

    frip = pl.read_parquet(
        snapshot / "atac_peak_enrichment_quality_metric.parquet",
        columns=["files", "frip_score_of_peaks", "status"],
    ).filter(pl.col("status") == "released").select(
        peak_id="files", frip="frip_score_of_peaks"
    )

    reads = (
        pl.read_parquet(
            snapshot / "atac_alignment_quality_metric.parquet",
            columns=["files", "usable_fragments", "status"],
        )
        .filter(pl.col("status") == "released")
        .join(files.select(files="id", experiment="experiment"), on="files")
        .group_by("experiment")
        .agg(reads=pl.col("usable_fragments").max())
    )

    duplicated = peaks.filter(pl.col("experiment").is_duplicated())
    if duplicated.height:
        raise SystemExit(f"Experiments with several default peaks:\n{duplicated}")

    out = (
        experiments.select(experiment="accession", biosample="biosample_summary")
        .join(peaks, on="experiment", how="left")
        .join(bigwigs, on="experiment", how="left")
        .join(frip, on="peak_id", how="left")
        .join(reads, on="experiment", how="left")
        .drop("peak_id")
        .sort("experiment")
    )
    incomplete = out.filter(pl.any_horizontal(pl.col("peak", "bigwig", "frip").is_null()))
    if incomplete.height:
        raise SystemExit(f"Experiments missing a peak, bigWig or FRiP:\n{incomplete}")
    keep = pl.col("frip") >= min_frip
    if min_reads is not None:
        keep &= pl.col("reads") >= min_reads
    print(f"{experiments.height} experiments; {out.filter(keep).height} pass filters", file=sys.stderr)
    return out.filter(keep)


def main(argv=None):
    """Select experiments, check local files, and write the sheet(s)."""
    args = build_parser().parse_args(argv)
    rows = select(args.snapshot, args.min_frip, args.min_reads, args.include_perturbed).with_columns(
        narrowpeak=pl.format("{}/{}/{}.bed.gz", pl.lit(args.peak_dir), "experiment", "peak"),
        bigwig_path=pl.format("{}/{}/{}.bigWig", pl.lit(args.bigwig_dir), "experiment", "bigwig"),
    )

    missing = [p for p in [*rows["narrowpeak"], *rows["bigwig_path"]] if not Path(p).exists()]
    if missing:
        raise SystemExit(f"{len(missing)} files missing locally:\n" + "\n".join(missing))

    rows.select(
        sample_id="experiment", narrowpeak="narrowpeak", bigwig="bigwig_path", biosample="biosample"
    ).write_csv(args.output, separator="\t")
    if args.legacy_output:
        rows.select("experiment", "peak", "bigwig", "biosample").write_csv(
            args.legacy_output, separator="\t", include_header=False
        )


if __name__ == "__main__":
    main()
