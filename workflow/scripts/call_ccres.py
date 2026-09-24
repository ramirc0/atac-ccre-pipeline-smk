"""Add ATAC cCREs to the registry cCREs and assign accessions.

An ATAC anchor (EH38A...) becomes an ATAC cCRE when its max z-score exceeds the
threshold. Regions whose coordinates match a registry cCRE keep that accession;
the rest get new accessions numbered after the registry's highest, in order.
"""

import argparse

import polars as pl

GENOME_IDS = {"hg19": "EH37", "hg38": "EH38", "mm10": "EM10"}
MODE_ID = "E"
COORDS = ["chrom", "start", "end"]


def build_parser():
    """Return the argument parser for call_ccres.py."""
    p = argparse.ArgumentParser(description="Call ATAC cCREs and assign accessions.")
    p.add_argument("--max-z", required=True, help="ATAC anchor/max_zscore TSV.")
    p.add_argument("--anchors", required=True, help="Anchors BED: chrom, start, end, anchor.")
    p.add_argument("--ccres", required=True, help="Registry cCRE BED (6 columns).")
    p.add_argument("--genome", required=True, choices=GENOME_IDS, help="Genome build.")
    p.add_argument("--threshold", type=float, default=1.64, help="Max z-score cutoff.")
    p.add_argument("-o", "--output", required=True, help="cCRE BED with header to write.")
    return p


def read_bed(path, names):
    """Read the first len(names) columns of a headerless BED."""
    return pl.read_csv(
        path,
        separator="\t",
        has_header=False,
        columns=list(range(len(names))),
        new_columns=names,
        schema_overrides={"chrom": pl.Utf8, "start": pl.Int64, "end": pl.Int64},
    )


def main(argv=None):
    """Select ATAC cCREs, merge with the registry, assign accessions, write sorted."""
    args = build_parser().parse_args(argv)
    new_prefix = GENOME_IDS[args.genome] + MODE_ID

    registry = read_bed(args.ccres, ["chrom", "start", "end", "anchor", "ccre", "type"])
    max_accession = registry["ccre"].str.split(MODE_ID).list.last().cast(pl.Int64).max()

    atac = (
        pl.read_csv(args.max_z, separator="\t", schema_overrides={"anchor": pl.Utf8})
        .filter(
            (pl.col("max_zscore") > args.threshold)
            & pl.col("anchor").str.contains("EH38A")
        )
        .join(
            read_bed(args.anchors, ["chrom", "start", "end", "anchor"]),
            on="anchor",
            maintain_order="left",
        )
        .select(*COORDS, "anchor", type=pl.lit("ATAC-cCRE"))
    )

    # Last registry entry wins on duplicate coordinates, as a dict would.
    previous = registry.unique(COORDS, keep="last", maintain_order=True).select(
        *COORDS, previous=pl.col("ccre").str.replace("EH37", "EH38")
    )
    regions = pl.concat([atac, registry.select(*COORDS, "anchor", "type")]).join(
        previous, on=COORDS, how="left", maintain_order="left"
    )

    is_new = pl.col("previous").is_null()
    number = max_accession + is_new.cum_sum()
    out = regions.select(
        *COORDS,
        "anchor",
        ccre=pl.when(is_new)
        .then(pl.lit(new_prefix) + number.cast(pl.Utf8).str.zfill(7))
        .otherwise(pl.col("previous")),
        type="type",
    ).sort(COORDS, maintain_order=True)

    out.write_csv(args.output, separator="\t")


if __name__ == "__main__":
    main()
