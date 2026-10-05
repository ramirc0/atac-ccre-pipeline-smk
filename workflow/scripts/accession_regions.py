"""Append a new accession to every region, numbered in input order.

Accessions are <genome id><mode id><7-digit number>, starting at 1: for hg38
ATAC regions, EH38A0000001, EH38A0000002, ...
"""

import argparse
import sys

import polars as pl

GENOME_IDS = {"hg19": "EH37", "hg38": "EH38", "mm10": "EM10"}
MODE_IDS = {
    "rDHS": "D",
    "cCRE": "E",
    "TF": "F",
    "MultiMap": "M",
    "ATAC": "A",
    "BETA": "B",
}


def build_parser():
    """Build the accession_regions.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the accession_regions.py flags.
    """
    p = argparse.ArgumentParser(description="Append accessions to regions.")
    p.add_argument("regions", help="Headerless TSV of regions.")
    p.add_argument("--genome", required=True, choices=GENOME_IDS, help="Genome build.")
    p.add_argument("--mode", required=True, choices=MODE_IDS, help="Region type.")
    return p


def main(argv=None):
    """Print each region with its accession as a last column.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    prefix = GENOME_IDS[args.genome] + MODE_IDS[args.mode]
    regions = pl.read_csv(
        args.regions,
        separator="\t",
        has_header=False,
        infer_schema_length=0,
        quote_char=None,
    )
    regions.with_columns(
        accession=pl.lit(prefix)
        + pl.int_range(1, regions.height + 1).cast(pl.Utf8).str.zfill(7)
    ).write_csv(sys.stdout, separator="\t", include_header=False, quote_style="never")


if __name__ == "__main__":
    main()
