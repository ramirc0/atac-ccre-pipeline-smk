# Config, sample sheet, and shared path constants.

from pathlib import Path

import polars as pl
from snakemake.exceptions import WorkflowError


configfile: "config/config.yaml"


RUN_ID = config.get("run_id", "default")
OUTDIR = f"{config['outdir']}/{RUN_ID}"
LOGDIR = f"logs/{RUN_ID}"
BENCHDIR = f"benchmarks/{RUN_ID}"
GENOME = config["genome"]
REFS = config["references"]

# The workspace path is relative to the rule file. Every rule file sits in rules/.
SOFTWARE_ENV = pixi(workspace="../envs", env=config["pixi_env"], locked=True)
QC_ENV = pixi(workspace="../envs", env="qc", locked=True)


# The pixi plugin installs the env inside each job. Parallel SLURM jobs then race
# on a fresh or relocked env, so install it once before any job starts.
onstart:
    shell("pixi install --manifest-path {workflow.basedir:q}/envs/pixi.toml --all --locked")


def _read_sheets(paths):
    """Concatenate one or more sample sheets.

    Parameters
    ----------
    paths : str or list of str
        Sample sheet TSVs.

    Returns
    -------
    polars.DataFrame
        Rows of every sheet, all cells as strings. A column missing from a
        sheet is null in its rows.

    Raises
    ------
    WorkflowError
        If a sheet does not exist.
    """
    paths = [paths] if isinstance(paths, str) else paths
    for p in paths:
        if not Path(p).exists():
            raise WorkflowError(
                f"Sample sheet not found: {p}\n"
                "Columns: sample_id, narrowpeak, bigwig, [biosample, fragments, "
                "fragment_format]; "
                "see config/samples.example.tsv."
            )
    return pl.concat(
        [pl.read_csv(p, separator="\t", infer_schema_length=0, comment_prefix="#") for p in paths],
        how="diagonal",
    )


_manifest = _read_sheets(config["samples"])
if _manifest.is_empty():
    raise WorkflowError(f"Sample sheets {config['samples']} have no rows.")
if _manifest["sample_id"].is_duplicated().any():
    raise WorkflowError(f"Duplicate sample_id in {config['samples']}.")


def _column(name):
    """Map each sample to its value in sheet column `name`.

    Parameters
    ----------
    name : str
        Sample sheet column.

    Returns
    -------
    dict of str to str or None
        Value per `sample_id`, with empty cells as None. Empty when the sheet
        has no such column.
    """
    if name not in _manifest.columns:
        return {}
    return {
        row["sample_id"]: (row[name] or None)
        for row in _manifest.iter_rows(named=True)
    }


SAMPLES = _manifest["sample_id"].to_list()
FRAGMENTS_OF = {s: f for s, f in _column("fragments").items() if f}
FRAGMENT_FORMAT_OF = _column("fragment_format")

# Fragment rows get their narrowPeak and bigWig from rules/fragments.smk.
NARROWPEAK_OF = {
    **_column("narrowpeak"),
    **{s: f"{OUTDIR}/fragments/{s}/{s}_peaks.narrowPeak" for s in FRAGMENTS_OF},
}
BIGWIG_OF = {
    **_column("bigwig"),
    **{s: f"{OUTDIR}/fragments/{s}/{s}_FE.bw" for s in FRAGMENTS_OF},
}
_incomplete = [s for s in SAMPLES if not (NARROWPEAK_OF.get(s) and BIGWIG_OF.get(s))]
if _incomplete:
    raise WorkflowError(f"Samples need fragments, or narrowpeak and bigwig: {_incomplete}")

_biosample = _column("biosample")
BIOSAMPLE_OF = {s: _biosample.get(s) or s for s in SAMPLES}

# DNase experiment -> bigWig; empty when `dnase` is null.
DNASE_BIGWIG_OF = {}
if config.get("dnase"):
    DNASE_BIGWIG_OF = {
        experiment: f"{config['dnase']['data_dir']}/{experiment}/{file}.bigWig"
        for experiment, file in pl.read_csv(
            config["dnase"]["list"], separator="\t", has_header=False, infer_schema_length=0
        )
        .select(pl.nth(0, 1))
        .iter_rows()
    }

# assay -> id -> bigWig, for the z-score rule.
BIGWIGS = {"ATAC": BIGWIG_OF, "DNase": DNASE_BIGWIG_OF}


def macs3_flags(sample):
    """Return the macs3 callpeak flag tokens of a fragment sample.

    Parameters
    ----------
    sample : str
        Sample ID with a `fragments` column.

    Returns
    -------
    list of str
        Flags from `config["macs3"]` and the sample's fragment format.
        Interpolate them with `:q`.

    Raises
    ------
    WorkflowError
        If the fragment format is not FRAG or BEDPE.
    """
    m = config["macs3"]
    fmt = (FRAGMENT_FORMAT_OF.get(sample) or "FRAG").upper()
    if fmt not in ("FRAG", "BEDPE"):
        raise WorkflowError(f"{sample}: fragment_format must be FRAG or BEDPE, not {fmt}")
    flags = ["-f", fmt, "-g", str(m["genome_size"]), "-q", str(m["qvalue"])]
    if fmt == "FRAG" and m["max_count"] != "all":
        flags += ["--max-count", str(m["max_count"])]
    return flags + list(m["extra"])


wildcard_constraints:
    assay="ATAC|DNase",
    id=r"[A-Za-z0-9][A-Za-z0-9_.+-]*",
