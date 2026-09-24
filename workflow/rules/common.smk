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

# Absolute so conda: resolves the same from any rule file.
CONDA_ENV = str((Path(workflow.basedir).parent / config["conda_env"]).resolve())


def _read_sheets(paths):
    """Concatenate one or more sample sheets; all cells read as strings."""
    paths = [paths] if isinstance(paths, str) else paths
    for p in paths:
        if not Path(p).exists():
            raise WorkflowError(
                f"Sample sheet not found: {p}\n"
                "Columns: sample_id, narrowpeak, bigwig, [biosample]; "
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
    """sample_id -> value for `name`, empty string treated as None."""
    if name not in _manifest.columns:
        return {}
    return {
        row["sample_id"]: (row[name] or None)
        for row in _manifest.iter_rows(named=True)
    }


SAMPLES = _manifest["sample_id"].to_list()
NARROWPEAK_OF = _column("narrowpeak")
BIGWIG_OF = _column("bigwig")
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


wildcard_constraints:
    assay="ATAC|DNase",
    id=r"[A-Za-z0-9][A-Za-z0-9_.+-]*",
