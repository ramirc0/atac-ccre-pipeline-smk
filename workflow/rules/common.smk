# Config, sample lists, and shared path constants.

import os
import pandas as pd
import polars as pl
from pathlib import Path

WORK_DIR = config["work_dir"]
RESULTS_DIR = config["ouput_dir"]   # keep if your config uses "ouput_dir"
LOG_DIR = os.path.normpath(os.path.join(WORK_DIR, config["log_dir"]))
TMP_DIR = config["tmp_dir"]
TOOLKIT = config["TOOLKIT"]

Path(RESULTS_DIR).mkdir(parents=True, exist_ok=True)
Path(LOG_DIR).mkdir(parents=True, exist_ok=True)

GENOME = config["genome"]
FILTER_FRIP = config["filter_FRIP"]

ENCODE_INCLUDED = config.get("ENCODE_data_included", False)
CUSTOM_INCLUDED = config.get("CUSTOM_data_included", False)

if ENCODE_INCLUDED and CUSTOM_INCLUDED:
    PREFIX = "MERGED"
elif ENCODE_INCLUDED:
    PREFIX = "ENCODE"
elif CUSTOM_INCLUDED:
    PREFIX = "CUSTOM"
else:
    raise ValueError("At least one of ENCODE_data_included or CUSTOM_data_included must be True")

# mini script to add the generation of the zscore bws 
intermediate_bw_files = []

if bool(ENCODE_INCLUDED):
    encode_atac = pd.read_csv(
        config["encode_atac_list"],
        sep="\t",
        header=None,
        usecols=[0, 2],
        names=["experiment", "file"]
    )

    for _, row in encode_atac.iterrows():
        intermediate_bw_files.append(
            f"{RESULTS_DIR}/bw_zscoring/encode-{row['experiment']}_{row['file']}.txt"
        )

if bool(CUSTOM_INCLUDED):
    custom_bw = pd.read_csv(
        config["CUSTOM_data_bw_info"],
        sep="\t",
        header=None,
        usecols=[0]
    )

    for sample in custom_bw[0]:
        intermediate_bw_files.append(
            f"{RESULTS_DIR}/bw_zscoring/custom-{sample}.txt"
        )

# sample -> narrowPeak path, in list order.
CUSTOM_PEAKS = {}
if CUSTOM_INCLUDED:
    CUSTOM_PEAKS = dict(
        pl.read_csv(config["CUSTOM_data_macs_info"], separator="\t", infer_schema_length=0)
        .select("sample", "path")
        .iter_rows()
    )

intermediate_DNase_files = []

encode_dnase = pd.read_csv(
    config["encode_DNase_list"],
    sep="\t",
    header=None,
    usecols=[0, 1],
    names=["experiment", "file"]
)

for _, row in encode_dnase.iterrows():
    intermediate_DNase_files.append(
        f"{RESULTS_DIR}/bw_zscoring-DNAse/encode-{row['experiment']}_{row['file']}.txt"
    )

