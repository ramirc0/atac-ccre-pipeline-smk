#configfile: "config.yaml"

import os
import pandas as pd
from pathlib import Path

WORK_DIR = config["work_dir"]
RESULTS_DIR = config["ouput_dir"]   # keep if your config uses "ouput_dir"
LOG_DIR = os.path.join(WORK_DIR, config["log_dir"])
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

include: "workflow/ENCODE_make_filtered_list.smk"
include: "workflow/ENCODE_prepare_peaks.smk"
include: "workflow/CUSTOM_prepare_peaks.smk"
include: "workflow/merge_preped_peaks.smk"

include: "workflow/ATAC_cluster_rPeaks.smk"
include: "workflow/ATAC_filter_rPeaks_no_rdhs_mappable.smk"
include: "workflow/ATAC_make_summary_and_accession.smk"
include: "workflow/add_new_anchors.smk"

include: "workflow/CUSTOM_Zscore_across_bw.smk"
include: "workflow/ENCODE_Zscore_across_bw.smk"
include: "workflow/ENCODE_Zscore_across_DNase.smk"

include: "workflow/call_maxZ.smk"


include: "workflow/call_ATAC_cCREs.smk"


include: "workflow/call_maxZ_ATAC_extended_data.smk"
include: "workflow/plot_extended_data.smk"


# plots=[]
#if config['plot_top_and_bottom_n']: 
    

rule all:
    input:
        f"{RESULTS_DIR}/{PREFIX}_{GENOME}-Anchors-ATAC.bed",
        # intermediate_bw_files,
        # intermediate_DNase_files,
        f"{RESULTS_DIR}/{PREFIX}_{GENOME}-ATAC-maxZ.txt",
        f"{RESULTS_DIR}/{PREFIX}_{GENOME}-DNase-maxZ.txt",
        f"{RESULTS_DIR}/{PREFIX}_{GENOME}-cCREs.bed",
        #f"{RESULTS_DIR}/plots/plot_extended_data.done"
