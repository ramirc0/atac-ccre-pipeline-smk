# Pipeline inputs

Applies to `Snakefile` with `config-02_25m.yaml`. The data flow diagram is `pipeline-dag.pdf`.

## Sample data

| Config key | File | Contents | Used by |
|---|---|---|---|
| `CUSTOM_data_macs_info` | `resources/CUSTOM_ATAC_List_narrowpeaks_*.txt` | Tab-separated, with a header: `sample`, `path` to a MACS3 narrowPeak. 122 samples. | `CUSTOM_prepare_peaks` |
| `CUSTOM_data_bw_info` | `resources/CUSTOM_ATAC_List_bw_*.txt` | Tab-separated, no header: sample, path to a fold-enrichment bigWig. Same 122 samples. | `CUSTOM_Zscore_across_bw` |
| `encode_DNase_list` | `resources/ENCODE-DNase-List.txt` | Tab-separated, no header: experiment, bigWig file accession, biosample. 1325 rows. Each bigWig is read from `/zata/data/zlab/projects/encode/data/{experiment}/{file}.bigWig`, or downloaded if it is missing. | `ENCODE_Zscore_across_DNase` |
| `encode_atac_list` | Only used when `ENCODE_data_included: True`. | ENCODE ATAC experiments, peaks and bigWigs. | `ENCODE_prepare_ATAC_peaks`, `ENCODE_Zscore_across_bw` |

The narrowPeak list has a header and the bigWig list does not.

## Reference files

| Config key | File | Role |
|---|---|---|
| `rdhs_path` | `resources/GRCh38-Anchors.bed` | Existing rDHS anchors. New ATAC peaks that overlap them are dropped. The rest are added to them. |
| `ccre_path` | `resources/GRCh38-cCREs.bed` | Existing Registry cCREs. Carried into the output with their accessions unchanged. |
| `mappable` | `resources/k100.Unique.Mappability.bed` | Declared as an input of `ATAC_filter_rPeaks` but not used by its shell command. |
| hard-coded | `/data/zusers/ramirezc/static/ENCFF356LFX.bed` | ENCODE hg38 blacklist. Used by `ATAC_filter_rPeaks`. |
| hard-coded | `/data/projects/encode/Registry/V4/GRCh38/GRCh38-MultiMap-cCREs.bed` | Declared as an input of `ATAC_filter_rPeaks` but not used by its shell command. |

Sources and checksums are in the README.

## Parameters

| Key | Effect |
|---|---|
| `ENCODE_data_included`, `CUSTOM_data_included` | Choose the ATAC peak sources. The output prefix is `ENCODE`, `CUSTOM` or `MERGED`. |
| `filter_FRIP`, `filter_min_reads` | Filter the ENCODE ATAC list. They don't affect CUSTOM samples, which are filtered before the list is made. |
| `ouput_dir`, `log_dir`, `tmp_dir` | Output, log and scratch locations. `ouput_dir` is misspelled in the code too. |

## Outputs

All outputs are in `ouput_dir`, prefixed `{PREFIX}_hg38-`:

- `Anchors-ATAC.bed`: rDHS anchors plus new ATAC anchors (`EH38A…`).
- `ATAC-maxZ.txt`, `DNase-maxZ.txt`: the maximum signal z-score of each anchor across all samples.
- `cCREs.bed`: existing Registry cCREs plus new `ATAC-cCRE` elements. A new ATAC anchor becomes a cCRE when its ATAC max-Z is above 1.64.
