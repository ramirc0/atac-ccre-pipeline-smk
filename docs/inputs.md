# Pipeline inputs

Keys are from `config/config.yaml.template`. The data flow diagram is `pipeline-dag.pdf`.

## Sample data

| Config key | File | Contents | Used by |
|---|---|---|---|
| `samples` | one TSV or a list, e.g. `config/samples-custom-2026-08-25.tsv`, `resources/encode-atac.2026-09-23.tsv` | One row per ATAC sample: `sample_id`, `narrowpeak`, `bigwig`, optional `biosample`. Fragment rows give `fragments` (and `fragment_format`) in place of narrowpeak and bigwig. | `prepare_peaks`, `zscore`, `fragment_call_peaks` |
| `dnase.list` | `resources/ENCODE-DNase-List.txt` | Tab-separated, no header: experiment, bigWig file accession, biosample. 1325 rows. Each bigWig is read from `{dnase.data_dir}/{experiment}/{file}.bigWig`. | `zscore` |

The custom and ENCODE sheets come from the preflight steps in the README.

## Reference files

| Config key | File | Role |
|---|---|---|
| `references.rdhs` | `resources/GRCh38-Anchors.bed` | Existing rDHS anchors. New ATAC rPeaks that overlap them are dropped. The rest are added to them. |
| `references.ccres` | `resources/GRCh38-cCREs.bed` | Existing Registry cCREs. Carried into the output with their accessions unchanged. |
| `references.blacklist` | `/data/zusers/ramirezc/static/ENCFF356LFX.bed` | ENCODE hg38 blacklist. Used by `filter_rpeaks`. |
| `references.chrom_sizes` | `/zata/data/zlab/common/genome/hg38.minimal.chrom.sizes` | Fragment rows only. Clips the MACS3 fold-enrichment bedGraph. |

Sources and checksums are in the README.

## Parameters

| Key | Effect |
|---|---|
| `run_id` | Output directory name under `outdir`, logs and benchmarks. |
| `genome` | Accession prefix (`EH38` for hg38). |
| `dnase` | `null` skips DNase z-scores and `DNase-maxZ.txt`. |
| `macs3` | `callpeak` options for fragment rows. |

## Outputs

All outputs are in `results/<run_id>/`:

- `Anchors-ATAC.bed`: rDHS anchors plus new ATAC anchors (`EH38A…`).
- `ATAC-maxZ.txt`, `DNase-maxZ.txt`: the maximum signal z-score of each anchor across all samples.
- `cCREs.bed`: existing Registry cCREs plus new `ATAC-cCRE` elements. A new ATAC anchor becomes a cCRE when its ATAC max-Z is above 1.64.
