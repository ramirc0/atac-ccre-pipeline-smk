# Glossary: manuscript ↔ this pipeline

Paper: Moore, J.E. et al. (2026). "An expanded registry of candidate *cis*-regulatory elements." *Nature*. https://doi.org/10.1038/s41586-025-09909-9
Paper code: https://github.com/weng-lab/ENCODE-cCREs, `Version-4/cCRE-Pipeline/` (below `P/`). Local clone: `/data/zusers/moorej3/GitHub/ENCODE-cCREs`.

## Versions

| Term | Meaning |
|---|---|
| `V7` in lab paths (`Encyclopedia/V7/Registry/V7-hg38`) | Registry **V4**, the version in the paper. `V5` = Registry V2, `V6` = Registry V3. |
| `hg38-cCREs-Unfiltered.bed` (our `GRCh38-cCREs.bed`) | The final cell type-agnostic V4 Registry: 2,348,854 cCREs. "Unfiltered" is a historical name. |
| `hg38-Anchors.bed` (our `GRCh38-Anchors.bed`) | All V4 anchors: rDHSs + TF rClusters, 2,950,228. |

## Regions

| Paper term | This pipeline | Notes |
|---|---|---|
| DHS (DNase hypersensitive site) | none | DNase peak. ATAC counterpart: OCR. |
| rDHS (representative DHS) | none | Anchor with an `EH38D` accession. ATAC counterpart: rOCR. |
| **OCR** (open chromatin region) | none | Not a term in the paper. ATAC counterpart of a DHS. Not used for the input narrowPeaks. |
| **rOCR** (representative OCR) | `cluster_rpeaks` output (`rpeaks/rpeaks.bed`) | Not a term in the paper. Built with the same iterative merge as rDHSs and TF rClusters (`pick_best_peak.py`, a port of `P/Toolkit/pick-best-peak.py`). |
| rPeaks | `rpeaks/rpeaks.bed` | Generic "representative peaks" from the iterative merge. Used for both rDHSs and TF rClusters in `P/`. |
| TF rClusters | none | TF ChIP-seq rPeaks supported by >= 5 experiments (`filter-tf-rpeaks.py`; here `filter_rpeaks.py`). Anchors with an `EH38F` accession. Our ATAC path copies this method. |
| rAPeaks | `filter_rpeaks` output (`rpeaks/no-rdhs-blacklist.tsv`) | Term from the lab slides: rOCRs with >= 5 experiments that don't overlap an anchor. |
| anchors | `references.rdhs`, `anchors` output (`Anchors-ATAC.bed`) | Paper: rDHSs + TF rClusters. The pipeline adds new ATAC anchors (`EH38A`). |
| multi-mapping cCREs | `GRCh38-MultiMap-cCREs.bed` (declared but unused in the original `ATAC_filter_rPeaks`; removed) | `EH38M` accessions. See "Located files". |

## Accessions (`make-region-accession.py`; here `accession_regions.py`)

`EH38` + type + 7 digits. `D` rDHS, `F` TF rCluster, `E` cCRE, `M` multi-map (all in the paper's pipeline). This repo adds `A` (ATAC anchor) and `B` (BETA). Accessions are reused only when coordinates are identical.

## Signal

| Paper term | This pipeline | Notes |
|---|---|---|
| Z-score of the log-transformed signal | `zscore.py` (was `log-zscore-AF-MC.py`) | Same logic as `P/Toolkit/log-zscore-normalization.py`, vectorized; values are bit-identical. log10 of `bigWigAverageOverBed` mean. Zero-signal anchors are excluded from the mean and SD and get z = -10. |
| max-Z | `max_z` → `ATAC-maxZ.txt`, `DNase-maxZ.txt` | Max z-score per anchor across all experiments. Keyed on the anchor ID, not the cCRE ID. |
| high signal | `max_zscore > 1.64` | The only threshold used in both codebases. |
| ">= 5 experiments" | `filter_rpeaks.py --min-experiments 5` | Counts distinct experiments. For CUSTOM data, experiment = sample. |

## Classes

| Paper class | Code label | This pipeline |
|---|---|---|
| none | `ATAC-cCRE` | **ATAC cCREs**: new ATAC anchors with ATAC max-Z > 1.64, from `call_ATAC_cCREs`. Preferred term. Not the paper's CA class, which also requires low H3K4me3, H3K27ac and CTCF. |
| PLS, pELS, dELS, CA-H3K4me3, CA-CTCF, CA-TF, TF | same strings | Copied unchanged from the Registry. |

## Located files

The paper-code glossary listed these as missing from the public repo. All found under `/data/zusers/moorej3` (below `M/`) on 2026-09-24. `V/` = `M/moorej.ghpcc.project/ENCODE/Encyclopedia/V7/Registry/V7-hg38`.

### Multi-map cCREs

| File | Location | Notes |
|---|---|---|
| Generating script | `M/Projects/ENCODE/Encyclopedia/Version7/cCRE-Pipeline/0.5_Call-MultiMap-DHSs.sh` | Step 0.5 in the lab's working pipeline, not in the public repo. Runs the DHS iterative-threshold calling on multimapper hotspot calls. Hard-coded to mm10. |
| Input list | `V/hg38-rDHS/hg38-Multimapper-List.txt` | Per-experiment output in `V/hg38-rDHS/MultiMappers/*.MM.bed`. |
| `hg38-MultiMap-cCREs.bed` (raw) | `V/hg38-MultiMap-cCREs.RAW.bed` | 24,160 regions. Same MD5 as `/data/projects/encode/Registry/V4/GRCh38/GRCh38-MultiMap-cCREs.bed`, the file our rule declares. |
| `hg38-MultiMap-cCREs.bed` (accessioned) | `V/hg38-cCREs-MultiMap.bed` | Same 24,160 regions with `EH38M` and `EH38E` accessions. |
| `hg38-Multimappers.bed` | `V/hg38-rDHS/hg38-Multimappers.bigBed`, `V/hg38-rDHS/hg38-All-MultiMap-rDHS.bed` | Only these forms exist. |

### Analysis data

| File | Location |
|---|---|
| `REST-Enhancers.bed`, `REST-Silencers.bed` | `V/Manuscript-Analysis/3_Silencers/`. Slightly different copies in `V/Manuscript-Analysis/REVISIONS/REST/Aggregation/`. |
| `STARR-Silencers.Stringent.bed`, `STARR-Silencers.Robust.bed` | `V/Manuscript-Analysis/3_Silencers/` |
| `Silencer-cCREs.bed` | `V/Manuscript-Analysis/3_Silencers/` |
| `biosample-classification.txt` | `V/Manuscript-Analysis/1_Updated-Registry/` (same size in `3_Silencers/`, `Classification-Expansion/`) |
| `K562-cCREs.bed` | Two different files: 18 MB in `V/Manuscript-Analysis/2_Functional-Characterization/cCRE-Enrichment/`, 258 MB in `V/Manuscript-Analysis/3_Silencers/Published-Silencers/Histone-Overlap/`. Check which one a script expects. |

### Helper scripts

| File | Location |
|---|---|
| `count.py`, `count-with-sum.py`, `plot-zoonomia-triangle.R`, `fisher-test-repeats-multimap.py`, `calculate-tf-motif-gc.py`, `gc-match-tf.py`, `starr-pca.py`, `assign-organ-tissue.py` | `M/GitHub/ENCODE-cCREs/Version-4/cCRE-Analysis/Toolkit/` (untracked files in the local clone, never committed) |
| `Aggregate-Signal-Batch-1.sh` | `M/Projects/ENCODE/Encyclopedia/Version7/cCRE-Analysis/` |
| `pull-HiC.py` | `V/Manuscript-Analysis/Encyclopedia-Overlap/Hi-C/` (latest of 3 copies) |
| `fisher.test.py`, `padjust.R` | `M/bin/`. Many older `fisher.test.py` copies exist. |
