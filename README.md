# ccre-pipeline-smk

Snakemake pipeline that calls ATAC-based cCREs.

![Pipeline DAG](docs/pipeline-dag.png)

## Run

From the repo root:

```sh
mamba env create -f environment.yaml
mamba activate ccre-pipeline
snakemake -n --configfile config-02_25m.yaml
snakemake --configfile config-02_25m.yaml --workflow-profile profile/slurm
```

`bigWigAverageOverBed` is read from `/zata/data/zlab/common/tools/ucsc.v385/`.

## Reference files

Not tracked in git. Place in `resources/`. MD5 checksums were verified against the sources on 2026-09-24.

`GRCh38-cCREs.bed` is the V7 `hg38-cCREs-Unfiltered.bed`, not `hg38-cCREs-Total.bed`.

| File | Source | MD5 |
|---|---|---|
| `GRCh38-cCREs.bed` | `/data/zusers/moorej3/moorej.ghpcc.project/ENCODE/Encyclopedia/V7/Registry/V7-hg38/hg38-cCREs-Unfiltered.bed` | `60c14d9f4903a7bd13b3db4f290ace48` |
| `GRCh38-Anchors.bed` | `/data/zusers/moorej3/moorej.ghpcc.project/ENCODE/Encyclopedia/V7/Registry/V7-hg38/hg38-Anchors.bed` | `9c5d57a3de8f09b27e9c00661801eddb` |
| `k100.Unique.Mappability.bed` | UCSC Umap, `https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k100.Unique.Mappability.bb`, converted with `bigBedToBed` | `eab6ef2b9961676b1f1cc0e84d88f5b3` |

```sh
D=/data/zusers/moorej3/moorej.ghpcc.project/ENCODE/Encyclopedia/V7/Registry/V7-hg38
cp $D/hg38-cCREs-Unfiltered.bed resources/GRCh38-cCREs.bed
cp $D/hg38-Anchors.bed resources/GRCh38-Anchors.bed
curl -O https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k100.Unique.Mappability.bb
bigBedToBed k100.Unique.Mappability.bb resources/k100.Unique.Mappability.bed
```

## Input data

Not copied. Read in place.

- CUSTOM ATAC narrowPeaks and bigWigs: `/zata/data/zlab/zusers/campbellm/t1d/`
- ENCODE data: `/data/projects/encode/`. The DNase rules use `/zata/data/zlab/projects/encode/data`, the same directory.

## Sample lists

`resources/CUSTOM_ATAC_List_{bw,narrowpeaks}_2026-08-25_0.2_25000000_orbonus_50per.txt` are copied from `/zata/zippy/campbellm/snakemake_atac_ccres/resources/` for `config-02_25m.yaml`. They are not tracked in git. They will be deprecated once test samples for pipeline development are chosen.
