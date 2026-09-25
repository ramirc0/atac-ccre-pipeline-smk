# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Snakemake pipeline that calls ATAC cCREs on top of the ENCODE Registry V7 (hg38). Read `README.md`
for setup, sample sheets, reference files and preflight scripts, and `docs/inputs.md` for config
keys and output formats.

## Commands

The driver env is `/zata/zippy/ramirezc/Miniforge3/envs/snakemake` (Snakemake 9, polars). Rule tools
come from `workflow/envs/env.yaml`, which the profiles build into `.conda/`.

```sh
snakemake -n -p --profile profiles/local        # dry run
snakemake --profile profiles/slurm              # full run on SLURM
snakemake --profile profiles/local results/<run_id>/ATAC.bed   # one target
snakemake --lint                                # style check; passes today
```

`config/config.yaml` and `config/*.tsv` are gitignored. There are no unit tests. Changes are
verified by running the pipeline and diffing outputs against a previous run.

## Architecture

Data flow (see `docs/pipeline-dag.png`):

1. `fragments.smk` (fragment rows only): MACS3 turns a fragments file into the narrowPeak and
   fold-enrichment bigWig that bulk rows provide directly.
2. `peaks.smk`: all narrowPeaks are pooled (`prepare_peaks`) and clustered into rPeaks by an
   iterative bedtools merge loop that keeps the strongest peak per cluster (`cluster_rpeaks`).
   rPeaks supported by >= 5 experiments and off the registry rDHSs and blacklist get new `EH38A...`
   accessions and are added to the rDHS anchors (`Anchors-ATAC.bed`).
3. `zscore.smk`: one `zscore` job per bigWig, ATAC samples and ENCODE DNase experiments alike
   (wildcards `assay`, `id`). This is the wide fan-out (~1700 jobs for ENCODE + custom).
4. `ccres.smk`: per-assay max-Z; an ATAC anchor with ATAC max-Z > 1.64 becomes an `ATAC-cCRE`.
   DNase max-Z is computed but `call_ccres` does not use it (original design, under review).
5. `qc.smk`: `ccre_support` tabulates the new anchors; `qc_plots` draws from that table in its
   own env (`envs/qc.yaml`). Plots follow `workflow/scripts/_style.py`: call `apply_style()`
   before importing pyplot, save with `save_figure()`, no per-script rcParams.

`workflow/rules/common.smk` is the only place config and sample sheets are read. It builds the
lookup dicts every rule uses (`NARROWPEAK_OF`, `BIGWIG_OF`, `BIGWIGS[assay][id]`,
`DNASE_BIGWIG_OF`). Fragment rows get their narrowPeak and bigWig paths pointed at
`fragments.smk` outputs here, so downstream rules never branch on input type.

## Conventions

- Layout and rule style follow `~/Projects/snakemake-template` (Nextstrain style guide): no `run:`
  blocks, logic in `workflow/scripts/` argparse scripts with `build_parser()`, every rule has
  `log:`, `benchmark:` and `conda: CONDA_ENV`, shell interpolations quoted with `:q`.
- Every Python rule exports `POLARS_MAX_THREADS={threads}`; profiles set threads via `set-threads`.
- Polars only for tabular work. Use polars-bio for interval ops; do not add bioframe or pandas.
- Outputs, logs and benchmarks are grouped by `run_id` (`results/<run_id>/`, `logs/<run_id>/`,
  `benchmarks/<run_id>/`).
- ENCODE inputs come from the encode-metadata snapshot via
  `workflow/scripts/encode_atac_samples.py`, never a live API pull.
- Scripts reproduce the original pipeline's numbers exactly (e.g. `zscore.py` uses `math.log`, not
  `np.log10`). Keep behavior changes separate from refactors, one logical change per commit.
- Editing `workflow/envs/env.yaml` changes its hash and marks existing results as outdated.
