# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Snakemake pipeline that calls ATAC cCREs on top of the ENCODE Registry V7 (hg38). Read `README.md`
for setup, sample sheets, reference files and preflight scripts, and `docs/inputs.md` for config
keys and output formats.

## Commands

The launcher is the root pixi workspace (Snakemake main, pixi and SLURM plugins). Rule tools come
from the pixi workspace in `workflow/envs/` (`default` env, `qc` for `qc_plots`). Run the commands
below inside `pixi shell`.

```sh
snakemake -n -p --profile profiles/local        # dry run
snakemake --profile profiles/slurm              # full run on SLURM
snakemake --profile profiles/local results/<run_id>/ATAC.bed   # one target
snakemake --lint                                # flags every rule for lacking conda next to pixi
pixi run lint                                   # ruff check, ruff format --check, numpydoc lint
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
   own pixi env (`qc`). Plots use [poikilos](https://github.com/ramirc0/poikilos) with the
   `plain` theme: `pk.use()` in `main()`, with any rcParam change passed as `rc=`.
   `pk.despine()` goes on every axes. `pk.save_figure()` writes the SVG + PNG pair.

`workflow/rules/common.smk` is the only place config and sample sheets are read. It builds the
lookup dicts every rule uses (`NARROWPEAK_OF`, `BIGWIG_OF`, `BIGWIGS[assay][id]`,
`DNASE_BIGWIG_OF`). Fragment rows get their narrowPeak and bigWig paths pointed at
`fragments.smk` outputs here, so downstream rules never branch on input type.

## Conventions

- Layout and rule style follow `~/Projects/snakemake-template` (Nextstrain style guide): no `run:`
  blocks, logic in `workflow/scripts/` argparse scripts with `build_parser()`, every rule has
  `log:`, `benchmark:` and `software: SOFTWARE_ENV` (`QC_ENV` for plots), shell interpolations
  quoted with `:q`.
- Script docstrings follow numpydoc and pass `pixi run lint`. The linters skip `.smk` files, so keep
  the `common.smk` helpers in the same style by hand.
- Every Python rule exports `POLARS_MAX_THREADS={threads}`; profiles set threads via `set-threads`.
- Polars only for tabular work. Use polars-bio for interval ops; do not add bioframe or pandas.
- Outputs, logs and benchmarks are grouped by `run_id` (`results/<run_id>/`, `logs/<run_id>/`,
  `benchmarks/<run_id>/`).
- ENCODE inputs come from the encode-metadata snapshot via
  `workflow/scripts/encode_atac_samples.py`, never a live API pull.
- Scripts reproduce the original pipeline's numbers exactly (e.g. `zscore.py` uses `math.log`, not
  `np.log10`). Keep behavior changes separate from refactors, one logical change per commit.
- Editing `workflow/envs/pixi.toml` or `pixi.lock` changes the env hash and marks existing results
  as outdated. So does moving the workdir, because the hash includes the workspace's absolute path.
