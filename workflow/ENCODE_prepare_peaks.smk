rule ENCODE_prepare_ATAC_peaks:
    input:
        filtered_list = f"{RESULTS_DIR}/{GENOME}-ATAC-List.Filtered.txt"
    output:
        tmp_bed = f"{RESULTS_DIR}/ENCODE_{GENOME}-ATAC-tmp.bed",
        sorted_bed = f"{RESULTS_DIR}/ENCODE_{GENOME}-ATAC-tmp.sorted.bed"
    params:
        genome = GENOME,
        data_type = "ATAC",
        datalocal = True,
        encode_data_dir = "/data/projects/encode/data"
    benchmark:
        f"{RESULTS_DIR}/benchmarks/ENCODE_prepare_ATAC_peaks.tsv"
    log:
        f"{LOG_DIR}/ENCODE_prepare_peaks.log"
    run:
        import os
        import glob
        import subprocess
        from pathlib import Path
        import polars as pl

        input_list = input.filtered_list
        out_file = output.tmp_bed
        sorted_file = output.sorted_bed
        log_file = log[0]

        Path(out_file).parent.mkdir(parents=True, exist_ok=True)
        Path(sorted_file).parent.mkdir(parents=True, exist_ok=True)
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)

        peak_df = (
            pl.read_csv(
                input_list,
                separator="\t",
                has_header=False,
                columns=[0, 1, 3],
                new_columns=["exp", "peak", "biosample"],
            )
            .with_columns(pl.lit("ATAC").alias("target"))
            .select(["exp", "peak", "biosample", "target"])
        )

        schema = {
            "chrom": pl.Utf8,
            "start": pl.Int64,
            "end": pl.Int64,
            "name": pl.Utf8,
            "score": pl.Float64,
            "strand": pl.Utf8,
            "signalValue": pl.Float64,
            "pValue": pl.Float64,
            "qValue": pl.Float64,
            "peak": pl.Int64,
        }

        for f in [out_file, sorted_file]:
            if os.path.exists(f):
                os.remove(f)

        n = peak_df.height

        with open(log_file, "w") as log_f, open(out_file, "w") as out_f:
            for x, row in enumerate(peak_df.iter_rows(named=True), start=1):
                exp = row["exp"]
                peak = row["peak"]
                biosample = row["biosample"]
                target = row["target"]

                matches = glob.glob(f"{params.encode_data_dir}/{exp}/{peak}.bed.gz")
                if not matches:
                    print(f"Missing peak file: {exp}/{peak}.bed.gz", file=log_f, flush=True)
                    continue

                bed_gz = matches[0]

                try:
                    df = pl.read_csv(
                        bed_gz,
                        separator="\t",
                        has_header=False,
                        new_columns=[
                            "chrom", "start", "end", "name", "score", "strand",
                            "signalValue", "pValue", "qValue", "peak"
                        ],
                        schema_overrides=schema,
                    )
                except Exception as e:
                    print(f"Failed reading {bed_gz}: {e}", file=log_f, flush=True)
                    continue

                width = pl.col("end") - pl.col("start")
                summit = pl.col("start") + pl.col("peak")

                df = (
                    df.with_row_index("row_nr", offset=1)
                    .filter(pl.col("chrom").str.contains(r"^chr([1-9]|1[0-9]|2[0-2]|X|Y)$"))
                    .with_columns([
                        pl.when(width < 150)
                        .then(summit - 75)
                        .when(width > 350)
                        .then(summit - 175)
                        .otherwise(pl.col("start"))
                        .alias("start"),

                        pl.when(width < 150)
                        .then(summit + 75)
                        .when(width > 350)
                        .then(summit + 175)
                        .otherwise(pl.col("end"))
                        .alias("end"),

                        pl.lit(exp).alias("exp"),
                        pl.lit(target).alias("target"),
                        pl.lit(biosample).alias("biosample"),
                        (pl.lit(f"{exp}-") + pl.col("row_nr").cast(pl.Utf8)).alias("unique_id"),
                    ])
                    .drop("row_nr")
                    .filter(pl.col("start") < pl.col("end"))
                )

                df.write_csv(out_f, separator="\t", include_header=False)

                if x % 500 == 0:
                    print(f"done with {x}/{n} ({100*x/n:.1f}%)", file=log_f, flush=True)

        shell("sort -k1,1V -k2,2n {output.tmp_bed} > {output.sorted_bed}")