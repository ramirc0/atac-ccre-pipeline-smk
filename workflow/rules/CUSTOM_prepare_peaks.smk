rule CUSTOM_prepare_peaks:
    input:
        custom_info=config["CUSTOM_data_macs_info"]
    output:
        tmp_bed=f"{config['ouput_dir']}/CUSTOM_{config['genome']}-ATAC-tmp.bed",
        sorted_bed=f"{config['ouput_dir']}/CUSTOM_{config['genome']}-ATAC-tmp.sorted.bed"
    log:
        f"{config['log_dir']}/CUSTOM_prepare_peaks.log"
    benchmark:
        f"{RESULTS_DIR}/benchmarks/CUSTOM_prepare_peaks.tsv"
    run:
        import os
        import glob
        import subprocess
        import polars as pl

        data_type = "ATAC"
        genome = config["genome"]

        high_qc_files = (
            pl.read_csv(
                input.custom_info,
                separator="\t",
                has_header=True,
                schema_overrides={
                    "sample": pl.Utf8,
                    "path": pl.Utf8,
                },
            )
            .select(["sample", "path"])
            .rows()
        )

        _datalocal = True
        out_file = output.tmp_bed
        sorted_file = output.sorted_bed

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

        # Remove previous outputs
        for f in [out_file, sorted_file]:
            if os.path.exists(f):
                os.remove(f)

        with open(out_file, "w") as out_f:
            for sample, file_path in high_qc_files:
                exp = sample
                target = "ATAC"
                biosample = sample

                bed_gz = None

                if _datalocal:
                    matches = glob.glob(file_path)
                    if matches:
                        bed_gz = matches[0]
                    else:
                        print(f"WARNING: No file found for {sample}: {file_path}")
                        continue
                else:
                    bed_gz = file_path

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

                width = pl.col("end") - pl.col("start")
                summit = pl.col("start") + pl.col("peak")

                df = (
                    df.with_row_index("row_nr", offset=1)
                    .filter(
                        pl.col("chrom").str.contains(
                            r"^chr([1-9]|1[0-9]|2[0-2]|X|Y)$"
                        )
                    )
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
                        (
                            pl.lit(f"{exp}-")
                            + pl.col("row_nr").cast(pl.Utf8)
                        ).alias("unique_id"),
                    ])
                    .drop("row_nr")
                    .filter(pl.col("start") < pl.col("end"))
                )

                df.write_csv(
                    out_f,
                    separator="\t",
                    include_header=False,
                )

        # Final genome sort
        subprocess.run(
            f"sort -k1,1V -k2,2n {out_file} > {sorted_file}",
            shell=True,
            check=True,
        )

        print(f"Finished writing: {out_file}")
        print(f"Finished sorting: {sorted_file}")