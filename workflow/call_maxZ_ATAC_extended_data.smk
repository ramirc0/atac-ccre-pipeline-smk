rule call_maxZ_ATAC_extended_data:
    input:
        files=intermediate_bw_files,
    output:
        maxZ=f"{RESULTS_DIR}/{PREFIX}_{GENOME}-ATAC-maxZ-summary_extended.txt",
    params:
        toolkit=TOOLKIT,
        output_dir=RESULTS_DIR,
    log:
        f"{LOG_DIR}/{PREFIX}_{GENOME}-ATAC-maxZ-ATAC.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/call_max-ATAC.tsv",
    run:
        import polars as pl
        import glob, os, shutil

        files = [f for f in input.files if os.path.exists(f) and os.path.getsize(f) > 0]

        if len(files) == 0:
            raise FileNotFoundError("No non-empty files found.")

        THRESHOLD = 1.64
        OUT_SUMMARY = output.maxZ
        HITS_DIR = f"{params.output_dir}/hits"
        os.makedirs(HITS_DIR, exist_ok=True)

        # -----------------------------
        # Optional: limit Polars threads
        # -----------------------------
        os.environ["POLARS_MAX_THREADS"] = "4"

        # -----------------------------
        # Initialize from first file
        # -----------------------------
        first_file = files[0]
        first_name = os.path.basename(first_file)
        first_is_custom = "custom" in first_name.lower()

        first = (
            pl.read_csv(
                first_file,
                separator="\t",
                has_header=False,
                columns=[0, 1],
                new_columns=["anchor", "zscore"],
                schema_overrides={"anchor": pl.String, "zscore": pl.String},
                truncate_ragged_lines=True,
                ignore_errors=True,
            )
            .with_columns([
                pl.col("anchor").str.strip_chars(),
                pl.col("zscore").str.strip_chars().cast(pl.Float64, strict=False),
            ])
            .drop_nulls("zscore")
        )

        anchor = first["anchor"]
        n = len(first)

        max_all = first["zscore"]
        src_all = pl.Series([first_name] * n)

        if first_is_custom:
            max_custom = first["zscore"]
            src_custom = pl.Series([first_name] * n)

            max_standard = pl.Series([float("-inf")] * n)
            src_standard = pl.Series([None] * n, dtype=pl.String)

        else:
            max_custom = pl.Series([float("-inf")] * n)
            src_custom = pl.Series([None] * n, dtype=pl.String)

            max_standard = first["zscore"]
            src_standard = pl.Series([first_name] * n)

        # -----------------------------
        # Save threshold hits
        # -----------------------------
        hits = (
            first.filter(pl.col("zscore") > THRESHOLD)
            .select("anchor")
            .with_columns(pl.lit(first_name).alias("filename"))
        )

        if hits.height > 0:
            hits.write_parquet(os.path.join(HITS_DIR, f"{first_name}.parquet"))

        # -----------------------------
        # Process remaining files
        # -----------------------------
        for f in files[1:]:
            fname = os.path.basename(f)
            is_custom = "custom" in fname.lower()

            print(f"Processing: {fname}")

            df = (
                pl.read_csv(
                    f,
                    separator="\t",
                    has_header=False,
                    columns=[0, 1],
                    new_columns=["anchor", "zscore"],
                    schema_overrides={"anchor": pl.String, "zscore": pl.String},
                    truncate_ragged_lines=True,
                    ignore_errors=True,
                )
                .with_columns([
                    pl.col("anchor").str.strip_chars(),
                    pl.col("zscore").str.strip_chars().cast(pl.Float64, strict=False),
                ])
                .drop_nulls("zscore")
            )

            if len(df) != n:
                raise ValueError(f"{fname} has {len(df)} rows, expected {n}.")

            z = df["zscore"]

            # -----------------------------
            # Save threshold hits
            # -----------------------------
            hits = (
                df.filter(pl.col("zscore") > THRESHOLD)
                .select("anchor")
                .with_columns(pl.lit(fname).alias("filename"))
            )

            if hits.height > 0:
                hits.write_parquet(os.path.join(HITS_DIR, f"{fname}.parquet"))

            # -----------------------------
            # Overall max
            # -----------------------------
            better = z > max_all
            max_all = pl.when(better).then(z).otherwise(max_all)
            src_all = pl.when(better).then(pl.lit(fname)).otherwise(src_all)

            # -----------------------------
            # Custom or standard max
            # -----------------------------
            if is_custom:
                better = z > max_custom
                max_custom = pl.when(better).then(z).otherwise(max_custom)
                src_custom = pl.when(better).then(pl.lit(fname)).otherwise(src_custom)

            else:
                better = z > max_standard
                max_standard = pl.when(better).then(z).otherwise(max_standard)
                src_standard = pl.when(better).then(pl.lit(fname)).otherwise(src_standard)

        # -----------------------------
        # Reconstruct per-anchor file list
        # -----------------------------
        parquet_files = glob.glob(os.path.join(HITS_DIR, "*.parquet"))

        if parquet_files:
            per_anchor = (
                pl.scan_parquet(os.path.join(HITS_DIR, "*.parquet"))
                .group_by("anchor")
                .agg(
                    pl.col("filename")
                    .sort()
                    .str.join(";")
                    .alias("files_above_1.64")
                )
                .collect()
            )
        else:
            per_anchor = pl.DataFrame(
                {"anchor": [], "files_above_1.64": []},
                schema={"anchor": pl.String, "files_above_1.64": pl.String},
            )

        # -----------------------------
        # Final output
        # -----------------------------
        out = pl.DataFrame({
            "anchor": anchor,
            "max_all": max_all,
            "max_all_file": src_all,
            "max_custom": max_custom,
            "max_custom_file": src_custom,
            "max_standard": max_standard,
            "max_standard_file": src_standard,
        })

        out = out.with_columns([
            pl.when(pl.col("max_custom") == float("-inf"))
            .then(None)
            .otherwise(pl.col("max_custom"))
            .alias("max_custom"),

            pl.when(pl.col("max_standard") == float("-inf"))
            .then(None)
            .otherwise(pl.col("max_standard"))
            .alias("max_standard"),

            pl.col("max_custom_file").fill_null("."),
            pl.col("max_standard_file").fill_null("."),
        ])

        out = out.join(per_anchor, on="anchor", how="left")

        out = out.with_columns(
            pl.col("files_above_1.64").fill_null("")
        )

        out.write_csv(OUT_SUMMARY, separator="\t")

        shutil.rmtree(HITS_DIR)

        print("Done.")
        print(f"Summary: {OUT_SUMMARY}")