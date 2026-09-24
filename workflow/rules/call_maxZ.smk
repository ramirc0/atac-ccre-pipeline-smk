rule call_maxZ_ATAC:
    input:
        files=intermediate_bw_files,
    output:
        maxZ=f"{RESULTS_DIR}/{PREFIX}_{GENOME}-ATAC-maxZ.txt",
    params:
        toolkit=TOOLKIT,
    log:
        f"{LOG_DIR}/{PREFIX}_{GENOME}-ATAC-maxZ-ATAC.log"
    benchmark:
        f"{RESULTS_DIR}/benchmarks/call_max-ATAC.tsv"
    run:
        import os, glob
        import numpy as np
        import polars as pl

        files=input.files
        
        os.environ["POLARS_MAX_THREADS"] = "8"
        out_file = output.maxZ

        def read_zscores(path):
            return (
                pl.read_csv(
                    path,
                    separator="\t",
                    has_header=False,
                    columns=[0, 1],
                    new_columns=["anchor", "zscore"],
                    schema_overrides={
                        "anchor": pl.Utf8,
                        "zscore": pl.Utf8,
                    },
                    ignore_errors=True,
                    truncate_ragged_lines=True,
                )
                .with_columns([
                    pl.col("anchor").str.strip_chars(),
                    pl.col("zscore").str.strip_chars().cast(pl.Float64, strict=False),
                ])
            )
        
        first = read_zscores(files[0])
        
        anchor = first["anchor"]
        n = len(anchor)
        
        max_zscore = first["zscore"].to_numpy()
        
        for i, f in enumerate(files[1:], start=2):
            print(f"Processing {i}/{len(files)}: {os.path.basename(f)}")
        
            df = read_zscores(f)
        
            if len(df) != n:
                raise ValueError(f"Row count mismatch in {f}: expected {n}, got {len(df)}")
        
            if not df["anchor"].equals(anchor):
                raise ValueError(f"Anchor order mismatch in {f}")
        
            z = df["zscore"].to_numpy()
        
            mask = np.nan_to_num(z, nan=-np.inf) > np.nan_to_num(max_zscore, nan=-np.inf)
            max_zscore[mask] = z[mask]

        out = pl.DataFrame({
            "anchor": anchor,
            "max_zscore": max_zscore,
        }).with_columns(
            pl.when(pl.col("max_zscore").is_nan())
              .then(None)
              .otherwise(pl.col("max_zscore"))
              .alias("max_zscore")
        )
        out.write_csv(out_file, separator="\t")
        print(f"Done: {out_file}")

rule call_maxZ_DNase:
    input:
        files=intermediate_DNase_files,
    output:
        maxZ=f"{RESULTS_DIR}/{PREFIX}_{GENOME}-DNase-maxZ.txt",
    params:
        toolkit=TOOLKIT,
    log:
        f"{LOG_DIR}/{PREFIX}_{GENOME}-ATAC-maxZ-DNase.log"
    benchmark:
        f"{RESULTS_DIR}/benchmarks/call_maxZ-DNASE.tsv"
    run:
        import os, glob
        import numpy as np
        import polars as pl

        files=input.files
        
        os.environ["POLARS_MAX_THREADS"] = "8"
        out_file = output.maxZ

        def read_zscores(path):
            return (
                pl.read_csv(
                    path,
                    separator="\t",
                    has_header=False,
                    columns=[0, 1],
                    new_columns=["anchor", "zscore"],
                    schema_overrides={
                        "anchor": pl.Utf8,
                        "zscore": pl.Utf8,
                    },
                    ignore_errors=True,
                    truncate_ragged_lines=True,
                )
                .with_columns([
                    pl.col("anchor").str.strip_chars(),
                    pl.col("zscore").str.strip_chars().cast(pl.Float64, strict=False),
                ])
            )
        
        first = read_zscores(files[0])
        
        anchor = first["anchor"]
        n = len(anchor)
        
        max_zscore = first["zscore"].to_numpy()
        
        for i, f in enumerate(files[1:], start=2):
            print(f"Processing {i}/{len(files)}: {os.path.basename(f)}")
        
            df = read_zscores(f)
        
            if len(df) != n:
                raise ValueError(f"Row count mismatch in {f}: expected {n}, got {len(df)}")
        
            if not df["anchor"].equals(anchor):
                raise ValueError(f"Anchor order mismatch in {f}")
        
            z = df["zscore"].to_numpy()
        
            mask = np.nan_to_num(z, nan=-np.inf) > np.nan_to_num(max_zscore, nan=-np.inf)
            max_zscore[mask] = z[mask]

        out = pl.DataFrame({
            "anchor": anchor,
            "max_zscore": max_zscore,
        }).with_columns(
            pl.when(pl.col("max_zscore").is_nan())
              .then(None)
              .otherwise(pl.col("max_zscore"))
              .alias("max_zscore")
        )
        out.write_csv(out_file, separator="\t")
        print(f"Done: {out_file}")