rule plot_extended_data:
    input:
        extended_summary=f"{RESULTS_DIR}/{PREFIX}_{GENOME}-ATAC-maxZ-summary_extended.txt",
        ANCHORS=f"{RESULTS_DIR}/{PREFIX}_{GENOME}-Anchors-ATAC.bed",
        ATAC_SUMMARY=f"{RESULTS_DIR}/{PREFIX}_{GENOME}-ATAC-Summary.txt",
    output:
        done=f"{RESULTS_DIR}/plots/plot_extended_data.done",
    params:
        toolkit=TOOLKIT,
        output_dir=f"{RESULTS_DIR}/plots",
        buffer=10000,
        top_n=50,
    log:
        f"{LOG_DIR}/{PREFIX}_{GENOME}-plot_extended_data.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/plot_extended_data.tsv",
    run:
        import os
        import re
        import shlex
        import pandas as pd
        import subprocess
        from pathlib import Path
        from io import StringIO

        os.makedirs(params.output_dir, exist_ok=True)
        os.makedirs("tmp", exist_ok=True)


        def extract_dnase_hits_df(
            input_anchor_id,
            index_file="index_cCRE.txt",
            matrix_file="/data/projects/encode/Registry/V4/GRCh38/Matrices/GRCh38.DNase-zscore.cCRE-V4.txt",
            threshold=1.64,
        ):
            grep_cmd = (
                f"grep -n -w {shlex.quote(str(input_anchor_id))} "
                f"{shlex.quote(str(index_file))} | cut -d: -f1"
            )

            index_result = subprocess.check_output(
                grep_cmd,
                shell=True,
                text=True
            ).strip()

            if not index_result:
                return pd.DataFrame(columns=["indx", "expID", "zscore"])

            matrix_index = index_result.splitlines()[0]

            awk_cmd = rf"""
            awk -F'\t' -v idx="{matrix_index}" -v thresh="{threshold}" '
            NR==1 {{
                for (i=1; i<=NF; i++) header[i]=$i
                next
            }}
            NR==idx {{
                for (i=2; i<=NF; i++) {{
                    if ($i >= thresh) {{
                        print i "\t" header[i] "\t" $i
                    }}
                }}
                exit
            }}
            ' {shlex.quote(str(matrix_file))}
            """

            output_txt = subprocess.check_output(
                awk_cmd,
                shell=True,
                text=True
            )

            if not output_txt.strip():
                return pd.DataFrame(columns=["indx", "expID", "zscore"])

            return pd.read_csv(
                StringIO(output_txt),
                sep="\t",
                header=None,
                names=["indx", "expID", "zscore"]
            )


        def get_dnase(anchor_value, atac_df, search_window=5, min_dist=10000, top_n=50):
            dictionary = pd.read_csv(
                "/data/projects/encode/Registry/V4/GRCh38/Biosample-Lists/DNase-List.txt",
                sep="\t",
                header=None,
                skiprows=1,
                usecols=[0, 1, 2],
                names=["expID", "fileID", "sample"],
            )

            atac_df = atac_df.copy()
            atac_df["anchor_order"] = range(len(atac_df))

            match = atac_df.loc[atac_df["anchor"] == anchor_value]

            if match.empty:
                return []

            row = match.iloc[0]

            idx_toget = int(row["anchor_order"])
            start = int(row["start"])
            end = int(row["end"])
            chrom = row["chrom"]

            nearby = atac_df.iloc[
                max(0, idx_toget - search_window): idx_toget + search_window + 1
            ]

            nearby_dnase_anchors = []

            for _, r in nearby.iterrows():
                anchor = str(r["anchor"])

                if "EH38A" in anchor:
                    continue

                if r["chrom"] != chrom:
                    continue

                r_start = int(r["start"])
                r_end = int(r["end"])

                if r_start <= end + min_dist and r_end >= start - min_dist:
                    nearby_dnase_anchors.append(anchor)

            signal_dir = Path("/zata/data/zlab/projects/encode/data")
            unique_signal_files = []
            seen = set()

            for anchor_lookup in nearby_dnase_anchors:
                dnase = extract_dnase_hits_df(anchor_lookup)

                if dnase.empty:
                    continue

                dnase = dnase.sort_values("zscore", ascending=False).head(top_n)
                dnase_to_pull = pd.merge(dnase, dictionary, on="expID", how="inner")

                for r in dnase_to_pull.itertuples(index=False):
                    file_path = str(signal_dir / f"{r.expID}/{r.fileID}.bigWig")
                    key = (file_path, r.sample)

                    if key not in seen:
                        seen.add(key)
                        unique_signal_files.append([file_path, r.sample])

            return unique_signal_files


        def extract_ENCODE_ATAC_overlap(files_above, ATAC_List):
            if pd.isna(files_above) or str(files_above).strip() == "":
                return []

            files = []
            seen = set()

            for name in str(files_above).split(";"):
                name = name.strip()

                if not name or "-" not in name:
                    continue

                bio, rest = name.split("-", 1)
                sample = rest.split(".")[0]

                possible_paths = [
                    f"/zata/data/zlab/projects/encode/data/{bio}/{sample}.bigWig",
                    f"/data/projects/encode/data/{bio}/{sample}.bigWig",
                ]

                file_path = next((p for p in possible_paths if os.path.exists(p)), None)

                if file_path is None:
                    continue

                sample_rows = ATAC_List.loc[ATAC_List["IDs"] == bio, "name"]
                sample_name = sample_rows.iloc[0] if len(sample_rows) > 0 else bio

                key = (file_path, sample_name)

                if key not in seen:
                    seen.add(key)
                    files.append([file_path, sample_name])

            return files


        def plot_genome_atac_only(
            x,
            anchor_id,
            buffer,
            output_loc,
            counter,
            extended_info,
            unique_signal_files=None,
            unique_signal_files_atac=None,
        ):
            final_region = (
                f"{x.iloc[0]['chrom']}:"
                f"{max(int(x['start'].min()) - buffer, 0)}-"
                f"{int(x['end'].max()) + buffer}"
            )

            safe_final_region = final_region.replace(":", "_").replace("-", "_")
            safe_info = re.sub(r"[^a-zA-Z0-9.,=_-]", "", str(extended_info))

            out_png = (
                f"{output_loc}/"
                f"{counter}_{anchor_id}_{safe_final_region}_{safe_info}.png"
            )

            # Add your track writing logic here
            # pyGenomeTracks command:
            cmd = [
                "pyGenomeTracks",
                "--tracks", "track.ini",
                "--region", final_region,
                "--trackLabelFraction", "0.15",
                "--fontSize", "8",
                "--outFileName", out_png,
            ]

            ret = subprocess.run(cmd)

            if ret.returncode != 0:
                raise RuntimeError(f"pyGenomeTracks failed for {anchor_id}")

            return out_png


        # -----------------------------
        # Main execution
        # -----------------------------
        output_df = pd.read_csv(input.extended_summary, sep="\t")

        atac_as_TF_v2 = pd.read_csv(
            input.ANCHORS,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "anchor"]
        )

        _ATAC_List = pd.read_csv(
            "/zata/data/zlab/projects/encode/Registry/V4/GRCh38/Biosample-Lists/ATAC-List.txt",
            sep="\t",
            header=None,
            names=["IDs", "file", "name"]
        )

        subset = (
            output_df.loc[
                (output_df["max_custom"] >= 1.64) &
                (output_df["anchor"].str.contains("EH38A"))
            ]
            .sort_values("max_all", ascending=False)
            .reset_index(drop=True)
        )

        for num, example in subset.iterrows():

            val = example["anchor"]

            extended_info = (
                f"BETA={round(example['max_custom'],2)},"
                f"_ENCODE_ONLY={round(example['max_standard'],2)},"
                f"_ALL={round(example['max_all'],2)}"
            )

            unique_dnase_signal_files = get_dnase(
                val,
                atac_as_TF_v2,
                min_dist=params.buffer,
                top_n=params.top_n
            )

            unique_atac_signal = extract_ENCODE_ATAC_overlap(
                example.get("files_above_1.64", ""),
                _ATAC_List
            )

            plot_genome_atac_only(
                atac_as_TF_v2[atac_as_TF_v2["anchor"] == val],
                anchor_id=val,
                buffer=params.buffer,
                output_loc=params.output_dir,
                counter=num,
                extended_info=extended_info,
                unique_signal_files=unique_dnase_signal_files,
                unique_signal_files_atac=unique_atac_signal,
            )

        Path(output.done).touch()