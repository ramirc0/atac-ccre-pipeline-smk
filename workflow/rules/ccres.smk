# Max z-scores and cCRE calls.


rule call_maxZ_ATAC:
    input:
        files=intermediate_bw_files,
    output:
        maxZ=f"{RESULTS_DIR}/{PREFIX}_{GENOME}-ATAC-maxZ.txt",
    log:
        f"{LOG_DIR}/{PREFIX}_{GENOME}-ATAC-maxZ-ATAC.log"
    benchmark:
        f"{RESULTS_DIR}/benchmarks/call_max-ATAC.tsv"
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/max_z.py \
            --inputs {input.files:q} \
            --output {output.maxZ:q}
        """


rule call_maxZ_DNase:
    input:
        files=intermediate_DNase_files,
    output:
        maxZ=f"{RESULTS_DIR}/{PREFIX}_{GENOME}-DNase-maxZ.txt",
    log:
        f"{LOG_DIR}/{PREFIX}_{GENOME}-ATAC-maxZ-DNase.log"
    benchmark:
        f"{RESULTS_DIR}/benchmarks/call_maxZ-DNASE.tsv"
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/max_z.py \
            --inputs {input.files:q} \
            --output {output.maxZ:q}
        """


rule call_ATAC_cCREs:
    input:
        ATAC_maxZ=f"{RESULTS_DIR}/{PREFIX}_{GENOME}-ATAC-maxZ.txt",
        ANCHORS=f"{RESULTS_DIR}/{PREFIX}_{GENOME}-Anchors-ATAC.bed",
        cCRES=config["ccre_path"],
    output:
        cCRES=f"{RESULTS_DIR}/{PREFIX}_{GENOME}-cCREs.bed",
    params:
        toolkit=TOOLKIT,
        genome=GENOME,
    log:
        f"{LOG_DIR}/{PREFIX}_{GENOME}-cCREs.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/call_cCREs_with_ATAC.tsv",
    run:
        import pandas as pd
        import bioframe

        def Set_IDs(genome, mode):
            genomeDict = {
                "hg19": "EH37",
                "hg38": "EH38",
                "mm10": "EM10",
            }

            modeDict = {
                "rDHS": "D",
                "ccRE": "E",
                "cCRE": "E",
                "TF": "F",
                "MultiMap": "M",
                "ATAC": "A",
                "BETA": "B",
            }

            if genome not in genomeDict:
                raise ValueError(f"Invalid genome: {genome}")

            if mode not in modeDict:
                raise ValueError(f"Invalid mode: {mode}")

            return genomeDict[genome], modeDict[mode]

        def Process_Previous(previousRegions, modeID, id_col="ccre"):
            previousDict = {}
            maxAccession = 0

            for row in previousRegions.itertuples(index=False):
                key = f"{row.chrom}{row.start}{row.end}"
                accession = getattr(row, id_col)

                previousDict[key] = accession

                number = int(str(accession).split(modeID)[-1])
                maxAccession = max(maxAccession, number)

            return previousDict, maxAccession

        def assign_accessions(regions, previousDict, genomeID, modeID, start_i):
            i = start_i
            rows = []

            for row in regions.itertuples(index=False):
                coords = f"{row.chrom}{row.start}{row.end}"

                if coords in previousDict:
                    accession = previousDict[coords].replace("EH37", "EH38")
                else:
                    accession = genomeID + modeID + str(i).zfill(7)
                    i += 1

                rows.append([
                    row.chrom,
                    row.start,
                    row.end,
                    row.anchor,
                    accession,
                    row.type,
                ])

            return pd.DataFrame(
                rows,
                columns=["chrom", "start", "end", "anchor", "ccre", "type"],
            )

        genome = params.genome
        mode = "ccRE"

        genomeID, modeID = Set_IDs(genome, mode)

        current_ccres = pd.read_csv(
            input.cCRES,
            sep="\t",
            header=None,
            usecols=[0, 1, 2, 3, 4, 5],
            names=["chrom", "start", "end", "anchor", "ccre", "type"],
        )

        previousRegions = current_ccres[["chrom", "start", "end", "ccre", "type"]].copy()
        previousRegions["version"] = "V4"

        previousDict, maxAccession = Process_Previous(previousRegions, modeID)
        start_i = maxAccession + 1

        ATAC_maxZ = pd.read_csv(input.ATAC_maxZ, sep="\t")

        ATAC_maxZ = ATAC_maxZ[
            (ATAC_maxZ["max_zscore"] > 1.64)
            & (ATAC_maxZ["anchor"].str.contains("EH38A", na=False))
        ].copy()

        ANCHORS = pd.read_csv(
            input.ANCHORS,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "anchor"],
        )

        new_cCREs = (
            pd.merge(ATAC_maxZ, ANCHORS, on="anchor")
            [["chrom", "start", "end", "anchor"]]
            .copy()
        )

        new_cCREs["type"] = "ATAC-cCRE"

        regions = pd.concat(
            [
                new_cCREs,
                current_ccres[["chrom", "start", "end", "anchor", "type"]],
            ],
            ignore_index=True,
        )

        out = assign_accessions(
            regions=regions,
            previousDict=previousDict,
            genomeID=genomeID,
            modeID=modeID,
            start_i=start_i,
        )

        out["chrom"] = out["chrom"].astype(str)
        out = bioframe.sort_bedframe(out)

        out.to_csv(output.cCRES, sep="\t", header=True, index=False)
