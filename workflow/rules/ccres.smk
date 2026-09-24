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
        genome=GENOME,
    log:
        f"{LOG_DIR}/{PREFIX}_{GENOME}-cCREs.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/call_cCREs_with_ATAC.tsv",
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/call_ccres.py \
            --max-z {input.ATAC_maxZ:q} \
            --anchors {input.ANCHORS:q} \
            --ccres {input.cCRES:q} \
            --genome {params.genome:q} \
            --output {output.cCRES:q}
        """
