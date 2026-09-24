# Max z-scores and cCRE calls.


rule max_z:
    input:
        lambda w: expand(f"{OUTDIR}/zscore/{w.assay}/{{id}}.parquet", id=BIGWIGS[w.assay]),
    output:
        f"{OUTDIR}/{{assay}}-maxZ.txt",
    log:
        f"{LOGDIR}/max_z/{{assay}}.txt",
    benchmark:
        f"{BENCHDIR}/max_z/{{assay}}.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})
        export POLARS_MAX_THREADS={threads}

        python workflow/scripts/max_z.py \
            --inputs {input:q} \
            --output {output:q}
        """


rule call_ccres:
    input:
        max_z=f"{OUTDIR}/ATAC-maxZ.txt",
        anchors=f"{OUTDIR}/Anchors-ATAC.bed",
        ccres=REFS["ccres"],
    output:
        f"{OUTDIR}/cCREs.bed",
    params:
        genome=GENOME,
    log:
        f"{LOGDIR}/call_ccres.txt",
    benchmark:
        f"{BENCHDIR}/call_ccres.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})
        export POLARS_MAX_THREADS={threads}

        python workflow/scripts/call_ccres.py \
            --max-z {input.max_z:q} \
            --anchors {input.anchors:q} \
            --ccres {input.ccres:q} \
            --genome {params.genome:q} \
            --output {output:q}
        """
