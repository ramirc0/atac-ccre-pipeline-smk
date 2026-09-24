# QC tables and plots over the finished outputs.


rule ccre_support:
    input:
        anchors=f"{OUTDIR}/Anchors-ATAC.bed",
        summary=f"{OUTDIR}/ATAC-Summary.txt",
        ccres=f"{OUTDIR}/cCREs.bed",
        atac_max_z=f"{OUTDIR}/ATAC-maxZ.txt",
        dnase_max_z=[f"{OUTDIR}/DNase-maxZ.txt"] if DNASE_BIGWIG_OF else [],
        zscores=expand(f"{OUTDIR}/zscore/ATAC/{{id}}.parquet", id=BIGWIG_OF),
    output:
        f"{OUTDIR}/qc/ccre-support.tsv",
    params:
        genome=GENOME,
    log:
        f"{LOGDIR}/ccre_support.txt",
    benchmark:
        f"{BENCHDIR}/ccre_support.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})
        export POLARS_MAX_THREADS={threads}

        python workflow/scripts/ccre_support.py \
            --anchors {input.anchors:q} \
            --summary {input.summary:q} \
            --ccres {input.ccres:q} \
            --atac-max-z {input.atac_max_z:q} \
            --dnase-max-z {input.dnase_max_z:q} \
            --zscores {input.zscores:q} \
            --genome {params.genome:q} \
            --output {output:q}
        """
