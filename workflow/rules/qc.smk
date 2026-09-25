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


# Page 3 needs DNase; keep in step with qc_plots.py.
QC_PAGES = [
    "1-funnel",
    "2-atac-max-z",
    *(["3-atac-vs-dnase"] if DNASE_BIGWIG_OF else []),
    "4-support",
    "5-nearest-registry",
    "6-sample-contribution",
]


rule qc_plots:
    input:
        support=f"{OUTDIR}/qc/ccre-support.tsv",
        atac_max_z=f"{OUTDIR}/ATAC-maxZ.txt",
        peaks=f"{OUTDIR}/peaks/peaks.bed",
        rpeaks=f"{OUTDIR}/rpeaks/rpeaks.bed",
        filtered=f"{OUTDIR}/rpeaks/filtered.tsv",
        no_rdhs=f"{OUTDIR}/rpeaks/no-rdhs.tsv",
    output:
        pdf=f"{OUTDIR}/qc/qc.pdf",
        figures=expand(f"{OUTDIR}/qc/figures/{{page}}.{{ext}}", page=QC_PAGES, ext=["svg", "png"]),
    params:
        figures=f"{OUTDIR}/qc/figures",
    log:
        f"{LOGDIR}/qc_plots.txt",
    benchmark:
        f"{BENCHDIR}/qc_plots.tsv"
    conda:
        "../envs/qc.yaml"
    shell:
        r"""
        exec &> >(tee {log:q})
        export POLARS_MAX_THREADS={threads}

        python workflow/scripts/qc_plots.py \
            --support {input.support:q} \
            --atac-max-z {input.atac_max_z:q} \
            --peaks {input.peaks:q} \
            --rpeaks {input.rpeaks:q} \
            --filtered {input.filtered:q} \
            --no-rdhs {input.no_rdhs:q} \
            --figures {params.figures:q} \
            --output {output.pdf:q}
        """
