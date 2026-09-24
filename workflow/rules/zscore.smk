# Per-bigWig z-scores over anchors.


# Anchor regions scored by every bigWig: deduplicated, widened by `width`.
rule anchor_regions:
    input:
        f"{OUTDIR}/Anchors-ATAC.bed",
    output:
        f"{OUTDIR}/zscore/anchor-regions.bed",
    params:
        width=0,
    log:
        f"{LOGDIR}/anchor_regions.txt",
    benchmark:
        f"{BENCHDIR}/anchor_regions.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        awk -v width={params.width:q} '
        BEGIN {{ FS=OFS="\t" }}
        {{
            start = $2 - width
            end   = $3 + width
            if (start < 0) start = 0
            print $1, start, end, $4
        }}
        ' {input:q} | sort -u > {output:q}
        """


# One job per bigWig: ATAC samples and DNase experiments alike.
rule zscore:
    input:
        regions=f"{OUTDIR}/zscore/anchor-regions.bed",
        bigwig=lambda w: BIGWIGS[w.assay][w.id],
    output:
        f"{OUTDIR}/zscore/{{assay}}/{{id}}.parquet",
    log:
        f"{LOGDIR}/zscore/{{assay}}/{{id}}.txt",
    benchmark:
        f"{BENCHDIR}/zscore/{{assay}}/{{id}}.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        workdir=$(mktemp -d)
        trap 'rm -rf "$workdir"' EXIT

        bigWigAverageOverBed \
            {input.bigwig:q} \
            {input.regions:q} \
            "$workdir/signal.tab"

        python workflow/scripts/zscore.py \
            --input "$workdir/signal.tab" \
            --output {output:q}
        """
