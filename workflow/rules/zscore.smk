# Per-bigWig z-scores over anchors.


# Anchor regions scored by every bigWig: deduplicated, widened by `width`.
rule anchor_regions:
    input:
        f"{RESULTS_DIR}/{PREFIX}_{GENOME}-Anchors-ATAC.bed",
    output:
        f"{RESULTS_DIR}/{PREFIX}_{GENOME}-Anchors-ATAC.regions.bed",
    params:
        width=0,
    log:
        f"{LOG_DIR}/{PREFIX}_anchor_regions.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/{PREFIX}_anchor_regions.tsv",
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


rule CUSTOM_Zscore_across_bw:
    input:
        regions=f"{RESULTS_DIR}/{PREFIX}_{GENOME}-Anchors-ATAC.regions.bed",
        bigwig=lambda wc: CUSTOM_BIGWIGS[wc.sample],
    output:
        zscore=f"{RESULTS_DIR}/bw_zscoring/custom-{{sample}}.parquet",
    params:
        bigwigAverageOverBed="/zata/data/zlab/common/tools/ucsc.v385/bigWigAverageOverBed",
    log:
        f"{LOG_DIR}/bw_zscoring/custom-{{sample}}.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {output.zscore})" "$(dirname {log})"

        workdir=$(mktemp -d "{TMP_DIR}/custom_bw_zscore_{wildcards.sample}.XXXXXX")
        trap 'rm -rf "$workdir"' EXIT

        echo "Processing custom sample: {wildcards.sample}" > {log}

        {params.bigwigAverageOverBed} \
            {input.bigwig:q} \
            {input.regions:q} \
            "$workdir/out2" >> {log} 2>&1

        python workflow/scripts/zscore.py \
            --input "$workdir/out2" \
            --output {output.zscore:q}

        """


rule ENCODE_Zscore_across_bw:
    input:
        regions=f"{RESULTS_DIR}/{PREFIX}_{GENOME}-Anchors-ATAC.regions.bed",
        bigwig=f"{ENCODE_DATA_DIR}/{{experiment}}/{{file}}.bigWig",
    output:
        zscore=f"{RESULTS_DIR}/bw_zscoring/encode-{{experiment}}_{{file}}.parquet",
    params:
        bigwigAverageOverBed="/zata/data/zlab/common/tools/ucsc.v385/bigWigAverageOverBed",
    log:
        f"{LOG_DIR}/bw_zscoring/encode-{{experiment}}_{{file}}.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {output.zscore})" "$(dirname {log})"

        workdir=$(mktemp -d "{TMP_DIR}/bw_zscore_{wildcards.experiment}_{wildcards.file}.XXXXXX")
        trap 'rm -rf "$workdir"' EXIT

        echo "Processing {wildcards.experiment} / {wildcards.file}" > {log}

        {params.bigwigAverageOverBed} \
            {input.bigwig:q} \
            {input.regions:q} \
            "$workdir/out2" >> {log} 2>&1

        python workflow/scripts/zscore.py \
            --input "$workdir/out2" \
            --output {output.zscore:q}
        """


rule ENCODE_Zscore_across_DNase:
    input:
        regions=f"{RESULTS_DIR}/{PREFIX}_{GENOME}-Anchors-ATAC.regions.bed",
        bigwig=f"{ENCODE_DATA_DIR}/{{experiment}}/{{file}}.bigWig",
    output:
        zscore=f"{RESULTS_DIR}/bw_zscoring-DNAse/encode-{{experiment}}_{{file}}.parquet",
    params:
        bigwigAverageOverBed="/zata/data/zlab/common/tools/ucsc.v385/bigWigAverageOverBed",
    log:
        f"{LOG_DIR}/bw_zscoring-DNAse/encode-{{experiment}}_{{file}}.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {output.zscore})" "$(dirname {log})"

        workdir=$(mktemp -d "{TMP_DIR}/bw_zscore_{wildcards.experiment}_{wildcards.file}.XXXXXX")
        trap 'rm -rf "$workdir"' EXIT

        echo "Processing {wildcards.experiment} / {wildcards.file}" > {log}

        {params.bigwigAverageOverBed} \
            {input.bigwig:q} \
            {input.regions:q} \
            "$workdir/out2" >> {log} 2>&1

        python workflow/scripts/zscore.py \
            --input "$workdir/out2" \
            --output {output.zscore:q}
        """
