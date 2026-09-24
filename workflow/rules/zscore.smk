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
    output:
        zscore=f"{RESULTS_DIR}/bw_zscoring/custom-{{sample}}.txt",
    params:
        bigwigAverageOverBed="/zata/data/zlab/common/tools/ucsc.v385/bigWigAverageOverBed",
        bw=lambda wc: (
            pd.read_csv(
                    config["CUSTOM_data_bw_info"],
                    sep="\t",
                    header=None,
                    usecols=[0, 1],
            ).set_index(0).loc[wc.sample, 1]
        ),
    log:
        f"{LOG_DIR}/bw_zscoring/custom-{{sample}}.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "$(dirname {output.zscore})" "$(dirname {log})"

        workdir=$(mktemp -d "{TMP_DIR}/custom_bw_zscore_{wildcards.sample}.XXXXXX")
        trap 'rm -rf "$workdir"' EXIT

        echo "Processing custom sample: {wildcards.sample}" > {log}
        echo "BigWig: {params.bw}" >> {log}

        if [[ ! -f "{params.bw}" ]]; then
            echo "ERROR: Missing bigWig: {params.bw}" >> {log}
            exit 1
        fi

        {params.bigwigAverageOverBed} \
            "{params.bw}" \
            {input.regions:q} \
            "$workdir/out2" >> {log} 2>&1

        python workflow/scripts/zscore.py \
            --input "$workdir/out2" \
            --output {output.zscore:q}

        head {output.zscore} >> {log}
        """


rule ENCODE_Zscore_across_bw:
    input:
        regions=f"{RESULTS_DIR}/{PREFIX}_{GENOME}-Anchors-ATAC.regions.bed",
    output:
        zscore=f"{RESULTS_DIR}/bw_zscoring/encode-{{experiment}}_{{file}}.txt",
    params:
        dataDir="/zata/data/zlab/projects/encode/data",
        toolkit=config['TOOLKIT'],
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

        bw="{params.dataDir}/{wildcards.experiment}/{wildcards.file}.bigWig"

        if [[ ! -f "$bw" ]]; then
            echo "BigWig not found locally. Downloading {wildcards.file}" >> {log}
            python "{params.toolkit}/download-portal-file.py" \
                "{wildcards.file}" bigWig "$workdir" >> {log} 2>&1
            bw="$workdir/{wildcards.file}.bigWig"
        fi

        if [[ ! -f "$bw" ]]; then
            echo "ERROR: Could not find bigWig: $bw" >> {log}
            exit 1
        fi

        {params.bigwigAverageOverBed} \
            "$bw" \
            {input.regions:q} \
            "$workdir/out2" >> {log} 2>&1

        python workflow/scripts/zscore.py \
            --input "$workdir/out2" \
            --output {output.zscore:q}
        """


rule ENCODE_Zscore_across_DNase:
    input:
        regions=f"{RESULTS_DIR}/{PREFIX}_{GENOME}-Anchors-ATAC.regions.bed",
    output:
        zscore=f"{RESULTS_DIR}/bw_zscoring-DNAse/encode-{{experiment}}_{{file}}.txt",
    params:
        dataDir="/zata/data/zlab/projects/encode/data",
        toolkit=config['TOOLKIT'],
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

        bw="{params.dataDir}/{wildcards.experiment}/{wildcards.file}.bigWig"

        if [[ ! -f "$bw" ]]; then
            echo "BigWig not found locally. Downloading {wildcards.file}" >> {log}
            python "{params.toolkit}/download-portal-file.py" \
                "{wildcards.file}" bigWig "$workdir" >> {log} 2>&1
            bw="$workdir/{wildcards.file}.bigWig"
        fi

        if [[ ! -f "$bw" ]]; then
            echo "ERROR: Could not find bigWig: $bw" >> {log}
            exit 1
        fi

        {params.bigwigAverageOverBed} \
            "$bw" \
            {input.regions:q} \
            "$workdir/out2" >> {log} 2>&1

        python workflow/scripts/zscore.py \
            --input "$workdir/out2" \
            --output {output.zscore:q}
        """
