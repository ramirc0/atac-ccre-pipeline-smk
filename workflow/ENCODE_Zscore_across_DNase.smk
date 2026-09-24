rule ENCODE_Zscore_across_DNase:
    input:
        anchors=f"{RESULTS_DIR}/{PREFIX}_{GENOME}-Anchors-ATAC.bed",
    output:
        zscore=f"{RESULTS_DIR}/bw_zscoring-DNAse/encode-{{experiment}}_{{file}}.txt",
    params:
        dataDir="/zata/data/zlab/projects/encode/data",
        width=0,
        toolkit=config['TOOLKIT'],
        zscoreScript=f"{TOOLKIT}/log-zscore-AF-MC.py",
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

        awk -v width="{params.width}" '
        BEGIN {{ FS=OFS="\t" }}
        {{
            start = $2 - width
            end   = $3 + width
            if (start < 0) start = 0
            print $1, start, end, $4
        }}
        ' {input.anchors} | sort -u > "$workdir/little.bed"

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
            -bedOut="$workdir/out2.bed" \
            "$bw" \
            "$workdir/little.bed" \
            "$workdir/out2" >> {log} 2>&1

        python "{params.zscoreScript}" "$workdir/out2" > "$workdir/anchor_signal.tsv"

        sort -k2,2rg "$workdir/anchor_signal.tsv" \
            | awk '
                BEGIN {{ FS=OFS="\t"; rank=0; before=""; running=1 }}
                {{
                    if ($2 != before) rank = running
                    print $1, $2, $3, rank
                    before = $2
                    running++
                }}
            ' \
            | sort -k1,1 \
            > {output.zscore}
        """