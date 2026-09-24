rule CUSTOM_Zscore_across_bw:
    input:
        anchors=f"{RESULTS_DIR}/{PREFIX}_{GENOME}-Anchors-ATAC.bed",
    output:
        zscore=f"{RESULTS_DIR}/bw_zscoring/custom-{{sample}}.txt",
    params:
        bigwigAverageOverBed="/zata/data/zlab/common/tools/ucsc.v385/bigWigAverageOverBed",
        zscoreScript=f"{TOOLKIT}/log-zscore-AF-MC.py",
        bw=lambda wc: (
            pd.read_csv(
                    config["CUSTOM_data_bw_info"],
                    sep="\t",
                    header=None,
                    usecols=[0, 1],
            ).set_index(0).loc[wc.sample, 1]
        ),
        width=0,
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

        awk -v width="{params.width}" '
        BEGIN {{ FS=OFS="\t" }}
        {{
            start = $2 - width
            end   = $3 + width
            if (start < 0) start = 0
            print $1, start, end, $4
        }}
        ' {input.anchors} | sort -u > "$workdir/little.bed"

        {params.bigwigAverageOverBed} \
            -bedOut="$workdir/out2.bed" \
            "{params.bw}" \
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

        head {output.zscore} >> {log}
        """