rule ATAC_cluster_rPeaks:
    input:
        peaks=f"{config['ouput_dir']}/{{prefix}}_{config['genome']}-ATAC-tmp.sorted.bed",
        pick_best=f"{config['TOOLKIT']}/pick-best-peak.py",
    output:
        rpeaks=f"{config['ouput_dir']}/{{prefix}}_tmp.rPeaks",
        summits=f"{config['ouput_dir']}/{{prefix}}_{config['genome']}-ATAC-tmp.summits",
    log:
        f"{config['log_dir']}/{{prefix}}_atac_cluster_rPeaks.log"
    benchmark:
        f"{RESULTS_DIR}/benchmarks/{{prefix}}_ATAC_cluster_rPeaks.tsv"
    shell:
        r"""
        set -euo pipefail

        workdir=$(mktemp -d)
        trap 'rm -rf "$workdir"' EXIT

        cp {input.peaks} "$workdir/remaining.bed"
        > "$workdir/rPeaks.bed"

        num=$(wc -l < "$workdir/remaining.bed")

        echo "Merging peaks..." > {log}

        while [ "$num" -gt 0 ]; do
            echo -e "\t$num" >> {log}

            bedtools merge \
                -i "$workdir/remaining.bed" \
                -c 14,7 \
                -o collapse,collapse \
                > "$workdir/tmp.merge"

            python {input.pick_best} "$workdir/tmp.merge" > "$workdir/tmp.peak-list"

            awk -F "\t" \
                'FNR==NR {{x[$1]; next}} ($14 in x)' \
                "$workdir/tmp.peak-list" \
                "$workdir/remaining.bed" \
                >> "$workdir/rPeaks.bed"

            sort -k1,1 -k2,2n "$workdir/rPeaks.bed" > "$workdir/rPeaks.sorted.bed"

            bedtools intersect \
                -v \
                -a "$workdir/remaining.bed" \
                -b "$workdir/rPeaks.sorted.bed" \
                > "$workdir/next.remaining.bed"

            mv "$workdir/next.remaining.bed" "$workdir/remaining.bed"

            num=$(wc -l < "$workdir/remaining.bed")
        done

        sort -k1,1V -k2,2n "$workdir/rPeaks.bed" > {output.rpeaks}

        awk -F "\t" 'BEGIN{{OFS="\t"}} {{print $1, $2+$10, $2+$10+1, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14}}' \
            {input.peaks} \
            > {output.summits}
        """