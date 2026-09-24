# rPeak clustering, filtering, accessioning, and anchors.


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


rule ATAC_filter_rPeaks:
    input:
        rpeaks=f"{config['ouput_dir']}/{{prefix}}_tmp.rPeaks",
        summits=f"{config['ouput_dir']}/{{prefix}}_{config['genome']}-ATAC-tmp.summits",
        filter_script=f"{config['TOOLKIT']}/filter-tf-rpeaks.py",
        rdhs=config["rdhs_path"],
        mappable=config["mappable"],
        multimap="/data/projects/encode/Registry/V4/GRCh38/GRCh38-MultiMap-cCREs.bed",
        blacklist="/data/zusers/ramirezc/static/ENCFF356LFX.bed",
    output:
        intersection=f"{config['ouput_dir']}/{{prefix}}_{config['genome']}-ATAC-tmp.intersection",
        filtered=f"{config['ouput_dir']}/{{prefix}}_{config['genome']}-ATAC-tmp.out",
        no_rdhs=f"{config['ouput_dir']}/{{prefix}}_{config['genome']}-ATAC-tmp.no-overlap",
        no_rdhs_mappable=f"{config['ouput_dir']}/{{prefix}}_{config['genome']}-ATAC-tmp.no-overlap-mappable",
    log:
        f"{config['log_dir']}/{{prefix}}_atac_filter_rPeaks_no_rdhs_mappable.log"
    benchmark:
        f"{RESULTS_DIR}/benchmarks/{{prefix}}_ATAC_filter_rPeaks_no_rdhs_mappable.tsv"
    shell:
        r"""
        set -euo pipefail

        bedtools intersect -wo \
            -a {input.rpeaks} \
            -b {input.summits} \
            > {output.intersection}

        python {input.filter_script} {output.intersection} \
            > {output.filtered}

        bedtools intersect -v \
            -a {output.filtered} \
            -b {input.rdhs} \
            > {output.no_rdhs}

        echo "No rDHS overlap:" > {log}
        wc -l {output.no_rdhs} >> {log}

        bedtools intersect -v \
            -a {output.no_rdhs} \
            -b {input.blacklist} \
            > {output.no_rdhs_mappable}

        echo "Final no rDHS + mappable + no MultiMap + no blacklist:" >> {log}
        wc -l {output.no_rdhs_mappable} >> {log}
        """


rule ATAC_make_summary_and_accession:
    input:
        no_overlap_mappable=f"{RESULTS_DIR}/{{prefix}}_{GENOME}-ATAC-tmp.no-overlap-mappable",
        make_region=f"{TOOLKIT}/make-region-accession.py",
    output:
        summary=f"{RESULTS_DIR}/{{prefix}}_{GENOME}-ATAC-Summary.txt",
        bed=f"{RESULTS_DIR}/{{prefix}}_{GENOME}-ATAC.bed",
    log:
        f"{LOG_DIR}/{{prefix}}_atac_make_summary_bed.log"
    params:
        genome=GENOME,
    benchmark:
        f"{RESULTS_DIR}/benchmarks/{{prefix}}_ATAC_make_summary_and_accession.tsv"
    shell:
        r"""
        set -euo pipefail

        python {input.make_region} \
            {input.no_overlap_mappable} \
            {params.genome} \
            ATAC \
            > {output.summary} \
            2> {log}

        awk 'BEGIN{{OFS="\t"}} {{print $1, $2, $3, $NF}}' \
            {output.summary} \
            > {output.bed} \
            2>> {log}
        """


rule add_new_anchors:
    input:
        atac_bed=f"{config['ouput_dir']}/{{prefix}}_{config['genome']}-ATAC.bed",
        rdhs=config["rdhs_path"],
    output:
        combined=f"{config['ouput_dir']}/{{prefix}}_{config['genome']}-Anchors-ATAC.bed",
    log:
        f"{LOG_DIR}/{{prefix}}_add_new_anchors.log"
    benchmark:
        f"{RESULTS_DIR}/benchmarks/{{prefix}}_add_new_anchors.tsv"
    shell:
        r"""
        set -euo pipefail

        cat {input.atac_bed} {input.rdhs} \
            | sort -k1,1V -k2,2n \
            > {output.combined}
        """
