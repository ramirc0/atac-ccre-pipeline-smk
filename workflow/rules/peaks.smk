# rPeak clustering, filtering, accessioning, and anchors.


rule ATAC_cluster_rPeaks:
    input:
        peaks=f"{config['ouput_dir']}/{{prefix}}_{config['genome']}-ATAC-tmp.sorted.bed",
    output:
        rpeaks=f"{config['ouput_dir']}/{{prefix}}_tmp.rPeaks",
        summits=f"{config['ouput_dir']}/{{prefix}}_{config['genome']}-ATAC-tmp.summits",
    log:
        f"{config['log_dir']}/{{prefix}}_atac_cluster_rPeaks.log"
    benchmark:
        f"{RESULTS_DIR}/benchmarks/{{prefix}}_ATAC_cluster_rPeaks.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        workdir=$(mktemp -d)
        trap 'rm -rf "$workdir"' EXIT

        cp {input.peaks:q} "$workdir/remaining.bed"
        > "$workdir/rPeaks.bed"

        num=$(wc -l < "$workdir/remaining.bed")

        echo "Merging peaks..."

        while [ "$num" -gt 0 ]; do
            echo -e "\t$num"

            bedtools merge \
                -i "$workdir/remaining.bed" \
                -c 14,7 \
                -o collapse,collapse \
                > "$workdir/tmp.merge"

            python workflow/scripts/pick-best-peak.py "$workdir/tmp.merge" > "$workdir/tmp.peak-list"

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

        sort -k1,1V -k2,2n "$workdir/rPeaks.bed" > {output.rpeaks:q}

        awk -F "\t" 'BEGIN{{OFS="\t"}} {{print $1, $2+$10, $2+$10+1, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14}}' \
            {input.peaks:q} \
            > {output.summits:q}
        """


rule ATAC_filter_rPeaks:
    input:
        rpeaks=f"{config['ouput_dir']}/{{prefix}}_tmp.rPeaks",
        summits=f"{config['ouput_dir']}/{{prefix}}_{config['genome']}-ATAC-tmp.summits",
        rdhs=config["rdhs_path"],
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
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bedtools intersect -wo \
            -a {input.rpeaks:q} \
            -b {input.summits:q} \
            > {output.intersection:q}

        python workflow/scripts/filter-tf-rpeaks.py {output.intersection:q} \
            > {output.filtered:q}

        bedtools intersect -v \
            -a {output.filtered:q} \
            -b {input.rdhs:q} \
            > {output.no_rdhs:q}

        echo "No rDHS overlap:"
        wc -l {output.no_rdhs:q}

        bedtools intersect -v \
            -a {output.no_rdhs:q} \
            -b {input.blacklist:q} \
            > {output.no_rdhs_mappable:q}

        echo "Final no rDHS + no blacklist:"
        wc -l {output.no_rdhs_mappable:q}
        """


rule ATAC_make_summary_and_accession:
    input:
        no_overlap_mappable=f"{RESULTS_DIR}/{{prefix}}_{GENOME}-ATAC-tmp.no-overlap-mappable",
    output:
        summary=f"{RESULTS_DIR}/{{prefix}}_{GENOME}-ATAC-Summary.txt",
        bed=f"{RESULTS_DIR}/{{prefix}}_{GENOME}-ATAC.bed",
    log:
        f"{LOG_DIR}/{{prefix}}_atac_make_summary_bed.log"
    params:
        genome=GENOME,
    benchmark:
        f"{RESULTS_DIR}/benchmarks/{{prefix}}_ATAC_make_summary_and_accession.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/make-region-accession.py \
            {input.no_overlap_mappable:q} \
            {params.genome:q} \
            ATAC \
            > {output.summary:q}

        awk 'BEGIN{{OFS="\t"}} {{print $1, $2, $3, $NF}}' \
            {output.summary:q} \
            > {output.bed:q}
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
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        cat {input.atac_bed:q} {input.rdhs:q} \
            | sort -k1,1V -k2,2n \
            > {output.combined:q}
        """
