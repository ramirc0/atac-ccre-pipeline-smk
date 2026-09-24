# Peaks to rPeaks to ATAC anchors: standardize, cluster, filter, accession.


rule prepare_peaks:
    input:
        [NARROWPEAK_OF[s] for s in SAMPLES],
    output:
        peaks=f"{OUTDIR}/peaks/peaks.bed",
        sorted=f"{OUTDIR}/peaks/peaks.sorted.bed",
    params:
        samples=SAMPLES,
        biosamples=[BIOSAMPLE_OF[s] for s in SAMPLES],
    log:
        f"{LOGDIR}/prepare_peaks.txt",
    benchmark:
        f"{BENCHDIR}/prepare_peaks.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/prepare_peaks.py \
            --peaks {input:q} \
            --samples {params.samples:q} \
            --biosamples {params.biosamples:q} \
            --output {output.peaks:q}

        sort -k1,1V -k2,2n {output.peaks:q} > {output.sorted:q}
        """


# Iteratively keep the strongest peak of each overlapping cluster.
rule cluster_rpeaks:
    input:
        f"{OUTDIR}/peaks/peaks.sorted.bed",
    output:
        rpeaks=f"{OUTDIR}/rpeaks/rpeaks.bed",
        summits=f"{OUTDIR}/peaks/summits.bed",
    log:
        f"{LOGDIR}/cluster_rpeaks.txt",
    benchmark:
        f"{BENCHDIR}/cluster_rpeaks.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        workdir=$(mktemp -d)
        trap 'rm -rf "$workdir"' EXIT

        cp {input:q} "$workdir/remaining.bed"
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

            python workflow/scripts/pick_best_peak.py "$workdir/tmp.merge" > "$workdir/tmp.peak-list"

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
            {input:q} \
            > {output.summits:q}
        """


# Keep rPeaks supported by >= 5 experiments, off rDHSs and the blacklist.
rule filter_rpeaks:
    input:
        rpeaks=f"{OUTDIR}/rpeaks/rpeaks.bed",
        summits=f"{OUTDIR}/peaks/summits.bed",
        rdhs=REFS["rdhs"],
        blacklist=REFS["blacklist"],
    output:
        intersection=f"{OUTDIR}/rpeaks/intersection.tsv",
        filtered=f"{OUTDIR}/rpeaks/filtered.tsv",
        no_rdhs=f"{OUTDIR}/rpeaks/no-rdhs.tsv",
        kept=f"{OUTDIR}/rpeaks/no-rdhs-blacklist.tsv",
    log:
        f"{LOGDIR}/filter_rpeaks.txt",
    benchmark:
        f"{BENCHDIR}/filter_rpeaks.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bedtools intersect -wo \
            -a {input.rpeaks:q} \
            -b {input.summits:q} \
            > {output.intersection:q}

        python workflow/scripts/filter_rpeaks.py {output.intersection:q} \
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
            > {output.kept:q}

        echo "Final no rDHS + no blacklist:"
        wc -l {output.kept:q}
        """


rule accession_rpeaks:
    input:
        f"{OUTDIR}/rpeaks/no-rdhs-blacklist.tsv",
    output:
        summary=f"{OUTDIR}/ATAC-Summary.txt",
        bed=f"{OUTDIR}/ATAC.bed",
    params:
        genome=GENOME,
    log:
        f"{LOGDIR}/accession_rpeaks.txt",
    benchmark:
        f"{BENCHDIR}/accession_rpeaks.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/accession_regions.py {input:q} \
            --genome {params.genome:q} \
            --mode ATAC \
            > {output.summary:q}

        awk 'BEGIN{{OFS="\t"}} {{print $1, $2, $3, $NF}}' \
            {output.summary:q} \
            > {output.bed:q}
        """


rule anchors:
    input:
        atac=f"{OUTDIR}/ATAC.bed",
        rdhs=REFS["rdhs"],
    output:
        f"{OUTDIR}/Anchors-ATAC.bed",
    log:
        f"{LOGDIR}/anchors.txt",
    benchmark:
        f"{BENCHDIR}/anchors.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        cat {input.atac:q} {input.rdhs:q} \
            | sort -k1,1V -k2,2n \
            > {output:q}
        """
