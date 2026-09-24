# Fragment inputs (e.g. snATAC pseudobulks) to the narrowPeak + fold-enrichment
# bigWig pair that bulk rows provide directly.


rule fragment_call_peaks:
    input:
        lambda w: FRAGMENTS_OF[w.sample],
    output:
        narrowpeak=f"{OUTDIR}/fragments/{{sample}}/{{sample}}_peaks.narrowPeak",
        pileup=f"{OUTDIR}/fragments/{{sample}}/{{sample}}_treat_pileup.bdg",
        lambda_bdg=f"{OUTDIR}/fragments/{{sample}}/{{sample}}_control_lambda.bdg",
    params:
        outdir=lambda w: f"{OUTDIR}/fragments/{w.sample}",
        flags=lambda w: macs3_flags(w.sample),
    log:
        f"{LOGDIR}/fragment_call_peaks/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/fragment_call_peaks/{{sample}}.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        macs3 callpeak \
            -t {input:q} \
            -n {wildcards.sample:q} \
            --outdir {params.outdir:q} \
            -B \
            --SPMR \
            {params.flags:q}
        """


rule fragment_bigwig:
    input:
        pileup=f"{OUTDIR}/fragments/{{sample}}/{{sample}}_treat_pileup.bdg",
        lambda_bdg=f"{OUTDIR}/fragments/{{sample}}/{{sample}}_control_lambda.bdg",
        chrom_sizes=lambda w: REFS["chrom_sizes"],
    output:
        f"{OUTDIR}/fragments/{{sample}}/{{sample}}_FE.bw",
    log:
        f"{LOGDIR}/fragment_bigwig/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/fragment_bigwig/{{sample}}.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        workdir=$(mktemp -d)
        trap 'rm -rf "$workdir"' EXIT

        macs3 bdgcmp \
            -t {input.pileup:q} \
            -c {input.lambda_bdg:q} \
            -m FE \
            -o "$workdir/FE.unsorted.bdg"

        # Clip to chromosome bounds; drop unknown contigs and empty intervals.
        awk 'BEGIN {{ OFS="\t" }}
            NR == FNR {{ size[$1] = $2; next }}
            ($1 in size) && $2 >= 0 && $2 < size[$1] && $3 > $2 {{
                end = ($3 > size[$1] ? size[$1] : $3)
                if (end > $2) print $1, $2, end, $4
            }}' {input.chrom_sizes:q} "$workdir/FE.unsorted.bdg" \
            | LC_ALL=C sort -k1,1 -k2,2n \
            > "$workdir/FE.bdg"

        bedGraphToBigWig "$workdir/FE.bdg" {input.chrom_sizes:q} {output:q}
        """
