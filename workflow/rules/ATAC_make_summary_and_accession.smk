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