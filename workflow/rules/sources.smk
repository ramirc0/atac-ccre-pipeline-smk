# Peak inputs: ENCODE peaks, custom peaks, and their merge.


rule ENCODE_prepare_ATAC_peaks:
    input:
        [f"{ENCODE_PEAK_DIR}/{exp}/{peak}.bed.gz" for exp, peak, _ in ENCODE_ATAC],
    output:
        tmp_bed=f"{RESULTS_DIR}/ENCODE_{GENOME}-ATAC-tmp.bed",
        sorted_bed=f"{RESULTS_DIR}/ENCODE_{GENOME}-ATAC-tmp.sorted.bed",
    params:
        samples=[exp for exp, _, _ in ENCODE_ATAC],
        biosamples=[biosample for _, _, biosample in ENCODE_ATAC],
    log:
        f"{LOG_DIR}/ENCODE_prepare_peaks.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/ENCODE_prepare_ATAC_peaks.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/prepare_peaks.py \
            --peaks {input:q} \
            --samples {params.samples:q} \
            --biosamples {params.biosamples:q} \
            --output {output.tmp_bed:q}

        sort -k1,1V -k2,2n {output.tmp_bed:q} > {output.sorted_bed:q}
        """


rule CUSTOM_prepare_peaks:
    input:
        list(CUSTOM_PEAKS.values()),
    output:
        tmp_bed=f"{config['ouput_dir']}/CUSTOM_{config['genome']}-ATAC-tmp.bed",
        sorted_bed=f"{config['ouput_dir']}/CUSTOM_{config['genome']}-ATAC-tmp.sorted.bed"
    params:
        samples=list(CUSTOM_PEAKS),
    log:
        f"{config['log_dir']}/CUSTOM_prepare_peaks.log"
    benchmark:
        f"{RESULTS_DIR}/benchmarks/CUSTOM_prepare_peaks.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/prepare_peaks.py \
            --peaks {input:q} \
            --samples {params.samples:q} \
            --output {output.tmp_bed:q}

        sort -k1,1V -k2,2n {output.tmp_bed:q} > {output.sorted_bed:q}
        """


rule merge_preped_peaks:
    input:
        ENCODE_peaks=f"{config['ouput_dir']}/ENCODE_{config['genome']}-ATAC-tmp.sorted.bed",
        CUSTOM_peaks=f"{config['ouput_dir']}/CUSTOM_{config['genome']}-ATAC-tmp.sorted.bed",
    output:
        MERGED_peaks=f"{config['ouput_dir']}/MERGED_{config['genome']}-ATAC-tmp.sorted.bed",
    log:
        f"{config['log_dir']}/MERGED_preped_peaks.log",
    benchmark:
        f"{RESULTS_DIR}/benchmarks/merge_preped_peaks.tsv"
    conda:
        CONDA_ENV
    shell:
        r"""
        set -euo pipefail

        cat {input.ENCODE_peaks} {input.CUSTOM_peaks} \
            | sort -k1,1 -k2,2n \
            > {output.MERGED_peaks} 2> {log}
        """
