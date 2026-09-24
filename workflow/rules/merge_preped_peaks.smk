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
    shell:
        r"""
        set -euo pipefail

        cat {input.ENCODE_peaks} {input.CUSTOM_peaks} \
            | sort -k1,1 -k2,2n \
            > {output.MERGED_peaks} 2> {log}
        """