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