rule ENCODE_make_filtered_ATAC_list:
    output:
        filtered_list = f"{config['ouput_dir']}/{config['genome']}-ATAC-List.Filtered.txt"

    params:
        genome=config["genome"],
        work_dir=config["work_dir"],
        filter_FRiP=config["filter_FRIP"],
        TOOLKIT=config["TOOLKIT"],
        filter_min_reads=config["filter_min_reads"],
    log:
        f"{LOG_DIR}/ENCODE_make_filtered_list.log"
    benchmark:
        f"{RESULTS_DIR}/benchmarks/ENCODE_make_filtered_ATAC_list.tsv"

    shell:
        r"""
        python -u {params.TOOLKIT}/pull-atac-experiments.py {params.genome} > {params.work_dir}/test.txt 2> {log}

        awk -F "\t" 'NR>1 && $7 >= {params.filter_FRiP} && $8 >= {params.filter_min_reads}  && $2 != "NA" {print $0}' {params.work_dir}/test.txt \
            > {output.filtered_list}

        """
