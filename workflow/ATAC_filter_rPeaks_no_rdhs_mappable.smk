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