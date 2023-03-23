import os

# Output directory
deeptools_output_path = os.path.join('output', config["output_directories"]["deeptools_output_dir"])


def get_bam_file_paths():
    """
    Get the paths to the bam files inside the bam directory defined in the config file
    :param wildcards:
    :return:
    """
    # Olego used as aligner
    if config["deeptools_settings"]["input_dir_of_bam_files"] == "olego":
        olego_output_path = os.path.join('output',config["output_directories"]["olego_output_dir"])
        wildcard_path = os.path.join(olego_output_path,"{sample_id}/{sample_id}.sorted.bam")
    # STAR used as aligner
    elif config["deeptools_settings"]["input_dir_of_bam_files"] == "star":
        star_output_path = os.path.join('output',config["output_directories"]["star_output_dir"])
        wildcard_path = os.path.join(star_output_path,"{sample_id}/{sample_id}.sorted.bam")
    # Own BAM files
    else:
        wildcard_path = main_helper_get_bam_file_for_sample(pep.sample_table["sample_name"])

    return wildcard_path


rule deeptools_create_multiBamSummary:
    input:
        bam_file_paths = expand(get_bam_file_paths(), sample_id=pep.sample_table["sample_name"]),
        bam_bai_file_paths = expand(get_bam_file_paths().replace(".bam", ".bam.bai"),
            sample_id=pep.sample_table["sample_name"])
    output:
        compressed_numpy_array = os.path.join(deeptools_output_path, "multiBamSummary.npz"),
        raw_read_counts_file = os.path.join(deeptools_output_path, "multiBamSummary_raw_read_counts.tsv")
    log:
        "logs/deepTools/create_multiBamSummary.log"
    params:
        min_map_quality = config["deeptools_settings"]["min_map_quality"],
        # Consider specific region/chromosome -> Only if set in config file
        region = "--region" + config["deeptools_settings"]["region"] if config["deeptools_settings"]["region"] else ""
    threads:
        32
    resources:
        mem_mb=32 * 4096,# total: assign 4096MB per CPU
        cpus=32,# uses 1 cpus
        time_min=300  # give max 5 hrs
    conda:
        "../envs/deeptools_env.yaml"
    shell:
        "multiBamSummary bins --minMappingQuality {params.min_map_quality} {params.region} "
        "--verbose --numberOfProcessors {threads} --bamfiles {input.bam_file_paths} "
        "--outFileName {output.compressed_numpy_array} --outRawCounts {output.raw_read_counts_file} "
        "2> {log}"



rule deeptools_plotPCA:
    input:
        rules.deeptools_create_multiBamSummary.output.compressed_numpy_array
    output:
        report(
            os.path.join(deeptools_output_path, "PCA_readCounts.png"),
            # Path relative to Snakefile
            caption="../../../report_source/module1_qc_preproc_alignment/deeptools_pca/caption.rst",
            category=config["output_dir_module1_qc_preprocessing_and_alignment"].replace("output/", ""),
            subcategory="Deeptools: PCA",
        )
    log:
        "logs/deepTools/plotPCA.log"
    params:
        plot_title="PCA of read counts",
        extras="--transpose"        # Perform the PCA on the transposed matrix -> Plot only makes sense that way...
    threads:
        1
    resources:
        mem_mb=1 * 4096,        # total: assign 4096MB per CPU
        cpus=1,                 # uses 1 cpus
        time_min=60             # give max 1 hrs
    conda:
        "../envs/deeptools_env.yaml"
    shell:
        "plotPCA -in {input} -o {output[0]} --plotTitle \"{params.plot_title}\" {params.extras} 2> {log}"
