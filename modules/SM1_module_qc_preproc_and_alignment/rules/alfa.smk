import os

# Output directory
alfa_output_path = os.path.join('output',config["output_directories"]["alfa_output_dir"])

def get_alfa_inputs(wildcards=None):
    """
    Get the input files for the alfa algorithm.
    Either use Star or Olego alignment results as input. -> This is determined by config file
    ATTENTION: This is only needed, in case multiple BAM files are analyzed at once by ALFA.
    :param wildcards:
    :return:
    """
    input_bam_files, input_bam_bai_files = [], []
    if config["switch_variables"]["run_alignment"]["use_star"]:
        star_output_path = os.path.join('output',config["output_directories"]["star_output_dir"])
        input_bam_files += expand(os.path.join(star_output_path, "{sample_id}/{sample_id}.sorted.bam"),
                                  sample_id=pep.sample_table["sample_name"])
        input_bam_bai_files += expand(os.path.join(star_output_path, "{sample_id}/{sample_id}.sorted.bam.bai"),
                                      sample_id=pep.sample_table["sample_name"])
    elif config["switch_variables"]["run_alignment"]["use_olego"]:
        olego_output_path = os.path.join('output', config["output_directories"]["olego_output_dir"])
        input_bam_files += expand(os.path.join(olego_output_path, "{sample_id}/{sample_id}.sorted.bam"),
                                  sample_id=pep.sample_table["sample_name"])
        input_bam_bai_files += expand(os.path.join(olego_output_path, "{sample_id}/{sample_id}.sorted.bam.bai"),
                                      sample_id=pep.sample_table["sample_name"])
    else:
        raise Exception("No alignment method selected")

    return {"bam_files": input_bam_files, "bam_bai_files": input_bam_bai_files}



def main_helper_return_required_gtf_file(wildcards=None):
    """
    Return the required gtf file path, depending on the whether a chromosome name transformation is needed or not.
    ATTENTION: This takes subsequently advantage of function helper_rule_transform_chromosome_names_in_gtf in
    the main snakemake directory's scripts/helper_functions.py.
    :return:
    """
    if config["alfa_settings"]["require_gtf_chromosome_name_conversion"]:
        gtf_file = os.path.join(alfa_output_path, "transformed_gtf_files",
                                os.path.basename(config["alfa_settings"]["ensemble_sourced_gtf_file"].replace(".gtf",
                                    ".transformed_chro_names.gtf")))
    else:
        gtf_file = config["alfa_settings"]["ensemble_sourced_gtf_file"]
    return gtf_file


rule alfa_transform_gtf_chromosome_names:
    """
    This rule transforms the chromosome names in the gtf file to the ones used in aligned BAM-files.
    I.e. 1 -> chr1, 2 -> chr2, etc.
    """
    input:
        config["alfa_settings"]["ensemble_sourced_gtf_file"]
    output:
        main_helper_return_required_gtf_file()
    log:
        "logs/alfa/transform_gtf_chromosome_names.log"
    params:
        script=workflow.source_path("../../../scripts/transform_chro_names_of_gtf_file.py")
    shell:
        "python {params.script} "
        "--input_gtf_file {input} "
        "--output_gtf_file {output} "
        "--log_file {log}"


rule alfa_sort_gtf_file:
    """
    This rule sorts the gtf file.
    """
    input:
        gtf_file = main_helper_return_required_gtf_file()
    output:
        sorted_gtf_file = os.path.join(main_helper_return_required_gtf_file()).replace(".gtf",".sorted.gtf")
    log:
        "logs/alfa/sort_gtf_file.log"
    shell:
        "sort -k1,1 -k4,4n -k5,5nr {input.gtf_file} > {output.sorted_gtf_file} 2> {log}"


rule alfa_generate_alfa_index:
    """
    This rule generates the alfa index.
    """
    input:
        gtf_file=rules.alfa_sort_gtf_file.output       # sorted GTF file
    output:
        output_index_dir = directory(os.path.join(alfa_output_path, "alfa_index")),
        index_file = os.path.join(alfa_output_path, "alfa_index", 
            config["alfa_settings"]["alfa_genome_index_name"] + ".unstranded.ALFA_index")
    log:
        "logs/alfa/generate_alfa_index.log"
    params:
        alfa_genome_index_name=config["alfa_settings"]["alfa_genome_index_name"]
    threads:
        10
    conda:
        "../envs/alfa_env.yaml"
    resources:
        mem_mb = 40960,  # total: assign 4096MB per CPU
        cpus = 10,
        time_min = 60  # max 1 hr
    shell:
        "alfa -a {input.gtf_file} -g {params.alfa_genome_index_name} -o {output.output_index_dir} -p {threads} 2> {log}"


# Analysis for single BAM files -> includes category analysis that is not needed for multiple BAM files
rule alfa_assess_feature_distributions_single_bam:
    """
    Assess the feature distributions (CDS, 5'UTR, 3'UTR...) in the aligned BAM files.
    """
    input:
        generated_alfa_gtf_index = rules.alfa_generate_alfa_index.output.index_file,
        input_bam_file="output/{chosen_aligner}/{sample_id}/{sample_id}.sorted.bam",
        input_indexed_bam_file="output/{chosen_aligner}/{sample_id}/{sample_id}.sorted.bam.bai"
    output:
        report(
            directory(os.path.join(alfa_output_path,"{chosen_aligner}/{sample_id}")),
            patterns=["{filename}.tsv"],
            # Path relative to Snakefile
            caption="../../../report_source/module1_qc_preproc_alignment/alfa/caption.rst",
            category=config["output_dir_module1_qc_preprocessing_and_alignment"].replace("output/", ""),
            subcategory="ALFA: Biotypes and Categories",
            labels={
                "Used Aligner": "{chosen_aligner}",
                "Sample": "{sample_id}",
                "Output": "{filename}.tsv"
            }
        ),
        biotypes_plot = os.path.join(alfa_output_path, "{chosen_aligner}/{sample_id}/ALFA_plots.Biotypes.pdf"),
        categories_plot = os.path.join(alfa_output_path, "{chosen_aligner}/{sample_id}/ALFA_plots.Categories.pdf"),
        feature_counts_table = os.path.join(alfa_output_path, "{chosen_aligner}/{sample_id}/{sample_id}.ALFA_feature_counts.tsv")
    log:
        "logs/alfa/assess_feature_distribution_{chosen_aligner}_{sample_id}.log"
    params:
        alfa_genome_index_name= rules.alfa_generate_alfa_index.output.index_file.replace(".unstranded.ALFA_index","")
    threads:
        10
    conda:
        "../envs/alfa_env.yaml"
    resources:
        mem_mb=40960,# total: assign 4096MB per CPU
        cpus=10,		# uses 32 cpus
        time_min=120
    shell:
        "alfa -g {params.alfa_genome_index_name} "
        "--bam {input.input_bam_file} {wildcards.sample_id} "
        "-o {output[0]} "
        "--processors {threads} "
        "2> {log}"


rule alfa_assess_feature_distributions_summary:
    """
    Assess the feature distributions (CDS, 5'UTR, 3'UTR...) for all the aligned BAM files.
    """
    input:
        unpack(get_alfa_inputs),  # Returns "bam_files"
        generated_alfa_gtf_index=rules.alfa_generate_alfa_index.output.index_file,
    output:
        report(
            directory(os.path.join(alfa_output_path, "{chosen_aligner}/summary")),
            patterns=["{filename}.pdf"],
            # Path relative to Snakefile
            caption="../../../report_source/module1_qc_preproc_alignment/alfa/caption.rst",
            category=config["output_dir_module1_qc_preprocessing_and_alignment"].replace("output/", ""),
            subcategory="ALFA: Biotypes and Categories",
            labels={
                "Used Aligner": "{chosen_aligner}",
                "Sample": "All",
                "Output": "{filename}.pdf"
            }
        ),
        biotypes_plot=os.path.join(alfa_output_path, "{chosen_aligner}/summary/ALFA_plots.Biotypes.pdf"),
        categories_plot=os.path.join(alfa_output_path, "{chosen_aligner}/summary/ALFA_plots.Categories.pdf"),
        feature_counts_tables=expand(os.path.join(alfa_output_path,
            "{{chosen_aligner}}/summary/{sample_id}.ALFA_feature_counts.tsv"),
            sample_id=pep.sample_table["sample_name"])
    log:
        "logs/alfa/assess_feature_distribution_{chosen_aligner}_summary.log"
    params:
        alfa_genome_index_name=rules.alfa_generate_alfa_index.output.index_file.replace(".unstranded.ALFA_index",""),
    threads:
        10
    conda:
        "../envs/alfa_env.yaml"
    resources:
        mem_mb=40960,# total: assign 4096MB per CPU
        cpus=10,	# uses 32 cpus
        time_min=120
    script:
        "../scripts/alfa_summary_analysis.py"
