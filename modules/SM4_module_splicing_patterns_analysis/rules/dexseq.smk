import os

# Output directory for all generated files of this Tool
dexseq_output_path = os.path.join('output', config["output_directories"]["dexseq_output_dir"])

# Conditions to query
all_conditions = pep.sample_table[
        (pep.sample_table["condition"] != "None") & (pep.sample_table["condition"] != "")
    ]["condition"].unique()

# ------------------------ Define helper functions ------------------------
def get_subread_stranded_info(sample_id):
    """
    Returns, depending on settings in config file, the strandedness.
    -> Strandedness: 0 (unstranded), 1 (stranded) and 2 (reversely stranded).
    :return: Strandedness (0, 1, or 2)
    """
    strandedness = pep.sample_table[pep.sample_table["sample_name"] == sample_id]["stranded"][0]
    if strandedness == "no":
        return "0"
    elif strandedness == "yes":
        return "1"      # forward stranded
    elif strandedness == "reverse":
        return "2"      # reverse stranded
    else:
        raise ValueError(f"Strandedness value [here: {strandedness}] not specified correctly in sample table.")


def get_control_and_specific_conditioned_samples(wildcards):
    """
    Returns table of samples that are controls or have a specific condition.
    -> Returns with 2 columns: sample_name and condition

    :param wildcards: Wildcards of the snakemake rule. -> Include wildcards.condition
    :return:
    """
    condition = wildcards.condition
    samples = pep.sample_table[(pep.sample_table["condition"] == condition)
                               | (pep.sample_table["control"] == "true")][["sample_name", "condition"]]
    return samples


def get_control_and_specific_conditioned_bam_files(wildcards):
    """
    Returns a list of bam file paths of samples that have a specific condition.
    :param wildcards: Wildcards from snakemake.
    :return: List of samples that have a specific condition.
    """
    samples = get_control_and_specific_conditioned_samples(wildcards)
    sample_ids = samples["sample_name"].tolist()

    # Using helper function from Main-Snakemodule
    bam_files = main_helper_get_all_bam_file_paths(sample_ids, config["input_dir_of_bam_files"],
                                       filename_extension=config["bam_files_attributes"]["filename_extension"])
    return bam_files

def get_bai_file_paths_with_specific_condition(wildcards):
    """
    Returns a list of bai file paths of samples that have a specific condition.
    :param wildcards:
    :return:
    """
    bam_files = get_bam_file_paths_with_specific_condition(wildcards)
    bai_files = [bam_file.replace(".bam", ".bam.bai") for bam_file in bam_files]
    return bai_files


# -------------------- Rules start here ----------------------
rule dexseq_transform_gtf_chromosome_names:
    """
    This rule transforms the chromosome names in the gtf file to the ones used in aligned BAM-files.
    I.e. 1 -> chr1, 2 -> chr2, etc.
    """
    input:
        config["dexseq_settings"]["ensemble_sourced_gtf_file"]
    output:
        # Uses Main-Snakemake Script
        # Output transformed GTF -> Force conversion with bool: True
        gtf_file=main_helper_return_required_gtf_file(config["dexseq_settings"]["require_gtf_chromosome_name_conversion"],
                                          config["dexseq_settings"]["ensemble_sourced_gtf_file"],
                                          os.path.join(dexseq_output_path, "transformed_gtf_files"))
    log:
        "logs/dexseq/transform_gtf_chromosome_names.log"
    params:
        script=workflow.source_path("../../../scripts/transform_chro_names_of_gtf_file.py")
    shell:
        "python {params.script} "
        "--input_gtf_file {input} "
        "--output_gtf_file {output} "
        "--log_file {log}"


rule dexseq_prepare_annotation_for_rsubread_create_flat_gtf:
    """
    Prepares annotation file (GTF) for the usage with Rsubread:
    Similar to DexSeq-preparation script but outputs also a SubRead-featureCounts-readable GTF file.
    
    General idea:
    Collapse exons that appear multiple times, e.g. due to multiple transcripts of the same gene.
    -> Exon counting bins are generated = list of intervals corresponding to exons.
    For overlapping genes with overlapping exons: Combine genes into a single "aggregate gene", referred to 
    with the IDs joined by a plus-sign.
    """
    input:
        # Get GTF-file (after transformation if needed)
        gtf_file=main_helper_return_required_gtf_file(config["dexseq_settings"]["require_gtf_chromosome_name_conversion"],
                                          config["dexseq_settings"]["ensemble_sourced_gtf_file"],
                                          os.path.join(dexseq_output_path, "transformed_gtf_files"))
    output:
        subread_compatible_flat_gtf_file="output/dexseq/transformed_gtf_files/featureCounts_flat.gtf",
        gff_file=main_helper_return_required_gtf_file(config["dexseq_settings"]["require_gtf_chromosome_name_conversion"],
                                          config["dexseq_settings"]["ensemble_sourced_gtf_file"],
                                          os.path.join(dexseq_output_path, "transformed_gtf_files")).replace(".gtf",".gff")
    params:
        script=workflow.source_path("../scripts/dexseq_helper_fcts/dexseq_prepare_annotation2.py")
    log:
        "logs/dexseq/prepare_annotation_for_subread.log"
    threads:
        1
    conda:
        "../envs/dexseq_env.yaml"
    shell:
        "python {params.script} -f {output.subread_compatible_flat_gtf_file} "
        "{input.gtf_file} {output.gff_file} 2>{log}"



rule dexseq_count_reads_with_subread_per_sample:
    """
    For each SAM/BAM file, we count the number of reads that overlap with each of the 
    exon counting bins defined in the flattened GFF file. 
    """
    input:
        subread_compatible_flat_gtf_file=rules.dexseq_prepare_annotation_for_rsubread_create_flat_gtf.output.subread_compatible_flat_gtf_file,
        bam_file=main_helper_get_bam_file_for_sample(config["input_dir_of_bam_files"],
                    filename_extension=config["bam_files_attributes"]["filename_extension"]),
        # indexed BAM files
        bam_index_files=main_helper_get_bam_file_for_sample(config["input_dir_of_bam_files"],
                    filename_extension=config["bam_files_attributes"]["filename_extension"]).replace(".bam", ".bam.bai")
    output:
        subread_exon_counting_bin_file=os.path.join(dexseq_output_path, "feature_counts_output/{sample_id}.featureCounts")
    log:
        "logs/dexseq/{sample_id}/dexseq_count_reads_with_subread.log"
    params:
        format="GTF",
        # -r pos: BAM files are sorted by position -> report counts as positions
        order="pos" if config["bam_files_attributes"]["sorted"] == "position" else "name",
        # --paired yes: paired-end reads
        paired="-p" if config["bam_files_attributes"]["paired"] else "",
        # --stranded: Stranded? yes/no/reverse
        stranded=lambda wildcards: get_subread_stranded_info(wildcards.sample_id)
    threads:
        4
    resources:
        cpus=4,
        mem_mb=4*4096,
        time_min=120
    conda:
        "../envs/subread_env.yaml"
    shell:
        "featureCounts -f "             # -f Option to count reads overlapping features
        "-O "                           # -O Option to count reads overlapping to multiple exons
        "-J "                           # -J: Count number of reads supporting each exon-exon junction -> Creates a separate file
        #"--fracOverlap 0.2 "            # --fracOverlap FLOAT: Minimum fraction of a read that must overlap a feature to be assigned to that feature
        "-s {params.stranded} "         # -s Strandedness: 0 (unstranded), 1 (stranded) and 2 (reversely stranded).
        "{params.paired} "              # -p If specified, fragments (or templates) will be counted instead of reads. 
        "-T {threads} "                 # Specify number of threads
        "-F {params.format} "           # Specify format of the provided annotation file
        "-a {input.subread_compatible_flat_gtf_file} "          # Name of annotation file
        "-o {output.subread_exon_counting_bin_file} "           # Output file including read counts
        "{input.bam_file} "
        "2> {log}"


rule dexseq_merge_count_files:
    """
    Merge the count files of all samples into one file.
    """
    input:
        count_files=expand(os.path.join(dexseq_output_path,"feature_counts_output/{sample_id}.featureCounts"),
                            sample_id=pep.sample_table["sample_name"]),
    output:
        total_counts_file=os.path.join(dexseq_output_path, "feature_counts_output", "total_counts.tsv")
    script:
        "../scripts/dexseq/merge_feature_count_files.py"


rule dexseq_run_analysis:
    """
    Execute the dexseq analysis for a specific annotated condition.

    Idea:
    For each gene DexSeq fits 2 generalized linear models (GLM) to the read counts in the exons of the gene.
    1. Null model: ~ sample + exon
    2. Alternative model: ~ sample + exon + exon:sample
    The deviances of both models are compared using a X^2-test, providing a p-value
    """
    input:
        exon_counting_bin_file=rules.dexseq_merge_count_files.output.total_counts_file,
        flattened_gtf_file=rules.dexseq_prepare_annotation_for_rsubread_create_flat_gtf.output.subread_compatible_flat_gtf_file,
    output:
        dexseq_results_object_file=os.path.join(dexseq_output_path,
            "dexseq_analysis_results/{condition}/dexseq_results_object.rds"),
        result_summary_csv_file=os.path.join(dexseq_output_path,
            "dexseq_analysis_results/{condition}/result_summary.csv")
    log:
        "logs/dexseq/{condition}/run_analysis.log"
    params:
        load_subreadOutput_script = workflow.source_path("../scripts/dexseq_helper_fcts/load_SubreadOutput.R"),
        sample_ids=lambda wildcards: get_control_and_specific_conditioned_samples(wildcards)["sample_name"].tolist(),
        sample_conditions=lambda wildcards: get_control_and_specific_conditioned_samples(wildcards)["condition"].tolist(),
        summary_report_fdr=config["dexseq_settings"]["summary_report"]["false_discovery_rate"]
    threads:
        25
    resources:
        mem_mb=25* 4096 *2, # total: assign 4096MB per CPU = 208400
        cpus=25,            # uses 25 cpus -> Attention: Do not use too many, since then the memory will be exhausted
        time_min=12 * 60    # give max 12 hrs
    conda:
        "../envs/dexseq_env.yaml"
    script:
        "../scripts/dexseq/dexseq_data_analysis.R"


rule dexseq_extract_significant_results:
    """
    Extract the top x significant results from the dexseq analysis.
    Additionally adds the gene name to the results.
    """
    input:
        summary_csv_file=rules.dexseq_run_analysis.output.result_summary_csv_file,
    output:
        filtered_results_file=os.path.join(dexseq_output_path,
            "dexseq_analysis_results/{condition}/significant_and_named_result_summary.csv")
    log:
        "logs/dexseq/{condition}/extract_significant_results.log"
    params:
        gene_mapping_file=config["dexseq_settings"]["ensembl_gene_id_to_gene_name_mapping_file"],
        top_x_results=1000
    threads:
        1
    conda:
        "../envs/dexseq_env.yaml"
    script:
        "../scripts/dexseq/extract_significant_results.py"


rule dexseq_create_summary_report_for_workflow_report:
    input:
        rules.dexseq_extract_significant_results.output.filtered_results_file,
    output:
        report(
            directory(os.path.join(dexseq_output_path,"html_reports", "{condition}")),
            patterns=["{filename}.csv.report.html"],
            # Path relative to Snakefile
            caption="../../../report_source/module4_splicing_patterns_analysis/dexseq/caption.rst",
            category=config["output_dir_module4_splicing_patterns_analysis"].replace("output/",""),
            subcategory="DexSeq",
            labels={  # Set labels manually
                "Condition": "{condition}",
                "File:": "{filename}"
            }
        ),
        os.path.join(dexseq_output_path,"html_reports", "{condition}", "significant_and_named_result_summary.csv.report.html")
    params:
        input_files = lambda wildcards, input: input,
        data_separators=["," for i in range(1)],
        data_titles=lambda wildcards: ["DEXSeq results for condition: {}".format(wildcards.condition)],
        # Cannot use lambda function with workflow.source_path -> ToDo Buggy
        info_texts=[workflow.source_path("../../../report_source/module4_splicing_patterns_analysis/dexseq/info.html")],
        html_output_dir=lambda wildcards, output: output[0],
        html_output_file_basenames=lambda wildcards, output: [os.path.basename(file) for file in output[1:]],
    conda:
        "../../../report_source/report_env.yaml"
    script:
        "../../../scripts/create_report_html_files.R"



# =============================== Not really needed, since it's not helpful... =================================
# However, this also includes nice SVG-plots
rule dexseq_create_html_summary_reports_with_dexseq_function:
    """
    Create html summary reports for each condition.
    
    Note: This is a functionality provided by DexSeq itself!!!
    ATTENTION: Not included in workflow report yet!
    # TODO check whether to inlcude in workflow report
    """
    input:
        dexseq_results_object_file_list=expand(os.path.join(dexseq_output_path,
            "dexseq_analysis_results/{condition}/dexseq_results_object.rds"), condition=all_conditions)
    output:
        result_html_summary_report_file_list=expand(os.path.join(dexseq_output_path,
                    "dexseq_analysis_results/html_reports/{condition}/testForDEU.html"), condition=all_conditions)
    params:
        summary_report_fdr=config["dexseq_settings"]["summary_report"]["false_discovery_rate"]
    threads:
        12
    resources:
        mem_mb=12* 4096,    # total: assign 4096MB per CPU
        cpus=12,            # uses 12 cpus -> Attention: Do not use too many, since then the memory will be exhausted
        time_min=1 * 60     # give max 1 hr
    conda:
        "../envs/dexseq_env.yaml"
    script:
        "../scripts/dexseq/dexseq_create_html_summary_reports.R"
