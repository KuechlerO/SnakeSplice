import os

# Output directory for all generated files of this Tool
summary_output_path = os.path.join('output', config["output_directories"]["summary_output_dir"])


def collect_result_files(wildcards):
    """
    Collects the result files from all module6 subdirectories

    :return:
    """
    current_condition = wildcards.condition
    result_files_set = {}

    # Leafcutter analysis
    if config["switch_variables"]["run_leafcutter_analysis"]:
        leafcutter_output_path = os.path.join("output",config["output_directories"]["leafcutter_output_dir"])

        result_files_set["leafcutter_results"] = os.path.join(
                leafcutter_output_path,
                config["leafcutter_settings"]["leafcutter_project_output_dir"],
                "diff_splicing_analysis", current_condition, "sample_analysis_cluster_significance_extracted_significant_p010.tsv"
            )


    # DexSeq Analysis
    if config["switch_variables"]["run_dexseq_analysis"]:
        dexseq_output_path = os.path.join('output',config["output_directories"]["dexseq_output_dir"])

        # Filtered results
        result_files_set["dexseq_results"] = os.path.join(
                dexseq_output_path,
                "dexseq_analysis_results", current_condition, "significant_and_named_result_summary.csv"
            )

    # FRASER
    if config["switch_variables"]["run_fraser_analysis"]:
        fraser_output_path = os.path.join("output",config["output_directories"]["fraser_output_dir"])

        # FRASER differential splicing analysis
        result_files_set["fraser_results"] = os.path.join(fraser_output_path, 'fraser_csv_summary_table.csv')


    # Private Junction Detection
    if config["switch_variables"]["run_private_junction_detection"]:
        # Output directory for tool
        pjd_output_path = os.path.join('output', config["output_directories"]["private_junction_detection_output_dir"])

        # Filtered results for control and condition
        result_files_set["pjd_control"] = os.path.join(pjd_output_path,
                current_condition, "filtered_junctions_0", "named_only_control_junctions_250.gene_symbol.junc"
            )
        result_files_set["pjd_condition"] = os.path.join(pjd_output_path,
                current_condition, "filtered_junctions_0", "named_only_condition_junctions_250.gene_symbol.junc"
            )

    # rMATS
    if config["switch_variables"]["run_rmats"]:
        rmats_output_path = os.path.join("output", config["output_directories"]["rmats_output_dir"])

        result_files_set["rmats_results_a3ss_jcec"] = os.path.join(rmats_output_path,
            current_condition, "filtered_A3SS.MATS.JCEC.txt")
        result_files_set["rmats_results_a5ss_jcec"] = os.path.join(rmats_output_path,
            current_condition, "filtered_A5SS.MATS.JCEC.txt")
        result_files_set["rmats_results_mxe_jcec"] = os.path.join(rmats_output_path,
            current_condition, "filtered_MXE.MATS.JCEC.txt")
        result_files_set["rmats_results_se_jcec"] = os.path.join(rmats_output_path,
            current_condition, "filtered_SE.MATS.JCEC.txt")
        result_files_set["rmats_results_ri_jcec"] = os.path.join(rmats_output_path,
            current_condition, "filtered_RI.MATS.JCEC.txt")

    return result_files_set


rule summarize_splicing_tools_results_into_one_table:
    input:
        unpack(collect_result_files)
    output:
        summary_output_table=os.path.join(summary_output_path, "{condition}", "splicing_tools_summary_table.tsv"),
    params:
        # Leafcutter
        leafcutter_gene_col_name="genes",
        leafcutter_adjusted_pval_col_name="p.adjust",
        # DexSeq
        dexseq_gene_col_name="gene_name",
        dexseq_adjusted_pval_col_name="padj",
        # FRASER
        fraser_gene_col_name="hgncSymbol",
        fraser_adjusted_pval_col_name="padjust",
        # Private Junction Detection
        pjd_gene_col_name="gene_name",
        # rMATS
        rmats_gene_col_name="geneSymbol",
        rmats_adjusted_pval_col_name="PValue",

        samples_with_condition=lambda wildcards: pep.sample_table[
            (pep.sample_table["condition"] == wildcards.condition)]["sample_name"].tolist(),
    threads:
        1
    resources:
        mem_mb=4096 * 8,
        cpus=1,
        time_min=60
    script:
        "../scripts/summary/create_splice_results_summary_file.py"


rule summarize_splicing_tools_results_create_html_reports:
    input:
        rules.summarize_splicing_tools_results_into_one_table.output.summary_output_table
    output:
        report(
            directory(os.path.join(summary_output_path, "html_reports", "{condition}")),
            patterns=["{filename}.report.html"],
            # Path relative to Snakefile
            caption="../../../report_source/module4_splicing_patterns_analysis/summary/caption.rst",
            category=config["output_dir_module4_splicing_patterns_analysis"].replace("output/",""),
            subcategory="Splicing Tools Summary",
            labels={  # Set labels manually
                "condition": "{condition}",
                "File:": "{filename}"
            }
        ),
        os.path.join(summary_output_path, "html_reports", "{condition}", "splicing_tools_summary_table.tsv.report.html"),
    params:
        input_files=lambda wildcards, input: input,
        data_separators=["\t" for i in range(2)],
        data_titles=lambda wildcards,output:
            [
                f"Summary: Results of Splicing Tools for '{wildcards.condition}' group",
            ],
        # Cannot use lambda function with workflow.source_path -> ToDo Buggy
        info_texts=[
            workflow.source_path("../../../report_source/module4_splicing_patterns_analysis/summary/info.html")
        ],
        html_output_dir=lambda wildcards, output: output[0],
        html_output_file_basenames=lambda wildcards, output: [os.path.basename(file) for file in output[1:]]
    conda:
        "../../../report_source/report_env.yaml"
    script:
        "../../../scripts/create_report_html_files.R"
