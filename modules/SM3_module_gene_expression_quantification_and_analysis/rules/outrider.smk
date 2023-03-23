# Identification of aberrant gene expression in RNA-seq data. Read count expectations are modeled
# by an autoencoder to control for confounders in the data. Given these expectations, the RNA-seq read
# counts are assumed to follow a negative binomial distribution with a gene-specific dispersion.
# Outliers are then identified as read counts that significantly deviate from this distribution. Furthermore,
# OUTRIDER provides useful plotting functions to analyze and visualize the results.

# Colab Tutorial: https://colab.research.google.com/drive/1OKT32eNIq7Cz839jjqz-GJlvoToPYbib

import os

# Output directory
outrider_output_path = os.path.join('output', config["output_directories"]["outrider_output_dir"])


rule outrider_transform_gtf_chromosome_names:
    """
    This rule transforms the chromosome names in the gtf file to the ones used in aligned BAM-files.
    I.e. 1 -> chr1, 2 -> chr2, etc.
    """
    input:
        config["outrider_settings"]["ensemble_sourced_gtf_file"]
    output:
        # Output transformed GTF -> Forced by bool: True
        gtf_file=main_helper_return_required_gtf_file(True,config["outrider_settings"]["ensemble_sourced_gtf_file"],
            os.path.join(outrider_output_path, "transformed_gtf_files"))
    log:
        "logs/outrider/transform_gtf_chromosome_names.log"
    params:
        script=workflow.source_path("../../../scripts/transform_chro_names_of_gtf_file.py")
    shell:
        "python {params.script} "
        "--input_gtf_file {input} "
        "--output_gtf_file {output} "
        "--log_file {log}"


def get_strandedness_type(sample_id):
    """
    Returns the library type for the given sample id.
    """
    strandedness = pep.sample_table[pep.sample_table["sample_name"] == sample_id]["stranded"][0]
    if strandedness in ["no", "yes", "reverse"]:
        return strandedness
    else:
        raise ValueError("Unknown strandedness: {}".format(strandedness))


rule outrider_create_count_table_per_sample_with_htseq:
    """
    This rule creates a count table per sample using HTSeq.
    Per sample is necessary even though HTseq can handle multiple BAM files at once, because
    the input files might be of different library types.
    """
    input:
        bam_file=main_helper_get_bam_file_for_sample(config["input_dir_of_bam_files"],
            filename_extension=config["bam_files_attributes"]["filename_extension"]),
        gtf_file=main_helper_return_required_gtf_file(config["outrider_settings"]["require_gtf_chromosome_name_conversion"],
            config["outrider_settings"]["ensemble_sourced_gtf_file"], "output/outrider/transformed_gtf_files")
    output:
        count_file = os.path.join(outrider_output_path, "count_files", "{sample_id}_feature_counts_table.tsv")
    log:
        "logs/outrider/create_count_table_{sample_id}.log"
    params:
        format="bam",   # Alignment files are BAM-files
        # Sorting order of BAM-files
        order="pos" if config["bam_files_attributes"]["sorted"] == "position"
                    else config["bam_files_attributes"]["sorted"],
        stranded=lambda wildcards: get_strandedness_type(wildcards.sample_id),
        add_attribute="gene_name"       # Also extract gene name
    threads:
        25
    resources:
        mem_mb = 25 * 4096,          # total: assign 4096MB per CPU
        cpus=25,                     # HT-Seq can use on thread per BAM-file
        time_min=6* 60                # give max 6 hrs
    conda:
        "../envs/htseq_env.yaml"
    shell:
        "htseq-count "
        "-n {threads} "
        "--format {params.format} "
        "--order {params.order} "
        "--stranded {params.stranded} "
        "--additional-attr {params.add_attribute} "
        "{input.bam_file} {input.gtf_file} >{output.count_file} 2>{log}"


rule outrider_merge_count_tables:
    """
    This rule merges the count tables of all samples into one count table.
    """
    input:
        sample_count_files=expand(os.path.join(outrider_output_path, "count_files", "{sample_id}_feature_counts_table.tsv"),
            sample_id=pep.sample_table["sample_name"])
    output:
        total_counts_table=os.path.join(outrider_output_path, "count_files", "total_count_table.tsv")
    params:
        sample_names=pep.sample_table["sample_name"],
        # additional_count_file=config["outrider_settings"]["additional_count_file"]
    log:
        "logs/outrider/merge_count_tables.log"
    conda:
        "../envs/outrider_env.yaml"
    script:
        "../scripts/outrider/merge_htseq_count_files.R"


rule outrider_modify_count_file:
    """
    Drops 1st column (to only keep gene names), sets header with sample-names as column names, and
    removes last 5 lines, which hold special cases: 
    __no_feature: Reads that could not be assigned to any feature
    __ambiguous: Reads that could be assigned to multiple features
    __too_low_aQual: Reads that were discarded because of low alignment quality
    __not_aligned: Reads in the SAM/BAM file without alignment
    __alignment_not_unique: Reads with more than one reported alignment
    """
    input:
        rules.outrider_merge_count_tables.output
    output:
        os.path.join(outrider_output_path, "count_files", "feature_counts_table_modified.tsv")
    log:
        "logs/outrider/outrider_modify_count_file"
    shell:
        # 0. Drop HTSeq-counts summary lines (e.g. __ambiguous)
        # 1. awk: Set second column to be combination of first & second: To make it unique!
        # 2. cut: Drop first column (which is Ensembl IDs, keep instead only Ensemble-ID-gene-name column)
        "grep -v \"^__\" {input} | "
        "awk 'BEGIN {{FS=OFS=\"\t\"}} {{$2=($1\"-\"$2)}} 1' - | "
        "cut -f 2- - >{output} 2> {log}"


rule outrider_create_expression_outlier_analysis_object:
    """
    Execute expression outlier identification via Outrider
    """
    input:
        counts_file=rules.outrider_modify_count_file.output,
    output:
        # Size Factors
        outrider_obj_with_estimated_size_factors_txt=os.path.join(outrider_output_path, "outrider_obj_with_estimated_size_factors.txt"),
        outrider_obj_with_estimated_size_factors_rds=os.path.join(outrider_output_path, "outrider_obj_with_estimated_size_factors.RDS"),

        # Final Outrider object
        outrider_object_file=os.path.join(outrider_output_path,"outrider_analysis_object.RData"),
    log:
        "logs/outrider/expression_outlier_analysis.log"
    conda:
        "../envs/outrider_env.yaml"
    threads:
        1
    resources:
        mem_mb=4*4096,      # total: assign 4096MB per CPU
        cpus=1,             # uses 1 cpus
        time_min=5*60        # give max 5 hrs
    script:
        "../scripts/outrider/outrider_create_analysis_object.R"


rule outrider_explore_analysis_results:
    """
    Explore Outrider analysis results
    """
    input:
        outrider_object_file=rules.outrider_create_expression_outlier_analysis_object.output.outrider_object_file,
    output:
        # Plots -> Genes of interest need to be defined in configuration file
        # 1. Per sample: Vulcano Plots of p-values (html & png)
        # 2. Per Gene: Expression ranking (html & jpg)
        # 3. Per Gene: Quantile-Quantile-Plots (Q-Q plots) (html & jpg)
        # 4. Per Gene: Observed versus expected Expression (html & jpg)
        report(
            directory(os.path.join(outrider_output_path, "result_exploration")),
            patterns=["{filename}.jpg"],
            # Path relative to Snakefile
            caption="../../../report_source/module3_gene_expression_quantification_and_analysis/outrider/caption.rst",
            category=config["output_dir_module3_gene_expression_quantification_and_analysis"].replace("output/",""),
            subcategory="OUTRIDER",
            labels={  # Set labels manually
                "Tool": "OUTRIDER",
                "Type": "Plot",
                "File:": "{filename}"
            }
        ),
        # Table data
        significant_results_p005_file=os.path.join(outrider_output_path,"result_exploration/significant_results_p005.tsv"),
        significant_results_p010_file=os.path.join(outrider_output_path,"result_exploration/significant_results_p010.tsv"),
        nr_aberrant_genes_per_sample=os.path.join(outrider_output_path,"result_exploration/nr_aberrant_genes_per_sample.txt"),
        nr_aberrant_samples_per_gene=os.path.join(outrider_output_path,"result_exploration/nr_aberrant_samples_per_gene.txt"),
    log:
        "logs/outrider/explore_analysis_results.log"
    params:
        sample_ids = pep.sample_table["sample_name"],
        genes_of_interest=config["outrider_settings"]["genes_of_interest"]
        # annotation_file = pep.sample_table
    conda:
        "../envs/outrider_env.yaml"
    threads:
        1
    resources:
        mem_mb=4 * 4096,# total: assign 4096MB per CPU
        cpus=1,# uses 1 cpus
        time_min=60  # give max 5 hrs
    script:
        "../scripts/outrider/outrider_explore_results.R"


rule outrider_create_html_reports:
    input:
        # rules.outrider_explore_analysis_results.output.significant_results_p005_file,
        rules.outrider_explore_analysis_results.output.significant_results_p010_file,
    output:
        report(
            directory(os.path.join(outrider_output_path, "html_reports")),
            patterns=["{filename}.tsv.report.html"],
            # Path relative to Snakefile
            caption="../../../report_source/module3_gene_expression_quantification_and_analysis/outrider/caption.rst",
            category=config["output_dir_module3_gene_expression_quantification_and_analysis"].replace("output/",""),
            subcategory="OUTRIDER",
            labels={  # Set labels manually
                "Tool": "OUTRIDER",
                "Type": "Report",
                "File:": "{filename}"
            }
        ),
        # os.path.join(outrider_output_path, "html_reports", "significant_results_p005.tsv.report.html"),
        os.path.join(outrider_output_path, "html_reports", "significant_results_p010.tsv.report.html"),
    params:
        input_files = lambda wildcards, input: input,
        data_separators=["\t" for i in range(4)],
        data_titles = ["Outrider: Significant results (p-adjusted<0.10)"],
        info_texts=[workflow.source_path("../../../report_source/module3_gene_expression_quantification_and_analysis/outrider/info.html")],
        html_output_dir=lambda wildcards, output: output[0],
        html_output_file_basenames = lambda wildcards, output: [os.path.basename(file) for file in output[1:]],
    conda:
        "../../../report_source/report_env.yaml"
    script:
        "../../../scripts/create_report_html_files.R"

