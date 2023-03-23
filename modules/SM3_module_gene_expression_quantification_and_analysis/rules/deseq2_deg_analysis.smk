import os

# ============= Due to not supported wildcards in reports kallisto and salmon rules are defined separately =============

# ------------------ 1. Run gene level analysis ------------------
rule deseq2_deg_analysis_salmon_run_gene_level_exp_analysis:
    input:
        deseq_dataset_r_obj=os.path.join('output', "salmon", "{condition}", "salmon_deseq_dataset_object.rds"),
    output:
        report(
            directory(os.path.join('output', "salmon", "{condition}", "deseq2_gene_level_analysis")),
            patterns=["{filename}.jpg"],
            # Path relative to Snakefile
            caption="../../../report_source/module3_gene_expression_quantification_and_analysis/salmon/caption.rst",
            category=config["output_dir_module3_gene_expression_quantification_and_analysis"].replace("output/",""),
            subcategory="Salmon: Differentially Expressed Genes (DEG)",
            labels={  # Set labels manually
                "Tool": "Salmon-DESeq2",
                "Condition": "{condition}",
                "Modus": "DEG",
                "Type": "Plot",
                "File:": "{filename}",
            }
        ),
        os.path.join("output", "salmon", "{condition}", "deseq2_gene_level_analysis", "statistical_analysis_summary.txt"),
        os.path.join("output", "salmon", "{condition}", "deseq2_gene_level_analysis", "ma_plot_not_shrinked.jpg"),
        os.path.join("output", "salmon", "{condition}", "deseq2_gene_level_analysis", "ma_plot_ashr_shrinked.jpg"),
        os.path.join("output", "salmon", "{condition}", "deseq2_gene_level_analysis", "pvalue_histogram.jpg"),
        os.path.join("output", "salmon", "{condition}", "deseq2_gene_level_analysis", "pvalue_adjusted_histogram.jpg"),
        os.path.join("output", "salmon", "{condition}", "deseq2_gene_level_analysis", "gene_cluster_heatmap.jpg"),
        os.path.join("output", "salmon", "{condition}", "deseq2_gene_level_analysis", "deseq2_results.csv"),
    log:
        "logs/salmon_run_gene_level_exp_analysis_with_deseq2/{condition}.log"
    params:
        output_file_paths = lambda wildcards,output: output[1:],
        used_algorithm = "salmon",
    conda:
        "../envs/deseq2_gene_level_env.yaml"
    threads:
        12
    resources:
        # TODO push it up...
        mem_mb=12 * 4096 *4,# total: assign 4096MB per CPU
        cpus=12,	# uses 12 cpus
        time_min=360		    # use 6 hrs, to make sure
    script:
        "../scripts/deseq2/gene_expression_analysis_with_deseq2.R"


rule deseq2_deg_analysis_kallisto_run_gene_level_exp_analysis:
    input:
        deseq_dataset_r_obj=os.path.join('output', "kallisto", "{condition}", "kallisto_deseq_dataset_object.rds"),
    output:
        report(
            directory(os.path.join('output', "kallisto", "{condition}", "deseq2_gene_level_analysis")),
            patterns=["{filename}.jpg"],
            # Path relative to Snakefile
            caption="../../../report_source/module3_gene_expression_quantification_and_analysis/kallisto/caption.rst",
            category=config["output_dir_module3_gene_expression_quantification_and_analysis"].replace("output/",""),
            subcategory="Kallisto: Differentially Expressed Genes (DEG)",
            labels={  # Set labels manually
                "Tool": "Kallisto-DESeq2",
                "Condition": "{condition}",
                "Modus": "DEG",
                "Type": "Plot",
                "File:": "{filename}",
            }
        ),
        os.path.join("output", "kallisto", "{condition}", "deseq2_gene_level_analysis", "statistical_analysis_summary.txt"),
        os.path.join("output", "kallisto", "{condition}", "deseq2_gene_level_analysis", "ma_plot_not_shrinked.jpg"),
        os.path.join("output", "kallisto", "{condition}", "deseq2_gene_level_analysis", "ma_plot_ashr_shrinked.jpg"),
        os.path.join("output", "kallisto", "{condition}", "deseq2_gene_level_analysis", "pvalue_histogram.jpg"),
        os.path.join("output", "kallisto", "{condition}", "deseq2_gene_level_analysis", "pvalue_adjusted_histogram.jpg"),
        os.path.join("output", "kallisto", "{condition}", "deseq2_gene_level_analysis", "gene_cluster_heatmap.jpg"),
        os.path.join("output", "kallisto", "{condition}", "deseq2_gene_level_analysis", "deseq2_results.csv"),
    log:
        "logs/kallisto_run_gene_level_exp_analysis_with_deseq2/{condition}.log"
    params:
        output_file_paths = lambda wildcards,output: output[1:],
        used_algorithm = "kallisto",
    conda:
        "../envs/deseq2_gene_level_env.yaml"
    threads:
        12
    resources:
        # TODO push it up...
        mem_mb=12 * 4096 *4,# total: assign 4096MB per CPU
        cpus=12,	# uses 12 cpus
        time_min=360		    # use 6 hrs, to make sure
    script:
        "../scripts/deseq2/gene_expression_analysis_with_deseq2.R"


# ----------------- 2. Extract significant results ----------------- #
rule deseq2_deg_analysis_extract_significant_results:
    input:
        deseq2_results=os.path.join("output", "{quant_tool}", "{condition}", "deseq2_gene_level_analysis", "deseq2_results.csv"),
    output:
        os.path.join("output", "{quant_tool}", "{condition}", "deseq2_gene_level_analysis", "deseq2_results_significant.csv"),
    run:
        import pandas as pd

        df = pd.read_csv(input.deseq2_results, sep=",")
        # filter for significant results (adjusted p-value < 0.10)
        df = df[df["padj"] < 0.10]
        # rename first column
        df = df.rename(columns={df.columns[0]: "subject"})

        # write to file
        df = df.sort_values(by=["padj"])
        df.to_csv(output[0], sep=",", index=False)


# ----------------- 3. Create HTML reports ----------------- #
rule deseq2_deg_analysis_salmon_create_html_reports:
    input:
        os.path.join("output", "salmon", "{condition}", "deseq2_gene_level_analysis", "deseq2_results_significant.csv"),
    output:
        report(
            directory(os.path.join("output", "salmon", "{condition}", "html_reports")),
            patterns=["{filename}.csv.report.html"],
            # Path relative to Snakefile
            caption="../../../report_source/module3_gene_expression_quantification_and_analysis/salmon/caption.rst",
            category=config["output_dir_module3_gene_expression_quantification_and_analysis"].replace("output/",""),
            subcategory="Salmon: Differentially Expressed Genes (DEG)",
            labels={  # Set labels manually
                "Tool": "Salmon-DESeq2",
                "Condition": "{condition}",
                "Modus": "DEG",
                "Type": "Report",
                "File:": "{filename}",
            }
        ),
        os.path.join("output", "salmon", "{condition}", "html_reports", "deseq2_results_significant.csv.report.html"),
    params:
        input_files=lambda wildcards, input: input,
        data_separators=["," for i in range(1)],
        data_titles=lambda wildcards: ["Salmon: Differentially Expressed Genes (DESeq2) for " + wildcards.condition],
        # Combining snakemake.source_path and lambda does not work -> ToDo Report bug
        info_texts=[
            workflow.source_path(
                "../../../report_source/module3_gene_expression_quantification_and_analysis/salmon/info.html")
        ],
        html_output_dir=lambda wildcards, output: output[0],
        html_output_file_basenames=lambda wildcards, output: [os.path.basename(out_file) for out_file in output[1:]]
    conda:
        "../../../report_source/report_env.yaml"
    script:
        "../../../scripts/create_report_html_files.R"


rule deseq2_deg_analysis_kallisto_create_html_reports:
    input:
        os.path.join("output", "kallisto", "{condition}", "deseq2_gene_level_analysis", "deseq2_results_significant.csv"),
    output:
        report(
            directory(os.path.join("output", "kallisto", "{condition}", "html_reports")),
            patterns=["{filename}.csv.report.html"],
            # Path relative to Snakefile
            caption="../../../report_source/module3_gene_expression_quantification_and_analysis/kallisto/caption.rst",
            category=config["output_dir_module3_gene_expression_quantification_and_analysis"].replace("output/",""),
            subcategory="Kallisto: Differentially Expressed Genes (DEG)",
            labels={  # Set labels manually
                "Tool": "Kallisto-DESeq2",
                "Condition": "{condition}",
                "Modus": "DEG",
                "Type": "Report",
                "File:": "{filename}",
            }
        ),
        os.path.join("output", "kallisto", "{condition}", "html_reports", "deseq2_results_significant.csv.report.html"),
    params:
        input_files=lambda wildcards, input: input,
        data_separators=["," for i in range(1)],
        data_titles=lambda wildcards: ["Kallisto: Differentially Expressed Genes (DESeq2) for " + wildcards.condition],
        # Combining snakemake.source_path and lambda does not work -> ToDo Report bug
        info_texts=[
            workflow.source_path(
                "../../../report_source/module3_gene_expression_quantification_and_analysis/kallisto/info.html")
        ],
        html_output_dir=lambda wildcards, output: output[0],
        html_output_file_basenames=lambda wildcards, output: [os.path.basename(out_file) for out_file in output[1:]]
    conda:
        "../../../report_source/report_env.yaml"
    script:
        "../../../scripts/create_report_html_files.R"