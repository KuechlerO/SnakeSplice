import os


# ============= Due to not supported wildcards in reports kallisto and salmon rules are defined separately =============

# ============================= Salmon: GSEA analysis =============================
rule clusterProfiler_salmon_deseq2_gene_set_enrichment_analysis:
    input:
        deseq_dataset_obj=os.path.join("output", "salmon", "{condition}", "salmon_deseq_dataset_object.rds")
    output:
        report(
            directory(os.path.join("output", "salmon", "{condition}", "gsea_analysis", "plots")),
            patterns=["gsea_{used_db}-{filename}.jpg"],
            # Path relative to Snakefile
            caption="../../../report_source/module3_gene_expression_quantification_and_analysis/gsea/caption.rst",
            category=config["output_dir_module3_gene_expression_quantification_and_analysis"].replace("output/", ""),
            subcategory="Salmon: Gene Set Enrichment Analysis (GSEA)",
            labels={  # Set labels manually
                "Tool": "Salmon-ClusterProfiler",
                "Database": "{used_db}",
                "Condition": "{condition}",
                "Type": "Plot",
                "File:": "gsea_{used_db}-{filename}.jpg",
            }
        ),
        gsea_go_obj_file_path = os.path.join("output", "salmon", "{condition}", "gsea_analysis", "gsea_go_obj.rds"),
        gsea_kegg_obj_file_path = os.path.join("output", "salmon", "{condition}", "gsea_analysis", "gsea_kegg_obj.rds"),
        gsea_wp_obj_file_path = os.path.join("output", "salmon", "{condition}", "gsea_analysis", "gsea_wp_obj.rds"),
        gsea_go_summary_file_path = os.path.join("output", "salmon", "{condition}", "gsea_analysis", "gsea_go_summary.csv"),
        gsea_kegg_summary_file_path = os.path.join("output", "salmon", "{condition}", "gsea_analysis", "gsea_kegg_summary.csv"),
        gsea_wp_summary_file_path = os.path.join("output", "salmon", "{condition}", "gsea_analysis", "gsea_wp_summary.csv"),

        # DotPlots
        dotplot_gsea_go_file_path = os.path.join("output", "salmon", "{condition}", "gsea_analysis", "plots", "gsea_go-dotplot.jpg"),
        dotplot_gsea_kegg_file_path = os.path.join("output", "salmon", "{condition}", "gsea_analysis", "plots", "gsea_kegg-dotplot.jpg"),
        dotplot_gsea_wp_file_path = os.path.join("output", "salmon", "{condition}", "gsea_analysis", "plots", "gsea_wp-dotplot.jpg"),
        # GSEA-plots top 10
        gsea_go_top10_plot_file_path_1=os.path.join("output", "salmon", "{condition}", "gsea_analysis", "plots", "gsea_go-top10_plot_1.jpg"),
        gsea_kegg_top10_plot_file_path_1=os.path.join("output", "salmon", "{condition}", "gsea_analysis", "plots", "gsea_kegg-top10_plot_1.jpg"),
        gsea_wp_top10_plot_file_path_1=os.path.join("output", "salmon", "{condition}", "gsea_analysis", "plots", "gsea_wp-top10_plot_1.jpg"),
        gsea_go_top10_plot_file_path_2=os.path.join("output", "salmon", "{condition}","gsea_analysis","plots","gsea_go-top10_plot_2.jpg"),
        gsea_kegg_top10_plot_file_path_2=os.path.join("output", "salmon", "{condition}","gsea_analysis","plots","gsea_kegg-top10_plot_2.jpg"),
        gsea_wp_top10_plot_file_path_2=os.path.join("output", "salmon", "{condition}","gsea_analysis","plots","gsea_wp-top10_plot_2.jpg"),
    params:
        input_algorithm="salmon"
    log:
        "logs/salmon_run_gsea_analysis/{condition}.log"
    threads:
        1
    resources:
        mem_mb=4096*8*2,         # total: assign 4096MB per CPU
        cpus=1,
        time_min=3*60		    # use 3 hrs, to make sure
    conda:
        "../envs/deseq2_gsea_env.yaml"
    script:
        "../scripts/deseq2/gene_set_enrichment_analysis.R"


rule clusterProfiler_salmon_create_gsea_html_report_files:
    input:
        os.path.join("output", "salmon", "{condition}", "gsea_analysis", "gsea_go_summary.csv"),
        os.path.join("output", "salmon", "{condition}", "gsea_analysis", "gsea_kegg_summary.csv"),
        os.path.join("output", "salmon", "{condition}", "gsea_analysis", "gsea_wp_summary.csv"),
    output:
        report(
            directory(os.path.join("output", "salmon", "{condition}", "gsea_analysis", "html_reports")),
            patterns=["gsea_{used_db}_{filename}.csv.report.html"],
            # Path relative to Snakefile
            caption="../../../report_source/module3_gene_expression_quantification_and_analysis/gsea/caption.rst",
            category=config["output_dir_module3_gene_expression_quantification_and_analysis"].replace("output/",""),
            subcategory="Salmon: Gene Set Enrichment Analysis (GSEA)",
            labels={  # Set labels manually
                "Tool": "Salmon-ClusterProfiler",
                "Database": "{used_db}",
                "Condition": "{condition}",
                "Type": "Report",
                "File:": "gsea_{used_db}_{filename}",
            }
        ),
        os.path.join("output", "salmon", "{condition}", "gsea_analysis", "html_reports", "gsea_go_summary.csv.report.html"),
        os.path.join("output", "salmon", "{condition}", "gsea_analysis", "html_reports", "gsea_kegg_summary.csv.report.html"),
        os.path.join("output", "salmon", "{condition}", "gsea_analysis", "html_reports", "gsea_wp_summary.csv.report.html"),
    params:
        input_files=lambda wildcards, input: input,
        data_separators=["," for i in range(3)],
        data_titles=lambda wildcards: [
            f"Gene Set Enrichment Analysis (clusterProfiler) for {wildcards.condition}" for i in range(3)],
        info_texts=[
            workflow.source_path(
                "../../../report_source/module3_gene_expression_quantification_and_analysis/gsea/gsea_salmon_info.html")
            for i in range(3)
        ],
        html_output_dir=lambda wildcards, output: output[0],
        html_output_file_basenames=lambda wildcards, input: [
            (os.path.basename(input_file) + ".report.html") for input_file in input
        ]
    conda:
        "../../../report_source/report_env.yaml"
    script:
        "../../../scripts/create_report_html_files.R"


# ============================= Kallisto: GSEA analysis =============================
rule clusterProfiler_kallisto_deseq2_gene_set_enrichment_analysis:
    input:
        deseq_dataset_obj=os.path.join("output", "kallisto", "{condition}", "kallisto_deseq_dataset_object.rds")
    output:
        report(
            directory(os.path.join("output", "kallisto", "{condition}", "gsea_analysis", "plots")),
            patterns=["gsea_{used_db}-{filename}.jpg"],
            # Path relative to Snakefile
            caption="../../../report_source/module3_gene_expression_quantification_and_analysis/gsea/caption.rst",
            category=config["output_dir_module3_gene_expression_quantification_and_analysis"].replace("output/", ""),
            subcategory="Kallisto: Gene Set Enrichment Analysis (GSEA)",
            labels={  # Set labels manually
                "Tool": "Kallisto-ClusterProfiler",
                "Database": "{used_db}",
                "Condition": "{condition}",
                "Type": "Plot",
                "File:": "gsea_{used_db}-{filename}.jpg",
            }
        ),
        gsea_go_obj_file_path = os.path.join("output", "kallisto", "{condition}", "gsea_analysis", "gsea_go_obj.rds"),
        gsea_kegg_obj_file_path = os.path.join("output", "kallisto", "{condition}", "gsea_analysis", "gsea_kegg_obj.rds"),
        gsea_wp_obj_file_path = os.path.join("output", "kallisto", "{condition}", "gsea_analysis", "gsea_wp_obj.rds"),
        gsea_go_summary_file_path = os.path.join("output", "kallisto", "{condition}", "gsea_analysis", "gsea_go_summary.csv"),
        gsea_kegg_summary_file_path = os.path.join("output", "kallisto", "{condition}", "gsea_analysis", "gsea_kegg_summary.csv"),
        gsea_wp_summary_file_path = os.path.join("output", "kallisto", "{condition}", "gsea_analysis", "gsea_wp_summary.csv"),

        # DotPlots
        dotplot_gsea_go_file_path = os.path.join("output", "kallisto", "{condition}", "gsea_analysis", "plots", "gsea_go-dotplot.jpg"),
        dotplot_gsea_kegg_file_path = os.path.join("output", "kallisto", "{condition}", "gsea_analysis", "plots", "gsea_kegg-dotplot.jpg"),
        dotplot_gsea_wp_file_path = os.path.join("output", "kallisto", "{condition}", "gsea_analysis", "plots", "gsea_wp-dotplot.jpg"),
        # GSEA-plots top 10
        gsea_go_top10_plot_file_path_1=os.path.join("output", "kallisto", "{condition}", "gsea_analysis", "plots", "gsea_go-top10_plot_1.jpg"),
        gsea_kegg_top10_plot_file_path_1=os.path.join("output", "kallisto", "{condition}", "gsea_analysis", "plots", "gsea_kegg-top10_plot_1.jpg"),
        gsea_wp_top10_plot_file_path_1=os.path.join("output", "kallisto", "{condition}", "gsea_analysis", "plots", "gsea_wp-top10_plot_1.jpg"),
        gsea_go_top10_plot_file_path_2=os.path.join("output", "kallisto", "{condition}","gsea_analysis","plots","gsea_go-top10_plot_2.jpg"),
        gsea_kegg_top10_plot_file_path_2=os.path.join("output", "kallisto", "{condition}","gsea_analysis","plots","gsea_kegg-top10_plot_2.jpg"),
        gsea_wp_top10_plot_file_path_2=os.path.join("output", "kallisto", "{condition}","gsea_analysis","plots","gsea_wp-top10_plot_2.jpg"),
    params:
        input_algorithm="kallisto"
    log:
        "logs/kallisto_run_gsea_analysis/{condition}.log"
    threads:
        1
    resources:
        mem_mb=4096*8*2,         # total: assign 4096MB per CPU
        cpus=1,
        time_min=3*60		    # use 3 hrs, to make sure
    conda:
        "../envs/deseq2_gsea_env.yaml"
    script:
        "../scripts/deseq2/gene_set_enrichment_analysis.R"


rule clusterProfiler_kallisto_create_gsea_html_report_files:
    input:
        os.path.join("output", "kallisto", "{condition}", "gsea_analysis", "gsea_go_summary.csv"),
        os.path.join("output", "kallisto", "{condition}", "gsea_analysis", "gsea_kegg_summary.csv"),
        os.path.join("output", "kallisto", "{condition}", "gsea_analysis", "gsea_wp_summary.csv"),
    output:
        report(
            directory(os.path.join("output", "kallisto", "{condition}", "gsea_analysis", "html_reports")),
            patterns=["gsea_{used_db}_{filename}.csv.report.html"],
            # Path relative to Snakefile
            caption="../../../report_source/module3_gene_expression_quantification_and_analysis/gsea/caption.rst",
            category=config["output_dir_module3_gene_expression_quantification_and_analysis"].replace("output/",""),
            subcategory="Kallisto: Gene Set Enrichment Analysis (GSEA)",
            labels={  # Set labels manually
                "Tool": "Kallisto-ClusterProfiler",
                "Database": "{used_db}",
                "Condition": "{condition}",
                "Type": "Report",
                "File:": "gsea_{used_db}_{filename}",
            }
        ),
        os.path.join("output", "kallisto", "{condition}", "gsea_analysis", "html_reports", "gsea_go_summary.csv.report.html"),
        os.path.join("output", "kallisto", "{condition}", "gsea_analysis", "html_reports", "gsea_kegg_summary.csv.report.html"),
        os.path.join("output", "kallisto", "{condition}", "gsea_analysis", "html_reports", "gsea_wp_summary.csv.report.html"),
    params:
        input_files=lambda wildcards, input: input,
        data_separators=["," for i in range(3)],
        data_titles=lambda wildcards: [
            f"Gene Set Enrichment Analysis (clusterProfiler) for {wildcards.condition}" for i in range(3)],
        info_texts=[
            workflow.source_path(
                "../../../report_source/module3_gene_expression_quantification_and_analysis/gsea/gsea_kallisto_info.html")
            for i in range(3)
        ],
        html_output_dir=lambda wildcards, output: output[0],
        html_output_file_basenames=lambda wildcards, input: [
            (os.path.basename(input_file) + ".report.html") for input_file in input
        ]
    conda:
        "../../../report_source/report_env.yaml"
    script:
        "../../../scripts/create_report_html_files.R"

