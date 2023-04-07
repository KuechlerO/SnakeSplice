# Detection of rare aberrant splicing events in transcriptome profiles.
# The workflow aims to assist the diagnostics in the field of rare diseases where RNA-seq is performed to
# identify aberrant splicing defects.

import os

fraser_output_path = os.path.join('output', config["output_directories"]["fraser_output_dir"])


rule fraser_create_annotation_file:
    """
    Create a sample annotation file for Fraser.
    Requires a "bamFile"-column, which contains the paths to the respective bam files.
    """
    input:
        bam_file_paths =
            expand(main_helper_get_bam_file_for_sample(config["input_dir_of_bam_files"],
                    filename_extension=config["bam_files_attributes"]["filename_extension"]),
                sample_id=pep.sample_table["sample_name"])
    output:
        annotation_table=os.path.join(fraser_output_path, 'sample_annotation_file.txt'),
    threads:
        1
    resources:
        mem_mb=4 * 4096,# total: assign 4096MB per CPU
        cpus=1,# uses 1 cpus
        time_min=300  # give max 5 hrs
    run:
        annotation_table = pep.sample_table.copy()
        # Add the bam file paths to the annotation table in column "bamFile"
        annotation_table["bamFile"] = input.bam_file_paths
        annotation_table["sampleID"] = annotation_table["sample_name"]
        annotation_table["pairedEnd"] = "TRUE"

        # Drop all columns except "sampleID", "bamFile" and "pairedEnd"
        annotation_table = annotation_table[["sampleID", "bamFile", "pairedEnd"]]

        # Save output table
        annotation_table.to_csv(output.annotation_table, sep="\t", index=False)


rule fraser_exploration_and_creation_of_fraser_dataset_object:
    """
    Run Fraser to detect differentially spliced genes.
    Utilizes the BAM-files annotated in the sample annotation file to derive the differential splicing events.
    """
    input:
        sample_annotation_file=rules.fraser_create_annotation_file.output.annotation_table
    output:
        report(
            directory(os.path.join(fraser_output_path, "fraser_exploration_plots")),
            patterns=["{filename}.jpg"],
            # Path relative to Snakefile
            caption="../../../report_source/module5_splicing_patterns_analysis/fraser/caption.rst",
            category=config["output_dir_module5_splicing_patterns_analysis"].replace("output/", ""),
            subcategory="FRASER",
            labels={  # Set labels manually
                "Tool": "FRASER",
                "Mode": "Exploration",
                "Type": "Plot",
                "File:": "{filename}"
            }
        ),
        # Plot files: After filtering, no normalization
        plot_filter_expression_file=os.path.join(fraser_output_path, "fraser_exploration_plots",
            'fraser_plot_filter_expression.jpg'),
        plot_cor_psi5_heatmap_file=os.path.join(fraser_output_path, "fraser_exploration_plots",
            'fraser_plot_cor_psi5_heatmap.jpg'),
        plot_cor_psi3_heatmap_file=os.path.join(fraser_output_path, "fraser_exploration_plots",
            'fraser_plot_cor_psi3_heatmap.jpg'),
        plot_cor_theta_heatmap_file=os.path.join(fraser_output_path, "fraser_exploration_plots",
            'fraser_plot_cor_theta_heatmap.jpg'),

        # Top 100 plots
        # plot_cor_psi5_top100_heatmap_file=os.path.join(fraser_output_path,"fraser_exploration_plots",
        #     'fraser_plot_cor_psi5_top100_heatmap.jpg'),
        # plot_cor_psi3_top100_heatmap_file=os.path.join(fraser_output_path,"fraser_exploration_plots",
        #     'fraser_plot_cor_psi3_top100_heatmap.jpg'),
        # plot_cor_theta_top100_heatmap_file=os.path.join(fraser_output_path, "fraser_exploration_plots",
        #     'fraser_plot_cor_theta_top100_heatmap.jpg'),

        # Plot files: After filtering, with normalization
        plot_normalized_cor_psi5_heatmap_file=os.path.join(fraser_output_path, "fraser_exploration_plots",
            'fraser_plot_normalized_cor_psi5_heatmap.jpg'),
        plot_normalized_cor_psi3_heatmap_file=os.path.join(fraser_output_path, "fraser_exploration_plots",
            'fraser_plot_normalized_cor_psi3_heatmap.jpg'),
        plot_normalized_cor_theta_heatmap_file=os.path.join(fraser_output_path, "fraser_exploration_plots",
            'fraser_plot_normalized_cor_theta_heatmap.jpg'),

        fraser_data_set_object_file=directory(os.path.join(fraser_output_path,
            "savedObjects", 'fraser_analysis_set_object.RData')),
    log:
        "logs/fraser/fraser_diff_splicing_analysis.log"
    conda:
        "../envs/fraser_env.yaml"
    threads:
        32
    resources:
        mem_mb=32 * 4096,       # total: assign 4096MB per CPU
        cpus=32,                # uses 32 cpus
        time_min=300            # give max 5 hrs
    script:
        "../scripts/fraser/fraser_dataset_exploration.R"


rule fraser_create_fraser_analysis_plots:
    """
    Create plots from the Fraser analysis set object.
    -> 1. Plot: Aberrant splicing events
    -> 2. Plot: QQ-Plot
    """
    input:
        fraser_analysis_set_object_file = rules.fraser_exploration_and_creation_of_fraser_dataset_object.output.fraser_data_set_object_file
    output:
        report(
            directory(os.path.join(fraser_output_path,"fraser_analysis_plots")),
            patterns=["{filename}.jpg"],
            # Path relative to Snakefile
            caption="../../../report_source/module5_splicing_patterns_analysis/fraser/caption.rst",
            category=config["output_dir_module5_splicing_patterns_analysis"].replace("output/", ""),
            subcategory="FRASER",
            labels={  # Set labels manually
                "Tool": "FRASER",
                "Mode": "Analysis",
                "Type": "Plot",
                "File:": "{filename}"
            }
        ),
        # Output files for differential splicing analysis
        csv_summary_table_file=os.path.join(fraser_output_path,'fraser_csv_summary_table.csv'),
        plot_aberrant_events_per_sample_file=os.path.join(fraser_output_path, "fraser_analysis_plots",
            'fraser_plot_aberrant_events.jpg'),
        plot_qq_plot_file=os.path.join(fraser_output_path, "fraser_analysis_plots",
            'fraser_plot_qq_plot.jpg'),
    log:
        "logs/fraser/create_fraser_analysis_plots.log"
    conda:
        "../envs/fraser_env.yaml"
    threads:
        1   # uses 1 cpus
    resources:
        mem_mb=4 * 4096,        # total: assign 4096MB per CPU
        cpus=1,                 # uses 1 cpus
        time_min=20             # give max 5 hrs
    script:
        "../scripts/fraser/create_fraser_analysis_plots.R"


rule fraser_create_html_reports:
    input:
        rules.fraser_create_fraser_analysis_plots.output.csv_summary_table_file,
    output:
        report(
            directory(os.path.join(fraser_output_path, "html_reports")),
            patterns=["{filename}.csv.report.html"],
            # Path relative to Snakefile
            caption="../../../report_source/module5_splicing_patterns_analysis/fraser/caption.rst",
            category=config["output_dir_module5_splicing_patterns_analysis"].replace("output/",""),
            subcategory="FRASER",
            labels={  # Set labels manually
                "Tool": "FRASER",
                "Mode": "Analysis",
                "Type": "Report",
                "File:": "{filename}"
            }
        ),
        # Just two output examples... -> Total of 26 files will be generated
        os.path.join(fraser_output_path, "html_reports", "fraser_csv_summary_table.csv.report.html"),
    params:
        input_files=lambda wildcards,input: input,
        data_separators=["," for i in range(1)],
        data_titles=["FRASER results" for i in range(1)],
        info_texts=[workflow.source_path("../../../report_source/module5_splicing_patterns_analysis/fraser/info.html")],
        html_output_dir=lambda wildcards, output: output[0],
        html_output_file_basenames=lambda wildcards, input: [
            (os.path.basename(input_file) + ".report.html") for input_file in input
        ]
    conda:
        "../../../report_source/report_env.yaml"
    script:
        "../../../scripts/create_report_html_files.R"
