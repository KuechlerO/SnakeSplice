# ========== Submodule: Checking stats of BAM-file ==========

import os

# Output directory
bamtools_output_path = os.path.join('output',config["output_directories"]["bamstats_output_dir"])


rule bamtools_stats_for_olego_output:
    input:
        rules.olego_samtools_sort.output		    # take sorted BAM-files as input
    output:
        report(
            os.path.join(bamtools_output_path,"olego","{sample_id}.bamstats"),
            caption="../../../report_source/module1_qc_preproc_alignment/bamtools/caption.rst",
            category=config["output_dir_module1_qc_preprocessing_and_alignment"].replace("output/",""),
            subcategory="Samtools Statistics",
            labels={  # Set labels manually
                "Aligner": "Olego",
                "File:": "{sample_id}.bamstats"
            }
        )
    params:
        ""
    log:
        "logs/bamtools/stats/olego_{sample_id}.log"
    wrapper:
        "v1.21.4/bio/bamtools/stats"


rule bamtools_stats_for_star_output:
    input:
        rules.star_samtools_sort.output		    # take sorted BAM-files as input
    output:
        os.path.join(bamtools_output_path, "star", "{sample_id}.bamstats"),
    params:
        ""
    log:
        "logs/bamtools/stats/star_{sample_id}.log"
    threads:
        1
    wrapper:
        "v1.21.4/bio/bamtools/stats"


# TODO also include percentage numbers?!
rule bamtools_stats_extract_overview:
    input:
        os.path.join(bamtools_output_path, "{aligner_tool}", "{sample_id}.bamstats"),
    output:
        os.path.join(bamtools_output_path,"{aligner_tool}","{sample_id}.bamstats.extracted.tsv"),
    log:
        "logs/bamtools/stats/extract_{aligner_tool}_{sample_id}.log"
    shell:
        "tail -n 13 {input} | head -n 12 | cut -f 1 > {output} 2> {log}"


rule bamtools_stats_merge_all_stats:
    input:
        expand(os.path.join(bamtools_output_path, "{{aligner_tool}}", "{sample_id}.bamstats.extracted.tsv"),
            sample_id=pep.sample_table["sample_name"]),
    output:
        os.path.join(bamtools_output_path, "{aligner_tool}", "merged_bamstats.tsv"),
    threads:
        1
    run:
        import pandas as pd

        # Open all TSV-files of input, transpose them and merge them into one dataframe
        # Keep filename as row name
        df_list = []
        for file in input:
            current_df = pd.read_csv(file, sep=":\s*", engine="python", header=None, index_col=0)
            current_df = current_df.transpose()
            current_df["BAM-file"] = os.path.basename(file)
            df_list.append(current_df)

        # Finalize column order and output merged dataframe
        output_df = pd.concat(df_list, axis=0)
        cols = output_df.columns.tolist()
        final_cols = cols[-1:] + cols[:-1]
        output_df = output_df[final_cols]
        output_df.to_csv(output[0], sep="\t", index=False)


rule bamtools_create_html_reports:
    input:
        rules.bamtools_stats_merge_all_stats.output
    output:
        report(
            directory(os.path.join(bamtools_output_path, "{aligner_tool}", "html_reports")),
            patterns=["{filename}.report.html"],
            # Path relative to Snakefile
            caption="../../../report_source/module1_qc_preproc_alignment/bamtools/caption.rst",
            category=config["output_dir_module1_qc_preprocessing_and_alignment"].replace("output/",""),
            subcategory="Samtools Statistics",
            labels={  # Set labels manually
                "Aligner": "{aligner_tool}",
                "File:": "{filename}"
            }
        ),
        os.path.join(bamtools_output_path,  "{aligner_tool}", "html_reports", "merged_bamstats.tsv.report.html"),
    params:
        input_files = lambda wildcards, input: [input[0]],
        data_separators=["\t"],
        data_titles = ["Samtools stats output"],
        info_texts=[workflow.source_path("../../../report_source/module1_qc_preproc_alignment/bamtools/info.html")],
        html_output_dir=lambda wildcards, output: output[0],
        html_output_file_basenames = lambda wildcards, output: [os.path.basename(output[1])
        ]
    conda:
        "../../../report_source/report_env.yaml"
    script:
        "../../../scripts/create_report_html_files.R"
