import os

# Output directory
arriba_output_path = os.path.join("output", config["output_directories"]["arriba_output_dir"])


def get_arriba_input_files(wildcards):
    """
    Identify arriba input files for a given sample_id.
    :param wildcards:
    :return:
    """
    given_sample = wildcards.sample_id
    input_bams_dir_choice = config["input_dir_of_bam_files"]
    # Use helper function, which is imported in Snakefile
    sample_file = main_helper_get_bam_file_for_sample(input_bams_dir_choice, given_sample,
                                          filename_extension=config["bam_files_attributes"]["filename_extension"])

    return {
            "bam": sample_file,
            "genome": config["arriba_settings"]["reference_genome_fasta_file"],
            "annotation": config["arriba_settings"]["reference_genome_annotation_file"]}


# Detect gene fusions from chimeric STAR output
rule arriba_detect_gene_fusion_events:
    """
    Detect gene fusions from chimeric STAR output by applying the arriba algorithm.
    """
    input:
        unpack(get_arriba_input_files)  # Input files for arriba
    output:
        # approved gene fusions
        fusions=os.path.join(arriba_output_path, "{sample_id}.fusions.tsv"),
        # discarded gene fusions
        discarded=os.path.join(arriba_output_path, "{sample_id}.fusions.discarded.tsv")
    params:
        # A tsv containing identified artifacts, such as read-through fusions of neighbouring genes, see https://arriba.readthedocs.io/en/latest/input-files/#blacklist
        blacklist=config["arriba_settings"]["arriba_blacklist_file"],
        extra=config["arriba_settings"]["arriba_optional_arguments"]   # Optional arguments for arriba
    log:
        "logs/arriba/{sample_id}.log"
    resources:
        mem_mb=131072,		# total: assign 4096MB per CPU = 131072 MB, 131GB
        cpus=32,			# uses 32 cpus
        time_min=120
    threads: 32
    wrapper:
        "v1.12.0/bio/arriba"


rule arriba_create_summary_table:
    """
    Create a summary table of all detected gene fusions.
    """
    input:
        fusions=expand(os.path.join(arriba_output_path, "{sample_id}.fusions.tsv"),
            sample_id=pep.sample_table["sample_name"]),
    output:
        summary=os.path.join(arriba_output_path, "arriba_gene_fusions_summary.tsv")
    resources:
        mem_mb=4096,		# total: assign 4096MB per CPU = 131072 MB, 131GB
        cpus=1,			# uses 32 cpus
        time_min=120
    threads: 1
    run:
        import pandas as pd
        import os

        # Concat all input fusion files, and add a column with the sample_id as first column
        df = pd.concat([pd.read_csv(f, sep="\t").assign(sample_id=os.path.basename(f).split(".")[0])
                        for f in input.fusions])
        # Place sample_id as first column
        cols = df.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        df = df[cols]

        # Save the summary table
        df.to_csv(output.summary, sep="\t", index=False)


rule arriba_create_html_reports:
    input:
        rules.arriba_create_summary_table.output[0],
    output:
        report(
            directory(os.path.join(arriba_output_path,"html_reports")),
            patterns=["{filename}.report.html"],
            # Path relative to Snakefile
            caption="../../../report_source/module2_gene_fusion_detection/arriba/caption.rst",
            category=config["output_dir_module2_gene_fusion_detection"].replace("output/",""),
            subcategory="Arriba",
            labels={  # Set labels manually
                "Tool": "Arriba",
                "File:": "{filename}"
            }
        ),
        os.path.join(arriba_output_path, "html_reports", "arriba_gene_fusions_summary.tsv.report.html"),
    params:
        input_files = lambda wildcards, input: [input[0]],
        data_separators=["\t", "\t"],
        data_titles = ["Arriba: Detected Fusion Events"],
        info_texts=[
            workflow.source_path("../../../report_source/module2_gene_fusion_detection/arriba/info.html")],
        html_output_dir=lambda wildcards, output: output[0],
        html_output_file_basenames = lambda wildcards, output: [os.path.basename(output[1])
        ]
    conda:
        "../../../report_source/report_env.yaml"
    script:
        "../../../scripts/create_report_html_files.R"