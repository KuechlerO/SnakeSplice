# ========== Submodule: Checking strandedness of reads ==========
# NOTE: Personal Experience: This tool in combination with Snakemake
#      is not very stable. It is recommended to run this tool repeatedly in
#      case of errors.

import os

cos_output_path = os.path.join('output', config["output_directories"]["check_strandedness_output_dir"])

ensembl_hs_ch38_fasta_link = "https://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
ensembl_hs_ch38_gtf_link = "https://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz"


# 1. Download transcriptome data as reference data
rule check_strandedness_download_transcript_reference_data:
    """
    Download transcriptome data as reference data
    """
    output:
        output_file=os.path.join(cos_output_path, "reference_files",
            ensembl_hs_ch38_fasta_link.split("/")[-1].replace(".gz", ""))
    params:
        output_dir=lambda wildcards, output: os.path.dirname(output[0]),
        ensembl_all_cdna_fasta_file_link=ensembl_hs_ch38_fasta_link,
        gz_file=ensembl_hs_ch38_fasta_link.split("/")[-1]
    shell:
        "mkdir -p {params.output_dir} && "
        "wget -O {output.output_file}.gz {params.ensembl_all_cdna_fasta_file_link};"
        "cd {params.output_dir};"
        "gunzip {params.gz_file}"


rule check_strandedness_download_annotation_reference_data:
    """
    Download annotation data as reference data
    """
    output:
        output_file=os.path.join(cos_output_path, "reference_files",
            ensembl_hs_ch38_gtf_link.split("/")[-1].replace(".gz",""))
    params:
        output_dir=lambda wildcards, output: os.path.dirname(output[0]),
        ensembl_gtf_file_link=ensembl_hs_ch38_gtf_link,
        gz_file=ensembl_hs_ch38_gtf_link.split("/")[-1]
    shell:
        "mkdir -p {params.output_dir} && "
        "wget -O {output.output_file}.gz {params.ensembl_gtf_file_link} && "
        "cd {params.output_dir};"
        "gunzip {params.gz_file}"


# 2. Run tool
rule check_strandedness_run_strandedness_check:
    """
    Run the strandedness check tool on the sample directory (for every pair of reads)
    """
    input:
        unpack(get_paired_read_file_paths),        # Defined in rules/lib_get_read_file_paths.smk, imported in Snakefile
        transcripts_file=rules.check_strandedness_download_transcript_reference_data.output.output_file,
        annotation_file=rules.check_strandedness_download_annotation_reference_data.output.output_file
    output:
        os.path.join(cos_output_path, "{sample_id}.txt")
    log:
        "logs/check_strandedness/{sample_id}.log"
    conda:
        "../envs/check_strandedness_env.yaml"
    resources:
        mem_mb=32000,	    # total: 32GB
        cpus=2				# uses 2 cpus
    shell:
        "check_strandedness "
        "--transcripts {input.transcripts_file} "
        "--gtf {input.annotation_file} "
        "--reads_1 {input.r1} --reads_2 {input.r2} > {output[0]} 2> {log};"


# 3. Summarize output
rule check_strandedness_summarize_strandedness_check:
    """
    Summarize strandedness check results
    """
    input:
        expand(os.path.join(cos_output_path, "{sample_id}.txt"), sample_id=pep.sample_table["sample_name"])
    output:
        os.path.join(cos_output_path, "summary_report.csv")
    run:
        # Iterate over all input files and read lines 9-13 and put them into a pandas table
        # Then save table into CSV-file
        import pandas as pd
        col0_key = "sample"
        col1_key = "Fraction of reads failed to determine strandedness"
        col2_key = "Fraction of reads explained by FR"
        col3_key = "Fraction of reads explained by RF"
        col4_key = "Summary"
        col5_key = "Conclusion"
        col6_key = "stranded-value"
        table = pd.DataFrame(columns=[col0_key, col1_key, col2_key, col3_key, col4_key, col5_key])

        for sample_file in input:
            sample_name = sample_file.split("/")[-1].split(".")[0]
            with open(sample_file, "r") as f:
                lines = f.readlines()

                # Add row to table
                row = pd.DataFrame([[sample_name, lines[-5].split(":")[1].strip(), lines[-4].split(":")[1].strip(),
                                     lines[-3].split(":")[1].strip(), lines[-2].strip(),
                                     lines[-1].strip()]],
                    columns=[col0_key, col1_key, col2_key, col3_key, col4_key, col5_key])
                table = table.append(row, ignore_index=True)

        # Add column with annotations for the sample configuration file
        table[col6_key] = "ERROR: strandedness could not be determined"
        # "no" for unstranded data
        table.loc[table[col5_key].str.contains("unstranded"), col6_key] = "no"
        # If the conclusion contains "RF/fr-firststrand" then set the value to "reverse", otherwise to "yes"
        table.loc[table[col5_key].str.contains("RF/fr-firststrand"), col6_key] = "reverse"
        table.loc[table[col5_key].str.contains("FR/fr-secondstrand"), col6_key] = "yes"

        table.to_csv(output[0], index=False)


rule check_strandedness_create_report:
    """
    Create HTML report
    """
    input:
        rules.check_strandedness_summarize_strandedness_check.output[0]
    output:
        report(
            os.path.join(cos_output_path, "html_reports", "summary_report.html"),
            # Path relative to Snakefile
            caption="../../../report_source/module1_qc_preproc_alignment/check_strandedness/caption.rst",
            category=config["output_dir_module1_qc_preprocessing_and_alignment"].replace("output/", ""),
            subcategory="Strandedness Check",
        ),
    params:
        input_files= lambda wildcards, input: [input[0]],
        data_separators=[","],
        data_titles=["Check Strandedness: Summary Report"],
        info_texts=[workflow.source_path("../../../report_source/module1_qc_preproc_alignment/check_strandedness/info.html")],
        html_output_dir=lambda wildcards, output: os.path.dirname(output[0]),
        html_output_file_basenames=lambda wildcards, output: [os.path.basename(output[0])]
    conda:
        "../../../report_source/report_env.yaml"
    script:
        "../../../scripts/create_report_html_files.R"


rule check_strandedness_clean_up:
    """
    Clean up:
    - Remove all files starting with "stranded_test_" in the Snakemake_Main directory
    - Deletes kallisto_index
    """
    input:
        rules.check_strandedness_summarize_strandedness_check.output[0]
    output:
        touch(os.path.join(cos_output_path, "clean_up.done"))
    shell:
        "rm -rf ./stranded_test_*; "
        "rm -rf ./kallisto_index; "
