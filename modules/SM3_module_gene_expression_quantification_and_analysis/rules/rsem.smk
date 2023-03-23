import os

# Output directory
rsem_output_path = os.path.join("output", config["output_directories"]["rsem_output_dir"])


rule rsem_gunzip_input_transcripts:
    """Unzip input transcripts"""
    input:
        config["rsem_settings"]["transcriptome_fasta"],
    output:
        temp(config["rsem_settings"]["transcriptome_fasta"].replace(".gz", "")),
    shell:
        "gunzip --keep --to-stdout {input} > {output}"


rule rsem_prepare_reference:
    """
    Creates index files for downstream analysis with RSEM
    """
    input:
        # reference FASTA with either the entire genome or transcript sequences
        reference_genome=config["rsem_settings"]["transcriptome_fasta"].replace(".gz", ""),
    output:
        # one of the index files created and used by RSEM (required)
        seq=os.path.join(rsem_output_path, "index", "reference.seq"),
        # RSEM produces a number of other files which may optionally be specified as output;
        # these may be provided so that snakemake is aware of them, but the wrapper doesn't do anything with this information other than to verify that the file path prefixes match that of output.seq.
        # for example,
        grp=os.path.join(rsem_output_path, "index", "reference.grp"),
        ti=os.path.join(rsem_output_path, "index", "reference.ti"),
    log:
        "logs/rsem/prepare-reference.log",
    params:
        # optional additional parameters, for example,
        #extra="--gtf annotations.gtf",
        # if building the index against a reference transcript set
        extra="",
    threads:
        32
    resources:
        mem_mb=32*4096,         # total: assign 4096MB per CPU = 131072 MB, 131GB
        cpus=32,		        # uses 32 cpus
        time_min=6*60		    # use 6 hrs, to make sure
    wrapper:
        "v1.21.3-3-gcb96cc40/bio/rsem/prepare-reference"


def get_reads_fasta_for_sample(wildcards):
    """
    Returns for a given sample (wildcards.sample_id) the path to the fasta file
    -> Read1-File
    :param wildcards:
    :return:
    """
    r1 = os.path.join(
        pep.sample_table[pep.sample_table["sample_name"] == wildcards.sample_id]["sample_directory"][0],
        pep.sample_table[pep.sample_table["sample_name"] == wildcards.sample_id]["read1"][0]
    )
    r2 = os.path.join(
        pep.sample_table[pep.sample_table["sample_name"] == wildcards.sample_id]["sample_directory"][0],
        pep.sample_table[pep.sample_table["sample_name"] == wildcards.sample_id]["read2"][0]
    )
    return {"fq_one": r1, "fq_two": r2}


rule rsem_create_bowtie_index_of_reference_transcriptome:
    input:
        reference_transcriptome=rules.rsem_gunzip_input_transcripts.output,
    output:
        multiext(os.path.join(rsem_output_path, "index/reference"),
            ".1.ebwt",".2.ebwt",".3.ebwt",".4.ebwt",".rev.1.ebwt",".rev.2.ebwt")
    log:
        "logs/rsem/create_bowtie_index_of_reference_transcriptome.log"
    params:
        output_base_name = lambda wildcards, output: output[0].split(".1.ebwt")[0]
    conda:
        "../envs/rsem_env.yaml"
    shell:
        "bowtie-build {input} {params.output_base_name} 2> {log}"


rule rsem_calculate_expression_fasta:
    input:
        # input.bam or input.fq_one must be specified (and if input.fq_one, optionally input.fq_two if paired-end)
        # an aligned to transcriptome BAM
        # bam="mapped/a.bam",
        unpack(get_reads_fasta_for_sample),     # return fq_one & fq_two
        # one of the index files created by rsem-prepare-reference; the file suffix is stripped and passed on to rsem
        reference=rules.rsem_prepare_reference.output.seq,
        # reference_bowtie: Additionally needed for FASTQ input; Index files created (by bowtie-build) from the reference transcriptome
        reference_bowtie=multiext(os.path.join(rsem_output_path, "index/reference"),
            ".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt"),
    output:
        # genes_results must end in .genes.results; this suffix is stripped and passed to rsem as an output name prefix
        # this file contains per-gene quantification data for the sample
        genes_results=os.path.join(rsem_output_path, "quantification", "quant_results_{sample_id}/{sample_id}.genes.results"),
        # isoforms_results must end in .isoforms.results and otherwise have the same prefix as genes_results
        # this file contains per-transcript quantification data for the sample
        isoforms_results=os.path.join(rsem_output_path, "quantification", "quant_results_{sample_id}/{sample_id}.isoforms.results"),
    log:
        "logs/rsem/calculate_expression/{sample_id}.log",
    params:
        # optional, specify if sequencing is paired-end
        paired_end=True,
        # additional optional parameters to pass to rsem, for example,
        extra="--seed 42 "
    threads: 32
    resources:
        mem_mb=32 * 4096,# total: assign 4096MB per CPU = 131072 MB, 131GB
        cpus=32,	# uses 32 cpus
        time_min= 6 *60		    # use 6 hrs, to make sure
    wrapper:
        "v1.21.3-3-gcb96cc40/bio/rsem/calculate-expression"


rule rsem_generate_gene_data_matrix:
    """
    Generate a gene data matrix from RSEM output
    -> Combine all genes.results/isoforms.results files into one matrix
    """
    input:
        # one or more expression files created by rsem-calculate-expression
        expand(os.path.join(rsem_output_path, "quantification",
            "quant_results_{sample_id}/{sample_id}.{{genes_or_isoforms}}.results"),
               sample_id=pep.sample_table["sample_name"]),
    output:
        # a tsv containing each sample in the input as a column
        os.path.join(rsem_output_path, "quantification", "{genes_or_isoforms}_total_expression_matrix.tsv"),
    params:
        # optional additional parameters
        extra="",
    log:
        "logs/rsem/generate_{genes_or_isoforms}_data_matrix.log",
    wrapper:
        "v1.20.0/bio/rsem/generate-data-matrix"


rule rsem_create_output_annotation_table:
    """
        Creates an annotation table that is required for the analysis of the Salmon output
        via the R-script
    """
    input:
        # sample.genes.results as Input
        quant_files=expand(os.path.join(rsem_output_path, "quantification",
            "quant_results_{sample_id}/{sample_id}.genes.results"), sample_id=pep.sample_table["sample_name"]),
    output:
        # Annotation table -> copying PEP-annotation table
        annotation_table=os.path.join(rsem_output_path, "{condition}", "rsem_output_annotation_table.tsv")
    params:
        control_samples=pep.sample_table[(pep.sample_table["control"] == "true")]["sample_name"],
        condition_samples=main_helper_extract_sample_ids_with_specific_condition
    threads:
        1
    run:
        all_controls_table = pep.sample_table[pep.sample_table["sample_name"].isin(params.control_samples)]
        all_condition_table = pep.sample_table[pep.sample_table["sample_name"].isin(params.condition_samples)]
        # Create a table with all samples
        all_samples_table = all_controls_table.append(all_condition_table)

        # Add the quant.sf files to the table
        results_dir_path = os.path.dirname(os.path.dirname(input.quant_files[0]))
        quant_file_path = os.path.join(results_dir_path, "quant_results_{sample_id}", "quant.sf")
        # For each sample, replace the sample name with the quant.sf file path in the all_samples_table
        all_samples_table["salmon_results_file"] = \
            all_samples_table.apply(lambda row: quant_file_path.replace("{sample_id}", row["sample_name"]), axis=1)

        # Save output table
        all_samples_table.to_csv(output.annotation_table, sep="\t", index=False)


rule rsem_create_deseq_dataset_obj:
    """
    Creates a DESeq2 object from the Kallisto output
    """
    input:
        annotation_table_file=os.path.join(rsem_output_path, "{condition}", "rsem_output_annotation_table.tsv")
    output:
        deseq_dataset_r_obj=os.path.join(rsem_output_path, "{condition}", "rsem_deseq_dataset_object.rds")
    params:
        count_algorithm="kallisto",
    conda:
        "../envs/deseq2_gene_level_env.yaml"
    threads:
        12
    resources:
        mem_mb=12 * 4096,       # total: assign 4096MB per CPU
        cpus=12,	            # uses 12 cpus
        time_min=6*60		    # use 6 hrs, to make sure
    script:
        "../scripts/deseq2/create_deseq_dataset_object.R"


rule rsem_explore_deseq_dataset_obj:
    """
    Explores the deseq dataset object
    -> Creates 2 Heatmaps & 3 PCA plots:
    1. Euclidian distance
    2. Poisson distance
    3. DESeq2 given plot function for PCAs
    4. Custom plot function for PCAs
    5. Generalized PCA using glmPCA
    """
    input:
        deseq_dataset_r_obj=rules.rsem_create_deseq_dataset_obj.output[0]
    output:
        report(
            directory(os.path.join(rsem_output_path, "{condition}", "explore_deseq_dataset_obj")),
            patterns=["{filename}.jpg"],
            # Path relative to Snakefile
            caption="../../../report_source/module3_gene_expression_quantification_and_analysis/kallisto/caption.rst",
            category=config["output_dir_module3_gene_expression_quantification_and_analysis"].replace("output/", ""),
            subcategory="Kallisto",
            labels={  # Set labels manually
                "Tool": "RSEM",
                "Condition": "{condition}",
                "Modus": "Exploration",
                "File:": "{filename}"
            }
        ),
        os.path.join(rsem_output_path, "{condition}", "explore_deseq_dataset_obj", "sampleEuclidianDistMatrix.jpg"),
        os.path.join(rsem_output_path, "{condition}", "explore_deseq_dataset_obj", "samplePoissonDistMatrix.jpg"),
        os.path.join(rsem_output_path, "{condition}", "explore_deseq_dataset_obj", "PCA_rlog_transformed.jpg"),
        os.path.join(rsem_output_path, "{condition}", "explore_deseq_dataset_obj", "PCA_glmPCA.jpg"),
    log:
        "logs/rsem/explore_deseq_dataset_obj_{condition}.log"
    params:
        output_file_paths = lambda wildcards, output: output[1:]
    conda:
        "../envs/deseq2_gene_level_env.yaml"
    threads:
        12
    resources:
        mem_mb=12*4096,         # total: assign 4096MB per CPU
        cpus=12,	            # uses 12 cpus
        time_min=360		    # use 6 hrs, to make sure
    script:
        "../scripts/deseq2/explore_deseq_dataset.R"


rule rsem_run_gene_level_exp_analysis_with_deseq2:
    input:
        deseq_dataset_r_obj=rules.rsem_create_deseq_dataset_obj.output[0]
    output:
        report(
            directory(os.path.join(rsem_output_path, "{condition}", "deseq2_gene_level_analysis")),
            patterns=["{filename}.jpg"],
            # Path relative to Snakefile
            caption="../../../report_source/module3_gene_expression_quantification_and_analysis/kallisto/caption.rst",
            category=config["output_dir_module3_gene_expression_quantification_and_analysis"].replace("output/",""),
            subcategory="Kallisto",
            labels={  # Set labels manually
                "Tool": "RSEM",
                "Condition": "{condition}",
                "Modus": "DESeq2",
                "File:": "{filename}",
            }
        ),
        os.path.join(rsem_output_path, "{condition}", "deseq2_gene_level_analysis", "statistical_analysis_summary.txt"),
        os.path.join(rsem_output_path, "{condition}", "deseq2_gene_level_analysis", "simple_ma_plot.jpg"),
        os.path.join(rsem_output_path, "{condition}", "deseq2_gene_level_analysis", "apeglm_shrinked_ma_plot.jpg"),
        os.path.join(rsem_output_path, "{condition}", "deseq2_gene_level_analysis", "pvalue_histogram.jpg"),
        os.path.join(rsem_output_path, "{condition}", "deseq2_gene_level_analysis", "pvalue_adjusted_histogram.jpg"),
        os.path.join(rsem_output_path, "{condition}", "deseq2_gene_level_analysis", "gene_cluster_heatmap.jpg"),
        os.path.join(rsem_output_path, "{condition}", "deseq2_gene_level_analysis", "deseq2_results.csv"),
    log:
        "logs/rsem/run_gene_level_exp_analysis_with_deseq2_{condition}.log"
    params:
        output_file_paths = lambda wildcards,output: output[1:],
        used_algorithm= "rsem",
    conda:
        "../envs/deseq2_gene_level_env.yaml"
    threads:
        12
    resources:
        mem_mb=12 * 4096 * 4,# total: assign 4096MB per CPU
        cpus=12,	# uses 12 cpus
        time_min=360		    # use 6 hrs, to make sure
    script:
        "../scripts/deseq2/gene_expression_analysis_with_deseq2.R"


rule rsem_extract_significant_results:
    input:
        deseq2_results=os.path.join(rsem_output_path, "{condition}", "deseq2_gene_level_analysis", "deseq2_results.csv"),
    output:
        os.path.join(rsem_output_path, "{condition}", "deseq2_gene_level_analysis", "deseq2_results_significant.csv"),
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


rule rsem_create_html_reports:
    input:
        rules.rsem_extract_significant_results.output[0],
    output:
        report(
            directory(os.path.join(rsem_output_path,"{condition}","html_reports")),
            patterns=["{filename}.csv.report.html"],
            # Path relative to Snakefile
            caption="../../../report_source/module3_gene_expression_quantification_and_analysis/rsem/caption.rst",
            category=config["output_dir_module3_gene_expression_quantification_and_analysis"].replace("output/",""),
            subcategory="RSEM",
            labels={  # Set labels manually
                "Tool": "RSEM",
                "Condition": "{condition}",
                "Modus": "DESeq2",
                "File:": "{filename}",
            }
        ),
        os.path.join(rsem_output_path,"{condition}", "html_reports", "deseq2_results_significant.csv.report.html"),
    params:
        input_files=lambda wildcards, input: input,
        data_separators=["," for i in range(1)],
        data_titles=["Title needs to be found" for i in range(1)],
        info_texts=["Info-text needs to be found" for i in range(1)],
        html_output_dir=lambda wildcards, output: output[0],
        html_output_file_basenames=lambda wildcards, input: [
            (os.path.basename(input_file) + ".report.html") for input_file in input
        ]
    conda:
        "../../../report_source/report_env.yaml"
    script:
        "../../../scripts/create_report_html_files.R"


# ==================================  Code that is not needed anymore.... ==================================
# rule rsem_convert_bam_files:
#     input:
#         bam_file = main_helper_get_bam_file_for_sample(config["input_dir_of_bam_files"],sample_id=None,
#         filename_extension=config["bam_files_attributes"]["filename_extension"])
#     output:
#         bam_file=os.path.join(rsem_output_path,"converted_bam_files","{sample_id}.bam"),
#     # TODO index might be needed for downstream analysis
#     # bam_file_bai = os.path.join(rsem_output_path, "bam", "{sample_id}.bam.bai")
#     params:
#         output_file_name=lambda wildcards, output: output.bam_file.replace(".bam",""),
#     conda:
#         "../envs/rsem_env.yaml"
#     log:
#         "logs/rsem/convert-bam-files-{sample_id}.log"
#     shell:
#         "convert-sam-for-rsem {input.bam_file} {params.output_file_name} 2> {log}"
#
#
# rule rsem_calculate_expression_bam:
#     input:
#         # input.bam or input.fq_one must be specified (and if input.fq_one, optionally input.fq_two if paired-end)
#         # an aligned to transcriptome BAM
#         bam=rules.rsem_convert_bam_files.output.bam_file,
#         # Index files created by rsem-prepare-reference
#         reference=rules.rsem_prepare_reference.output.seq,
#         # reference_bowtie: Additionally needed for FASTQ input; Index files created (by bowtie-build) from the reference transcriptome
#         # reference_bowtie=multiext("index/reference", ".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt"),
#     output:
#         # genes_results must end in .genes.results; this suffix is stripped and passed to rsem as an output name prefix
#         # this file contains per-gene quantification data for the sample
#         genes_results=os.path.join(rsem_output_path, "quantification", "quant_results_{sample_id}/{sample_id}.genes.results"),
#         # isoforms_results must end in .isoforms.results and otherwise have the same prefix as genes_results
#         # this file contains per-transcript quantification data for the sample
#         isoforms_results=os.path.join(rsem_output_path,"quantification", "quant_results_{sample_id}/{sample_id}.isoforms.results"),
#     log:
#         "logs/rsem/calculate_expression/{sample_id}.log",
#     params:
#         # optional, specify if sequencing is paired-end
#         paired_end=True,
#         # additional optional parameters to pass to rsem, for example,
#         extra="--seed 42 "
#     threads:
#         32
#     resources:
#         mem_mb=32 * 4096,# total: assign 4096MB per CPU = 131072 MB, 131GB
#         cpus=32,	# uses 32 cpus
#         time_min= 6 *60		    # use 6 hrs, to make sure
#     wrapper:
#         "v1.21.3-3-gcb96cc40/bio/rsem/calculate-expression"


