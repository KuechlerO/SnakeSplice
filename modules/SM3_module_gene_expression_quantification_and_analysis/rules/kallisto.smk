import os

# Output directory
kallisto_output_path = os.path.join('output', config["output_directories"]["kallisto_output_dir"])


rule kallisto_index:
    input:
        fasta=config["kallisto_settings"]["transcriptome_fasta"]
    output:
        index=os.path.join(kallisto_output_path, "reference", "transcriptome.idx")
    params:
        extra="",  # optional parameters
    log:
        "logs/kallisto/create_index_transcriptome.log",
    threads: 32
    resources:
        mem_mb=131072,	        # total: assign 4096MB per CPU = 131072 MB, 131GB
        cpus=32,			    # uses 32 cpus
        time_min=360		    # use 6 hrs, to make sure
    wrapper:
        "v1.19.2/bio/kallisto/index"


def get_kallisto_stranded_argument(sample_id):
    """
    Returns the kallisto strandedness argument for a given sample.
    :param sample_id:
    :return:
    """
    strandedness = pep.sample_table[pep.sample_table["sample_name"] == sample_id]["stranded"][0]
    if strandedness == "yes":
        return "--fr-stranded"
    elif strandedness == "reverse":
        return "--rf-stranded"
    elif strandedness == "no":
        return ""
    else:
        raise ValueError(f"Strandedness not specified for sample {sample_id}.")


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
    return [r1, r2]


rule kallisto_quant:
    """
    Pseudoalign reads and quantify transcripts using kallisto.
    """
    input:
        # Read-Files (R1 and R2)
        fastq=get_reads_fasta_for_sample,
        # Index
        index=rules.kallisto_index.output.index
    output:
        # Only one argument is allowed by the wrapper !!!
        # Output directory
        directory(os.path.join(kallisto_output_path, "quantification", "quant_results_{sample_id}"))
        # os.path.join(kallisto_output_path, "quant_results_{sample_id}", "abundance.h5")
    params:
        # --plaintext: produces abundance.tsv -> A plain tsv-file with the estimated abundances for each transcript
        # --bias: Perform sequence based bias correction
        extra=lambda wildcards: "--bias " + get_kallisto_stranded_argument(wildcards.sample_id),
    log:
        "logs/kallisto_quant_{sample_id}.log",
    threads: 32
    resources:
        mem_mb=131072,          # total: assign 4096MB per CPU = 131072 MB, 131GB
        cpus=32,			    # uses 32 cpus
        time_min=360		    # use 6 hrs, to make sure
    wrapper:
        "v1.19.2/bio/kallisto/quant"


rule kallisto_create_output_annotation_table:
    """
        Creates an annotation table that is required for the analysis of the Salmon output
        via the R-script
    """
    input:
        quant_dirs=expand(os.path.join(kallisto_output_path, "quantification", "quant_results_{sample_id}"),
                          sample_id=pep.sample_table["sample_name"])
    output:
        # Annotation table -> copying PEP-annotation table
        annotation_table=os.path.join(kallisto_output_path, "{condition}",
                                      "kallisto_output_annotation_table.tsv")
    params:
        control_samples = pep.sample_table[(pep.sample_table["control"] == "true")]["sample_name"],
        condition_samples = main_helper_extract_sample_ids_with_specific_condition
    threads:
        1
    run:
        all_controls_table = pep.sample_table[pep.sample_table["sample_name"].isin(params.control_samples)]
        all_condition_table = pep.sample_table[pep.sample_table["sample_name"].isin(params.condition_samples)]
        # Create a table with all samples
        all_samples_table =all_controls_table.append(all_condition_table)

        # Add the quant.sf files to the table
        results_dir_path = os.path.dirname(input.quant_dirs[0])
        abundance_file_path = os.path.join(results_dir_path, "quant_results_{sample_id}", "abundance.h5")
        all_samples_table["kallisto_results_file"] = \
            all_samples_table.apply(lambda row: abundance_file_path.replace("{sample_id}", row["sample_name"]), axis=1)

        # Save output table
        all_samples_table.to_csv(output.annotation_table, sep="\t", index=False)