import os

# Output directory
salmon_output_path = os.path.join('output', config["output_directories"]["salmon_output_dir"])

rule salmon_prepare_files_to_be_indexed:
    """
    Prepares genome files for salmon indexing
    """
    input:
        transcriptome=config["salmon_settings"]["transcriptome_fasta"],
        genome=config["salmon_settings"]["genome_fasta"]
    output:
        gentrome=os.path.join(salmon_output_path, "reference_transcriptome", "gentrome.fasta.gz"),
        decoys=os.path.join(salmon_output_path, "reference_transcriptome", "decoys.txt"),
    threads: 32
    resources:
        mem_mb=131072,	        # total: assign 4096MB per CPU = 131072 MB, 131GB
        cpus=32,			    # uses 32 cpus
        time_min=360		    # use 6 hrs, to make sure
    log:
        "logs/salmon/decoys.log"
    wrapper:
        "v1.19.2/bio/salmon/decoys"


rule salmon_create_transcriptome_index:
    """
    Creates a salmon index for the transcriptome
    """
    input:
        sequences=rules.salmon_prepare_files_to_be_indexed.output.gentrome,
        decoys=rules.salmon_prepare_files_to_be_indexed.output.decoys
    output:
        multiext(
            os.path.join(salmon_output_path, "transcriptome_index/"),
            "complete_ref_lens.bin",
            "ctable.bin",
            "ctg_offsets.bin",
            "duplicate_clusters.tsv",
            "info.json",
            "mphf.bin",
            "pos.bin",
            "pre_indexing.log",
            "rank.bin",
            "refAccumLengths.bin",
            "ref_indexing.log",
            "reflengths.bin",
            "refseq.bin",
            "seq.bin",
            "versionInfo.json",
        ),
    log:
        "logs/salmon/transcriptome_index.log",
    threads: 32                 # needs at least 2 threads to run
    resources:
        mem_mb=131072,          # total: assign 4096MB per CPU = 131072 MB, 131GB
        cpus=32,		        # uses 32 cpus
        time_min=360		    # use 6 hrs, to make sure
    params:
        # optional parameters
        extra="",
    wrapper:
        "v1.19.2/bio/salmon/index"


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
    return {"r1": r1, "r2": r2}


# Integrate replicate consideration?
# -> Or just include replicate considerations in the downstream analysis?
rule salmon_quant_reads:
    input:
        # If you have multiple fastq files for a single sample (e.g. technical replicates)
        # use a list for r1 and r2.
        # Read-Files -> Function defined in Main-Snakemake/scripts/helper_functions.smk
        unpack(get_reads_fasta_for_sample),
        # index created by salmon_create_transcriptome_index
        index=rules.salmon_create_transcriptome_index.output
    output:
        quant=os.path.join(salmon_output_path, "quantification", "quant_results_{sample_id}", "quant.sf"),
        lib=os.path.join(salmon_output_path, "quantification", "quant_results_{sample_id}", "lib_format_counts.json"),
    log:
        "logs/salmon/{sample_id}.log",
    params:
        # optional parameters
        libtype="A",        # Allow Salmon to infer the library type automatically
        # Extras
        # --gcBias: Perform fragment GC bias correction: correct for biases in how likely a sequence is to be observed based on its internal GC content.
        # --numGibbsSamples:    Allows to estimate the variance in abundance estimates. However, in this case the samples are generated using posterior Gibbs sampling over the fragment equivalence classes rather than bootstrapping.
        extra="--gcBias --numGibbsSamples 20",
    threads: 32
    resources:
        mem_mb=131072,      # total: assign 4096MB per CPU = 131072 MB, 131GB
        cpus=32,	        # uses 32 cpus
        time_min=360		# use 6 hrs, to make sure
    wrapper:
        "v1.19.2/bio/salmon/quant"


rule salmon_create_output_annotation_table:
    """
    Creates an annotation table that is required for the analysis of the Salmon output
    via the Deseq2
    """
    input:
        quant_files=expand(os.path.join(salmon_output_path, "quantification", "quant_results_{sample_id}", "quant.sf"),
                            sample_id=pep.sample_table["sample_name"])
    output:
        # Annotation table -> copying PEP-annotation table
        annotation_table=os.path.join(salmon_output_path, "{condition}", "salmon_output_annotation_table.tsv")
    params:
        control_samples=pep.sample_table[(pep.sample_table["control"] == "true")]["sample_name"],
        condition_samples=main_helper_extract_sample_ids_with_specific_condition
    threads:
        1
    run:
        all_controls_table = pep.sample_table[pep.sample_table["sample_name"].isin(params.control_samples)]
        all_condition_table = pep.sample_table[pep.sample_table["sample_name"].isin(params.condition_samples)]
        # Create a table with all samples
        all_samples_table =all_controls_table.append(all_condition_table)

        # Add the quant.sf files to the table
        results_dir_path = os.path.dirname(os.path.dirname(input.quant_files[0]))
        quant_file_path = os.path.join(results_dir_path, "quant_results_{sample_id}", "quant.sf")
        all_samples_table["salmon_results_file"] = \
            all_samples_table.apply(lambda row: quant_file_path.replace("{sample_id}", row["sample_name"]), axis=1)

        # Save output table
        all_samples_table.to_csv(output.annotation_table, sep="\t", index=False)
