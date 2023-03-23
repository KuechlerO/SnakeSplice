import os

# Output directory for all generated files of this Tool
irfinder_output_path = os.path.join('output', config["output_directories"]["irfinder_output_dir"])

rule irfinder_sort_bam_by_name:
    """
    This rule sorts the bam file by name.
    Reason: IRFinder needs input bam files to be sorted by read name
    """
    input:
        # Returns path with {sample_id} as wildcard
        main_helper_get_bam_file_for_sample(config["input_dir_of_bam_files"], None, config["bam_files_attributes"]["filename_extension"])
    output:
        os.path.join(irfinder_output_path, "sorted_bams/{sample_id}.sorted.bam")
    log:
        "logs/irfinder/bam_name_sorting_{sample_id}.log",
    params:
        extra="-m 4G",
    threads: 8
    resources:
        cpus=8,
        mem_mb=8 * 4096,
        time_min=60
    wrapper:
        "v1.18.3/bio/samtools/sort"


rule irfinder_create_irfinder_reference_by_local_instance:
    """
    This rule creates the IRFinder reference (extracting the potential introns).
    """
    input:
        gtf_file=config["irfinder_settings"]["ensemble_sourced_gtf_file"],
        fasta_file=config["irfinder_settings"]["ensemble_sourced_fasta_file"]
    output:
        ref_dir=directory(os.path.join(irfinder_output_path, "irfinder_reference")),
        bed_file=os.path.join(irfinder_output_path, "irfinder_reference/introns.unique.bed")
    log:
        "logs/irfinder/irfinder_reference_creation.log"
    params:
        bed_file_consensus_excludable = "../data/irfinder/Human_hg19_wgEncodeDacMapabilityConsensusExcludable.bed.gz",
        bed_file_non_polya = "../data/irfinder/Human_hg19_nonPolyA_ROI.bed"
    threads:
        8
    resources:
        cpus=8,
        mem_mb=15 * 4096,   # Allocate 60 GB of memory
        time_min=120
    conda:
        "../envs/irfinder_env.yaml"
    shell:
        "mkdir -p {output.ref_dir};"
        "ln -s {input.gtf_file} {output.ref_dir}/transcripts.gtf;"
        "ln -s {input.fasta_file} {output.ref_dir}/genome.fa;"
        "IRFinder -m BuildRefProcess -t {threads} "
        "-r {output.ref_dir} "
        "-b {params.bed_file_consensus_excludable} "
        "-R {params.bed_file_non_polya} "
        "2> {log};"


# rule irfinder_create_irfinder_reference_by_download:
#     """
#     This rule creates the IRFinder reference (extracting the potential introns).
#     """
#     output:
#         ref_dir=directory("output/irfinder/irfinder_reference"),
#         bed_file="output/irfinder/irfinder_reference/introns.unique.bed"
#     log:
#         "logs/irfinder/irfinder_reference_creation.log"
#     params:
#         # ensembl_link=config["irfinder_settings"]["ensemble_source_link"],
#         ensembl_link="ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz"
#     threads:
#         8
#     resources:
#         cpus=8,
#         mem_mb=15 * 4096,   # Allocate 60 GB of memory
#         time_min=120
#     conda:
#         "../envs/irfinder_env.yaml"
#     shell:
#         "if [ -d \"{output.ref_dir}\" ]; then rm -Rf {output.ref_dir}; fi && "     # Remove the directory if it already exists
#         "IRFinder -m BuildRef -t {threads} "
#         "-r {output.ref_dir} "
#         "{params.ensembl_link} "
#         "2> {log}"


rule irfinder_quantify_irf:
    input:
        unique_introns_bed_file=rules.irfinder_create_irfinder_reference_by_local_instance.output.bed_file,
        bam_file=
            # Returns path with {sample_id} as wildcard
            rules.irfinder_sort_bam_by_name.output if config["bam_files_attributes"]["sorted"] != "name"
            else
                main_helper_get_bam_file_for_sample(config["input_dir_of_bam_files"], None,
                                        config["bam_files_attributes"]["filename_extension"]),
        reference_dir=rules.irfinder_create_irfinder_reference_by_local_instance.output.ref_dir
    output:
        output_dir=directory(os.path.join(irfinder_output_path, "{sample_id}")),
        non_directional_rna_seq_file=os.path.join(irfinder_output_path, "{sample_id}/IRFinder-IR-nondir.txt"),
        directional_rna_seq_file=os.path.join(irfinder_output_path, "{sample_id}/IRFinder-IR-dir.txt")
    log:
        "logs/irfinder/irf_quantification_{sample_id}.log"
    params:
        dummy="TODO: Remove"
        # This will cause a bug
        # reference_dir = directory("output/irfinder/irfinder_reference")
    conda:
        "../envs/irfinder_env.yaml"
    threads:
        8
    resources:
        cpus=8,
        mem_mb=8 * 4096,
        time_min=60
    shell:
        "IRFinder -m BAM "
        "-r {input.reference_dir} "
        "-d {output.output_dir} "
        "-t {threads} "
        "{input.bam_file} 2> {log}"


def get_irf_sample_result_file_paths(wildcards):
    """
    Return a list of file paths of the IRFinder results for the given sample.
    Either the non-directional or directional results can be returned.
    -> depends on settings in "irfinder_settings/strandedness_of_input" in config_module4_splicing_patterns_analysis.yaml

    param wildcards:    Wildcards object, holds the condition information
    """
    control = "None"
    stranded = config["irfinder_settings"]["strandedness_of_input"]

    # Extract needed sample names (according to condition)
    condition_and_control_samples = \
        (pep.sample_table[(pep.sample_table["condition"] == wildcards.condition)
                          | (pep.sample_table["condition"] == control)]["sample_name"]).tolist()

    # Get template urls
    if stranded == "no":
        sample_path_with_wildcard = os.path.join(irfinder_output_path, "{sample_id}/IRFinder-IR-nondir.txt")
    else:
        sample_path_with_wildcard = os.path.join(irfinder_output_path, "{sample_id}/IRFinder-IR-dir.txt")

    # Plug sample-IDs into file-paths
    all_paths = []
    for sample in condition_and_control_samples:
        all_paths.append(sample_path_with_wildcard.format(sample_id=sample))

    return all_paths


rule irfinder_collect_sample_result_file_paths_in_file:
    """
    Collects all IRF sample result file paths in a file.
    Each path is written in a new line.
    
    Attention: Needs to be paired with the output of rule "irfinder_create_sample_condition_mapping_file".
    """
    input:
        get_irf_sample_result_file_paths
    output:
        os.path.join(irfinder_output_path, "{condition}/irf_sample_result_file_paths.txt")
    threads:
        1
    run:
        with open(output[0], "w") as f:
            for path in input:
                f.write(path + "\n")


rule irfinder_create_sample_condition_mapping_file:
    """
    Creates a mapping file for the IRF sample result files:
    Each row has the following format: <sample_name> <condition>
    
    Attention: Needs to be paired (i.e. samples in same order) with the 
    output of rule "irfinder_collect_sample_result_file_paths_in_file".
    """
    output:
        os.path.join(irfinder_output_path, "{condition}/sample_condition_mapping.txt")
    params:
        samples = lambda wildcards:
                        pep.sample_table[pep.sample_table["condition"] == wildcards.condition]["sample_name"].tolist(),
        conditions = lambda wildcards:
                        pep.sample_table[pep.sample_table["condition"] == wildcards.condition]["condition"].tolist()
    threads:
        1
    run:
        with open(output[0], "w") as f:
            assert(len(params.samples) == len(params.conditions))
            f.write("SampleNames\tCondition\n")
            for i in range(len(params.samples)):
                f.write(params.samples[i] + "\t" + params.conditions[i] + "\n")


rule irfinder_dif_ir_analysis_with_glm:
    """
    Differential Intron-Retention analysis.
    
    This rule performs differential IR analysis with GLM.
    Attention: input.irfinder_results_file_paths_collection & input.sample_condition_mapping_file need
    to be paired (i.e. having same sample listing order) with the output of 
    rule "irfinder_collect_sample_result_file_paths_in_file".
    """
    input:
        irfinder_results_file_paths_collection=rules.irfinder_collect_sample_result_file_paths_in_file.output,
        sample_condition_mapping_file=rules.irfinder_create_sample_condition_mapping_file.output,
        dexseq2_constructor_script=workflow.source_path("../scripts/irfinder/DESeq2Constructor.R"),
        irfinder_analysis_script=workflow.source_path("../scripts/irfinder/run_glm_analysis.R")
    output:
        output_plot_file=os.path.join(irfinder_output_path, "{condition}/glm_analysis.jpeg"),
        output_results_csv_file=os.path.join(irfinder_output_path, "{condition}/glm_analysis_results.csv")
    threads:
        1
    conda:
        "../envs/irfinder_env.yaml"
    script:
        "../scripts/irfinder/run_glm_analysis.R"
