import os

# Output directory
cummerbund_output_path = os.path.join('output', config["output_directories"]["cummerbund_output_dir"])

# Original data directory -> Cufflinks output
cufflinks_output_path = os.path.join('output', config["output_directories"]["cufflinks_output_dir"])

rule cummerbund:
    """
    Use cummerbund to summarizes and visualize the data produced by cufflinks.
    """
    input:
        merged_cufflinks_transcriptome_assemblies_gtf=rules.cufflinks_merge_assemblies.output[0],
        # wildcards: condition
        analysis_output=rules.cufflinks_complete_analysis.output
    output:
        cummerbund_output_dir = directory(os.path.join(cummerbund_output_path, "{condition}")),
        cummerbund_summary_results = os.path.join(cummerbund_output_path, "{condition}", "analysis_output.txt")
    params:
        cufflinks_output_files_dir=lambda wildcards, input:
            os.path.dirname(input.analysis_output[0]),
        original_genome_build = config["cummerbund_settings"]["genome_build"],
        chosen_genes_of_interest = config["cummerbund_settings"]["genes_for_detailed_overview"]
    log:
        "logs/cummerbund_analysis_{condition}.log"
    threads:
        1
    resources:
        mem_mb=80000		# total: assign 80GB
    conda:
        "../envs/conda_cummerbund_env.yaml"
    script:
        "../scripts/cummerbund_script.R"