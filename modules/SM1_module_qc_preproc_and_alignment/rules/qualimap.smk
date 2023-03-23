# ========== Submodule: QualiMap - Assess quality of alignment ==========
import os

# Output directory
qualimap_output_path = os.path.join('output',config["output_directories"]["qualimap_output_dir"])

# A: Assess quality of alignment for Olego alignment
rule qualimap_statistics_for_olego_output:
    """
    A: Qualimap Statistics for Olego Alignment Results
    """
    input:
        rules.olego_samtools_sort.output		    # take sorted BAM-files as input
    output:
        output_dir=directory(os.path.join(qualimap_output_path, "olego/{sample_id}")),
        output_report=os.path.join(qualimap_output_path, "olego/{sample_id}/qualimapReport.html")
    params:
        java_heap_size="131G"
    log:
        "logs/qualimap_statistics/olego/{sample_id}.log"
    resources:
        mem_mb=131072,		# total: assign 4096MB per CPU = 131072
        cpus=32,			# uses 32 cpus
        time_min=180		# 3 hrs
    threads:
        32
    conda:
        "../envs/qualimap_env.yaml"
    shell:
        # nt: number of Threads
        "qualimap bamqc -nt {threads} "
        "--java-mem-size={params.java_heap_size} "
        "-bam {input} "
        "-outdir {output.output_dir} "
        "2>{log}"


# B: Qualimap Statistics for STAR Alignment Results
rule qualimap_statistics_for_star_output:
    """
    B: Qualimap Statistics for STAR Alignment Results
    """
    input:
        rules.star_samtools_sort.output		    # take sorted BAM-files as input
    output:
        output_dir=directory(os.path.join(qualimap_output_path, "star/{sample_id}")),
        output_report=os.path.join(qualimap_output_path, "star/{sample_id}/qualimapReport.html")
    params:
        java_heap_size=config["qualimap_settings"]["qualimap_java_heap_size"]
    log:
        "logs/qualimap_statistics/star/{sample_id}.log"
    resources:
        mem_mb=131072,		# total: assign 4096MB per CPU = 131072
        cpus=32,			# uses 32 cpus
        time_min=180		# 3 hrs
    threads:
        32
    conda:
        "../envs/qualimap_env.yaml"
    shell:
        # nt: number of Threads
        "qualimap bamqc -nt {threads} "
        "--java-mem-size={params.java_heap_size} "
        "-bam {input} "
        "-outdir {output.output_dir} "
        "2>{log}"
