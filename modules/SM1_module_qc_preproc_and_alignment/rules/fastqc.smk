# ========== Submodule: FastQC - Creates Quality Control Overview DashBoard ==========
import os 

# Output directory
fastqc_output_path = os.path.join('output',config["output_directories"]["fastqc_output_dir"])

# ------------------- 1.) FastQC before trimming --------------------
rule fastqc_before_trimmomatic:
    """
    FastQC before trimming
    """
    input:
        get_single_read_file_paths      # Defined in rules/lib_get_read_file_paths.smk, imported in Snakefile
    output:
        html=os.path.join(fastqc_output_path, "non_trimmed/{sample_id}.html"),
        zip=os.path.join(fastqc_output_path, "non_trimmed/{sample_id}_fastqc.zip")  # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    threads:
        1			# Only one for one input
    log:
        "logs/fastqc_before_trimming/{sample_id}.log"
    wrapper:
        "0.79.0/bio/fastqc"


# ------------------- 2.) FastQC after trimming --------------------
rule fastqc_after_trimmomatic:
    """
    FastQC after trimming
    """
    input:
        trimmed_reads_file="output/trimmomatic/{sample_id}.fastq.gz",
    output:
        html=os.path.join(fastqc_output_path, "trimmed/{sample_id}.html"),
        zip=os.path.join(fastqc_output_path, "trimmed/{sample_id}_fastqc.zip")     # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    threads:
        1			# Only one for one input
    log:
        "logs/fastqc_after_trimming/{sample_id}.log"
    wrapper:
        "0.79.0/bio/fastqc"
