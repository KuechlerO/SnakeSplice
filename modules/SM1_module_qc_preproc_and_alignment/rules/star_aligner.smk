# ========== Submodule: STAR - Alignment tool ==========
import os

# Output directory
star_output_path = os.path.join('output',config["output_directories"]["star_output_dir"])
# trimmomatic output directory
trimmomatic_output_path = os.path.join('output',config["output_directories"]["trimmomatic_output_dir"])

# 1. Build BWT-index for reference genome
rule star_build_index:
    input:
        fasta=config["alignment_settings"]["star_alignment_settings"]["reference_genome_annotation_file"],    # required by wrapper
        gtf=config["alignment_settings"]["star_alignment_settings"]["reference_genome_fasta_file"]	        # required by wrapper
    output:
        directory(os.path.join(star_output_path, "ref_genome_index"))
    resources:
        mem_mb=131072,		# total: assign 4096MB per CPU = 131072 MB, 131GB
        cpus=32,			# uses 32 cpus
        time_min=120		# use 2 hrs (with this settings normally only 30 mins are needed)
    threads:
        32
    params:
        sjdbOverhang=config["alignment_settings"]["star_alignment_settings"]["index_build"]["sjdbOverhang"],
        extra=config["alignment_settings"]["star_alignment_settings"]["index_build"]["extra_settings"]
    log:
        "logs/star_index_genome.log"
    wrapper:
        "0.80.2/bio/star/index"


def get_star_alignment_options():
    """
    Extracts alignment options from the YAML-module1-configuration file.
    :return:
    """
    extra = config["alignment_settings"]["star_alignment_settings"]["alignment"]["extra_settings"]
    if config["alignment_settings"]["star_alignment_settings"]["alignment"]["use_sjdb_file"]:
        extra += " --sjdbGTFfile " + config["alignment_settings"]["star_alignment_settings"]["reference_genome_annotation_file"]

    return extra


# 2. Aligning reads
# No specific STAR options are necessary for stranded Data
rule star_align:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1=os.path.join(trimmomatic_output_path, "{sample_id}_R1.fastq.gz"),
        fq2=os.path.join(trimmomatic_output_path, "{sample_id}_R2.fastq.gz"), #optional
        idx=rules.star_build_index.output
    output:
        # see STAR manual for additional output files
        aln=os.path.join(star_output_path, "{sample_id}/Aligned.out.bam"),		# output only BAM, not SAM -> Use samtools view to convert to SAM
        log=os.path.join(star_output_path, "{sample_id}/extra_log.out"),
        sj=os.path.join(star_output_path, "{sample_id}/ReadsPerGene.out.tab")
    log:
        "logs/star/{sample_id}_alignment.log"
    params:
        # path to STAR reference genome index
        index=rules.star_build_index.output,
        # optional parameters
        extra=get_star_alignment_options()
    resources:
        mem_mb=32*4096,		# total: assign 4096MB per CPU = 131072 MB, 131GB
        cpus=32,			# uses 32 cpus
        time_min=120		# use 2 hrs (with this setting normally only 30 mins are needed)
    threads: 32
    wrapper:
        "v1.21.0/bio/star/align"


# ----------------------- Converting output -------------------------
# 3. Rename output-BAM-files
rule star_rename_output_bam_files:
    input:
        rules.star_align.output.aln
    output:
        # Mark as temporary file in order to keep only sorted & indexed BAM files as final output
        # Otherwise: rm */*[!.sorted].bam
        temp(os.path.join(star_output_path, "{sample_id}/{sample_id}.bam"))
    shell:
        "mv {input} {output}"


# 4. SAMtools sort -> Sort by position
# Sort bam-files -> Also might give a better compression ratio (because similar sequences are grouped together)
# Resources: 32 Threads & 131 GB
rule star_samtools_sort:
    input:
        rules.star_rename_output_bam_files.output[0]
    output:
        os.path.join(star_output_path, "{sample_id}/{sample_id}.sorted.bam")
    log:
        "logs/samtools_sort/star_" + "{sample_id}.log"
    params:
        extra="-m 4G",			# max mem per thread
        tmp_dir="/tmp/"			# should be automatically included in shadow-prefix
    resources:
        mem_mb=32*4096,		# total: assign 4096MB per CPU = 131072
        cpus=25				# uses only 25 CPUS to leave space for overhead
    threads:  # Samtools takes additional threads through its option -@
        25	 # This value will be sent to -@.
    wrapper:
        "v1.21.0/bio/samtools/sort"


# 5. Index BAM-files (e.g. for Leafcutter)
# Resources: 32 Threads & 131 GB
rule star_samtools_index_bam_files:
    input:
        rules.star_samtools_sort.output[0]
    output:
        os.path.join(star_output_path, "{sample_id}/{sample_id}.sorted.bam.bai")
    log:
        "logs/samtools_index/{sample_id}.log"
    params:
        extra=""		    # Optional params string
    resources:
        mem_mb=131072,		# total: assign 4096MB per CPU = 131072
        cpus=32				# uses 32 cpus
    threads:  # Samtools takes additional threads through its option -@
        32	 # This value-1 will be sent to -@
    wrapper:
        "v1.14.0/bio/samtools/index"
