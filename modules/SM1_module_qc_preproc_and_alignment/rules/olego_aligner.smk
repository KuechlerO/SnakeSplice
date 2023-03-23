# ========== Submodule: Olego - Alignment tool ==========
import os

# Output directory
olego_output_path = os.path.join('output',config["output_directories"]["olego_output_dir"])
# trimmomatric output directory
trimmomatic_output_path = os.path.join('output',config["output_directories"]["trimmomatic_output_dir"])


# Catch whether to install Olego
# Capture olego installation directory
olego_installation_dir=os.path.join(config["alignment_settings"]["olego_alignment_settings"]["olego_installation_dir"],
    "olego")

# Load regression_model
olego_regression_model = config["alignment_settings"]["olego_alignment_settings"]["olego_regression_model"]   # regression model
if olego_regression_model is None or olego_regression_model == "":
    olego_regression_model = os.path.join(olego_installation_dir, "models/hg.cfg")


# 0. Install Olego:
rule olego_install_olego:
    output:
        os.path.join(config["alignment_settings"]["olego_alignment_settings"]["olego_installation_dir"], "olego")
    params:
        olego_dir = config["alignment_settings"]["olego_alignment_settings"]["olego_installation_dir"],
        olego_url= "https://github.com/chaolinzhanglab/olego.git"
    shell:
        "mkdir -p {params.olego_dir};"
        "cd {params.olego_dir};"
        "git clone {params.olego_url};"
        "cd olego;"
        "make"


# 1. Build the BWT index in order to run OLego successfully
# ATTENTION: A lot of computing power is needed
# Just take 32 CPUs & 60 GB
rule olego_create_bwt_index_for_reference_genome:
    input:
        ref_genome_fasta=config["alignment_settings"]["olego_alignment_settings"]["reference_genome_fasta_file"],
        olego_installed=rules.olego_install_olego.output
    output:
        config["alignment_settings"]["olego_alignment_settings"]["reference_genome_fasta_file"] + ".bwt"		# Does not collide with STAR, since STAR creates own Index-DIR
    resources:
        mem_mb=60000,		# total: assign 2006MB per CPU = 60GB
        cpus=32,			# uses 32 cpus
        time_min=120		# use 2 hrs
    params:
        olego_installation_dir=olego_installation_dir
    threads:
        32
    log:
        "logs/olego/olego_create_index_for_reference_genome.log"
    shell:
        "{params.olego_installation_dir}/olegoindex {input} 2> {log}"


# 2. Alignment
# -> Aligning firstly without considering pairing mode
# Use 32 CPUs & 131 GB & 3 hrs
rule olego_rna_alignment:
    input:
        os.path.join(trimmomatic_output_path, "{sample_id}_R{read_file_nr}.fastq.gz"),							# take trimmed reads
        rules.olego_create_bwt_index_for_reference_genome.output,	    # BWT index of reference genome
        olego_installed=rules.olego_install_olego.output
    output:
        os.path.join(olego_output_path, "{sample_id}/{sample_id}_R{read_file_nr}.sam")						            # alignment results
    params:
        olego_installation_dir=olego_installation_dir,
        # Also possible: -j: Annotation file for known exon junctions
        r=config["alignment_settings"]["olego_alignment_settings"]["olego_regression_model"],   # regression model
        M=config["alignment_settings"]["olego_alignment_settings"]["olego_allowed_missmatches"],    # default: 4
        e=config["alignment_settings"]["olego_alignment_settings"]["olego_min_exon_size"],	# Minimum micro-exon size to be searched (INT) [ default: 9 ].
        ref_index=config["alignment_settings"]["olego_alignment_settings"]["reference_genome_fasta_file"]			# Indexed reference genome
    resources:
        mem_mb=131072,		# total: assign 4096MB per CPU = 131072
        cpus=32,				# uses 32 cpus
        time_min=180		# use 3 hrs
    threads:
        32
    log:
        "logs/olego/alignment/{sample_id}_R{read_file_nr}.log"
    shell:
        "{params.olego_installation_dir}/olego -v -t {threads} -r {params.r} -M {params.M} -o {output}  "
        "{params.ref_index} {input[0]} 2> {log}"


# 3. Olego: Merging both read-ends into merge.sam
# Uses 30 min on 1 CPU
# Try allocating 8 CPUs, 8GB, 2hrs -> Does it work?
rule olego_merge_both_ends:
    input:
        os.path.join(olego_output_path, "{sample_id}/{sample_id}_R1.sam"),
        os.path.join(olego_output_path, "{sample_id}/{sample_id}_R2.sam"),
        olego_installed=rules.olego_install_olego.output
    output:
        os.path.join(olego_output_path, "{sample_id}/{sample_id}_merged.sam")
    log:
        "logs/olego/merging/{sample_id}.log"
    resources:
        mem_mb=8000,		# total: assign 4096MB per CPU = 8GB
        cpus=8,				# uses 8 cpus
        time_min=120		# uses 2 hrs
    params:
        olego_installation_dir=olego_installation_dir
    threads:
        8
    shell:
        "perl {params.olego_installation_dir}/mergePEsam.pl -v {input[0]} {input[1]} {output} 2> {log}"



# ----------------------------------- Converting output ----------------------------------------
# 4. Convert SAM-files to BAM-files
# Takes 15 min on 1 CPU
# Default resources: 1 CPU, 2GB
rule olego_samtools_convert_sams_to_bams:
    input:
        rules.olego_merge_both_ends.output,
    output:
        os.path.join(olego_output_path, "{sample_id}/{sample_id}_merged.bam")
    log:
        "logs/olego/convert_sam_to_bam_{sample_id}.log"
    params:
        extra="" 			# optional params string
    wrapper:
        "v1.14.0/bio/samtools/view"


# 5. SAMtools sort -> Sort by position
# Sort bam-files -> Also might give a better compression ratio (because similar sequences are grouped together)
# Resources: 32 Threads & 131 GB
rule olego_samtools_sort:
    input:
        rules.olego_samtools_convert_sams_to_bams.output
    output:
        os.path.join(olego_output_path, "{sample_id}/{sample_id}.sorted.bam")
    log:
        "logs/samtools_sort/olego_" + "{sample_id}.log"
    params:
        extra="-m 4G",			# max mem per thread
        tmp_dir="/tmp/"			# should be automatically included in shadow-prefix
    resources:
        mem_mb=131072,		# total: assign 4096MB per CPU = 131072
        cpus=32				# uses 32 cpus
    threads:  # Samtools takes additional threads through its option -@
        32	 # This value will be sent to -@.
    wrapper:
        "v1.14.0/bio/samtools/sort"


# 6. Index BAM-files (e.g. for Leafcutter)
# Resources: 32 Threads & 131 GB
rule olego_samtools_index_bam_files:
    input:
        rules.olego_samtools_sort.output
    output:
        os.path.join(olego_output_path, "{sample_id}/{sample_id}.sorted.bam.bai")
    log:
        "logs/samtools_index/{sample_id}.log"
    params:
        extra=""  # Optional params string
    resources:
        mem_mb=131072,	# total: assign 4096MB per CPU = 131072
        cpus=32				# uses 32 cpus
    threads:  # Samtools takes additional threads through its option -@
        32	 # This value-1 will be sent to -@
    wrapper:
        "v1.14.0/bio/samtools/index"


#  --------------------------- Not needed ------------------------------------
# 7. Olego convert SAM file to BED file
rule olego_convert_sam_to_bed:
    input:
        rules.olego_merge_both_ends.output,
        olego_installed=rules.olego_install_olego.output
    output:
        os.path.join(olego_output_path, "{sample_id}_merged.bed")
    log:
        "logs/olego/convert_sam_to_bed/{sample_id}.log"
    resources:
        mem_mb=8000,		# total: assign 4096MB per CPU = 8GB
        cpus=8,				# uses 32 cpus
        time_min=120		# uses 2 hrs
    params:
        olego_installation_dir=olego_installation_dir
    threads:
        8
    shell:
        "perl {params.olego_installation_dir}/sam2bed.pl -v --use-RNA-strand {input} {output} 2> {log}"


# 8. Find the junctions in the BED file (transforms BED to JUNC)
# Uses 2min on 1 CPU
# Use default resources
rule olego_find_junctions:
    input:
        rules.olego_convert_sam_to_bed.output,
        olego_installed=rules.olego_install_olego.output
    output:
        os.path.join(olego_output_path, "{sample_id}_merged.junc")
    log:
        "logs/olego/find_junctions/{sample_id}.log"
    params:
        olego_installation_dir=olego_installation_dir
    shell:
        "perl {params.olego_installation_dir}/bed2junc.pl {input} {output} 2> {log}"
