# ========== Submodule: Kraken2 - Check for potential contamination ==========
import os

# Output directory
kraken_output_path = os.path.join('output',config["output_directories"]["kraken2_output_dir"])
# Trimmomatic directory
trimmomatic_output_path = os.path.join('output',config["output_directories"]["trimmomatic_output_dir"])

# Kraken2 -> a taxonomic classification system using exact k-mer matches
# to achieve high accuracy and fast classification speeds.

minikraken2_download_link = "ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/old/minikraken2_v1_8GB_201904.tgz"

# Check whether kraken2-database has been specified in config-file. If not, then download it (with following rule).
if config["contamination_check_settings"]['file_kraken2_db']:
    kraken2_db = config["contamination_check_settings"]['file_kraken2_db']
else:
    kraken2_db = "data/" + minikraken2_download_link.split("/")[-1].split(".")[0]


rule download_minikraken2_v1_databse:
    """
    Download the miniKraken2 database.
    """
    output:
        db_dir=directory(kraken2_db),
        db_dir_content=os.path.join(kraken2_db, "taxo.k2d")
    params:
        minikraken2_v1_db_link=minikraken2_download_link,
        download_file=minikraken2_download_link.split("/")[-1].split(".")[0]
    shell:
        "wget {params.minikraken2_v1_db_link};"
        "tar -xvzf {params.download_file};"
        "mv {params.download_file} {output.db_dir};"


# Resources: 32 Threads & 131 GB
rule kraken2_mapping:
    """
    Run Kraken2 on the trimmed read files. 
    This wil perform a check for potential contamination.
    """
    input:
        # Kraken2-DB: If not specified in config, then download will be called before
        db_dir_content=os.path.join(kraken2_db, "taxo.k2d"),

        # take paired trimmed reads as input
        r1=os.path.join(trimmomatic_output_path, "{sample_id}_R1.fastq.gz"),
        r2=os.path.join(trimmomatic_output_path, "{sample_id}_R2.fastq.gz"),

        # reads where trimming entirely removed the mate
        r1_unpaired=os.path.join(trimmomatic_output_path, "{sample_id}_R1.unpaired.fastq.gz"),
        r2_unpaired=os.path.join(trimmomatic_output_path, "{sample_id}_R2.unpaired.fastq.gz")
    output:
        kraken2_kmer_mapping = os.path.join(kraken_output_path, "{sample_id}.report"),
        kraken2_report = os.path.join(kraken_output_path, "{sample_id}_report_summary.report")
    params:
        kraken2_db=kraken2_db
    log:
        "logs/kraken2/{sample_id}.log"
    resources:
        mem_mb=131072,		# total: assign 4096MB per CPU = 131072 = 130 GB
        cpus=32,			# uses 32 cpus
        time_min=120		# uses 2 hrs
    threads:
        32
    conda:
        "../envs/kraken2_env.yaml"
    shell:
        "kraken2 --use-names --threads {threads} --db {params.kraken2_db} "
        "--report {output.kraken2_report} "
        "--paired {input.r1} {input.r2} "
        "> {output.kraken2_kmer_mapping} 2> {log}"
