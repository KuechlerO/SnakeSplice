# ========== Submodule: Quality trimming of reads ==========
import os

trim_output_path = os.path.join('output', config["output_directories"]["trimmomatic_output_dir"])


def get_adapter_file_path(wildcards):
    adapter_file_path = pep.sample_table[pep.sample_table["sample_name"] == wildcards.sample_id]["adaptors_file"][0]
    return adapter_file_path


# 1. Perform trimming on fastq-files
rule run_trimmomatic_pe:
    """
    Run trimmomatic on paired-end fastq-files
    """
    input:
        # returns "r1" and "r2" as variables
        unpack(get_paired_read_file_paths),     # Defined in rules/lib_get_read_file_paths.smk, imported in Snakefile
        adaptors_file=get_adapter_file_path
    output:
        r1=os.path.join(trim_output_path, "{sample_id}_R1.fastq.gz"),
        r2=os.path.join(trim_output_path, "{sample_id}_R2.fastq.gz"),

        # reads where trimming entirely removed the mate
        r1_unpaired=os.path.join(trim_output_path, "{sample_id}_R1.unpaired.fastq.gz"),
        r2_unpaired=os.path.join(trim_output_path, "{sample_id}_R2.unpaired.fastq.gz")
    log:
        "logs/trimmomatic/{sample_id}.log"
    params:
        # list of trimmers (see manual)
        # 1. ILLUMINACLIP: Adaptor-file
        # -> 2 seed mismatches are allowed
        # -> 30 palindromeClipThreshold: specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment.
        # -> 7 simpleClipThreshold: specifies how accurate the match between any adapter etc. sequence must be against a read.
        #      each matching base adds just over 0.6, while each mismatch reduces the alignment score by Q/10. Therefore, a perfect match of a 12 base sequence will score just over 7
        # -> 2 minAdapterLength: minimum length of adapter
        # -> TRUE keepBothReads: if true, both reads are written to the output files
        # 2. TRAILING: drop trailing low quality bases
        # 3. LEADING: drop leading low quality bases
        # 4. SLIDINGWINDOW: Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15
        # 5. MINLEN: Drop reads which are less than 36 bases long
        trimmer=lambda wildcards, input:
        [
            "ILLUMINACLIP:" +
                input.adaptors_file +
                ":2:30:6:2:TRUE",
            "TRAILING:3",
            "LEADING:3",
            "SLIDINGWINDOW:4:15",
            "MINLEN:36"
        ],
        # optional parameters
        extra="",                                               # "-trimlog {sample_id}.extra.log",
        compression_level="-9"
    resources:
        mem_mb=1024,		# total memory usage: 1GB
        cpus=32				# uses 32 cpus
    threads:
        32
    wrapper:
        "v1.12.2/bio/trimmomatic/pe"