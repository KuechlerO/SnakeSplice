# This python script executes the ALFA analysis for multiple bam files in parallel.
# It was set up to make sure the file pathways are correctly imported (considering the module output structure).
# (otherwise problems can occur due to the module structure)

import os

if __name__ == '__main__':

    # input
    input_bam_files = snakemake.input.bam_files
    # output
    output_dir = snakemake.output[0]
    # log file
    log_file = snakemake.log[0]
    # params
    alfa_genome_index_name = snakemake.params.alfa_genome_index_name
    nr_processes = snakemake.threads

    # Creates a string of bam-files and their respective labels.
    #     Format: BAM_FILE1 LABEL1 [BAM_FILE2 LABEL2 …]
    bam_files_and_labels = " ".join(["{0} {1}".format(bam_file, bam_file.replace(".sorted.bam", "").split("/")[-1]) for
                                     bam_file in input_bam_files])

    command = "alfa -g {alfa_genome_index_name} "\
              "--bam {bam_files_with_labels} "\
              "-o {output_dir} "\
              "--processors {threads} "\
              "2> {log}".format(alfa_genome_index_name=alfa_genome_index_name,
                                bam_files_with_labels=bam_files_and_labels,
                                output_dir=output_dir,
                                threads=nr_processes,
                                log=log_file)

    os.system(command)  # execute command
