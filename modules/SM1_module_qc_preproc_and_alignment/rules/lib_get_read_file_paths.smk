import os


def get_single_read_file_paths(wildcards):
    """"
    Get the paths to the single read files.
    :params wildcards:  wildcards from the input file
                        Expected: {"sample": "76660_R1.fastq.gz"}
    :return: path to single read file
    """
    given_sample = wildcards.sample_id                             # e.g. {sample_id}_R1.fq
    read_extension = given_sample.split("_")[-1].split(".")[0]      # result: e.g. R1
    # Attention: Do also consider names like: sample_1_2_R1.fq.gz! [multiple underscores]
    sample_name = given_sample[:-len(given_sample.split("_")[-1])][:-1]  # results: e.g. {sample} / sample_1_2
    sample_dir = pep.sample_table["sample_directory"][sample_name]
    if read_extension == "R1":
        read1_file = pep.sample_table["read1"][sample_name]
        output_file = os.path.join(sample_dir,read1_file)
    elif read_extension == "R2":
        read2_file = pep.sample_table["read2"][sample_name]
        output_file = os.path.join(sample_dir,read2_file)
    else:
        raise Exception("Read extension not recognized")

    return output_file


def get_paired_read_file_paths(wildcards):
    """
    Returns both read file paths for the given sample.
    :param wildcards:   Wildcards from the rule.
                        Expected: {"sample": "76660"}
    :return:            Dictionary of read file paths.
    """
    sample_dir = pep.sample_table["sample_directory"][wildcards.sample_id]
    sample_read1 = pep.sample_table["read1"][wildcards.sample_id]
    sample_read2 = pep.sample_table["read2"][wildcards.sample_id]

    read1_file = os.path.join(sample_dir, sample_read1)
    read2_file = os.path.join(sample_dir, sample_read2)

    return {"r1": read1_file, "r2": read2_file}
