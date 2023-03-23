#
# Note:
# To make use of these functions and rules, you must include the
# following line (before any other includes!) in your Snakefile:
# include:    "./../../scripts/main_helper_functions.smk"
# -----------------------------------------------------------------------------

import os
from pathlib import Path
import warnings


def main_helper_get_bam_file_for_sample(input_dir_choice, sample_id=None, filename_extension='.sorted.bam'):
    """
    Identify input BAM file for a given sample_id.
    This includes a check, whether the file exists.
    BAM files can be either saved directly in the given input_dir_choice, or in a subdirectory of the input input_dir_choice
    E.g.: input_dir_choice = '/path/to/input_dir' and sample_id = 'sample_1' can result in:
    '/path/to/input_dir/sample_1.bam' or '/path/to/input_dir/sample_1/sample_1.bam'

    :param input_dir_choice:    Chosen input dir. Either star, olego, or real absolute path
    :param sample_id:           Sample ID, if None, then return wildcard: '{sample_id}'
    :param filename_extension:  File extension of the BAM file
    :return: Absolute URL leading to BAM-file of given sample_id
    """

    current_file_dir = Path(workflow.snakefile).parent.absolute()           # Snakemake_Main/scripts
    # current_rule_directory = Path(workflow.current_basedir).absolute()    # File of rule, which is calling this fct

    # 1. Identify directory, while still including wildcard for sample_id
    base_path = os.path.join(current_file_dir, "..", config["output_dir_module1_qc_preprocessing_and_alignment"])
    if input_dir_choice == "star":
        base_path = os.path.join(base_path, os.path.join("output", "star"))
    elif input_dir_choice == "olego":
        base_path = os.path.join(base_path, os.path.join("output", "olego"))
    else:
        base_path = input_dir_choice

    # 2. Finally, add suffix leading to BAM-file for sample
    wildcard_file_path = os.path.join(base_path,"{sample_id}/{sample_id}" + filename_extension)
    if sample_id is not None:
        # Return absolute path to BAM-file -> replace wildcard with sample_id and check, whether file exists
        if os.path.isfile(wildcard_file_path.replace("{sample_id}", sample_id)):
            final_file_path = wildcard_file_path.replace("{sample_id}", sample_id)
        elif os.path.isfile(wildcard_file_path.replace("{sample_id}/{sample_id}", sample_id)):
            final_file_path = wildcard_file_path.replace("{sample_id}/{sample_id}", sample_id)
        else:
            # If external directory for BAM-inputs is given, then BAM-files have to exist
            if input_dir_choice not in ["star", "olego"]:
                raise RuntimeError("Could not find BAM file for sample_id: " + sample_id + ". Path: " + wildcard_file_path)
            else:
                # For first module, where BAM-files are not yet created...
                warnings.warn("Could not find BAM file for sample_id: " + sample_id + ". Path: " + wildcard_file_path)
                final_file_path = wildcard_file_path.replace("{sample_id}", sample_id)
                warnings.warn("Return standard path for BAM file: " + final_file_path)
    else:
        final_file_path = wildcard_file_path    # Return wildcards for sample_id

    return final_file_path


def main_helper_get_all_bam_file_paths(sample_id_array, bam_files_dir, filename_extension='.sorted.bam'):
    """
    Collects paths for all BAM-files of all Samples
    :param sample_id_array:         Array of sample IDs
    :param bam_files_dir:           Directory, where BAM-files are stored
    :param filename_extension:     File extension of the BAM-files
    :return:
    """
    # sample_id_array = pep.sample_table["sample_name"]
    # bam_files_dir = config["input_dir_of_bam_files"]
    all_bam_file_paths = []
    for sample_id in sample_id_array:
        # references help-script in "Snakemake_Main/scripts"
        all_bam_file_paths.append(main_helper_get_bam_file_for_sample(bam_files_dir, sample_id=sample_id,
            filename_extension=filename_extension))

    return all_bam_file_paths


def main_helper_return_required_gtf_file(require_conversion_bool, original_gtf_file_path, output_dir, wildcards=None):
    """
    Return the required gtf file path, depending on the whether a chromosome name transformation is needed or not.
    ATTENTION: This takes subsequently advantage of function helper_rule_transform_chromosome_names_in_gtf in
    the main snakemake directory's scripts/helper_functions.py.
    :param require_conversion_bool:     True, False -> Whether name conversion shall be done or not
    :param original_gtf_file_path:      Original GTF-file of interest
    :param output_dir:                  Target output directory in case of needed transformation
    :param wildcards:
    :return:                            Path to (transformed) GTF-file
    """
    if require_conversion_bool:
        transformed_gtf_file_name = os.path.basename(
            original_gtf_file_path.replace(".gtf", ".transformed_chro_names.gtf"))
        gtf_file = os.path.join(output_dir, transformed_gtf_file_name)
    else:
        gtf_file = original_gtf_file_path
    return gtf_file


def main_helper_extract_sample_ids_with_specific_condition(wildcards):
    """
    Depending on the condition-wildcard, this function returns a list of sample ids
    :param wildcards:
    :return:
    """
    condition = wildcards.condition
    return pep.sample_table[pep.sample_table["condition"] == condition]["sample_name"]


# TODO can be deleted? -> Or replaced by function...
rule main_helper_export_pep_table_with_annotations:
    output:
        temp("output/temporary_output/pep_table_with_annotations.csv")
    run:
        pep.sample_table.to_csv(output[0], index=False)

