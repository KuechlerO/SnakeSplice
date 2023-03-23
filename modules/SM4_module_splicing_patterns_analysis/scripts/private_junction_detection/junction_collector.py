import pandas as pd


def merge_junctions_naivly(input_file_list, input_sample_list, output_file):
    """
    Merges the junctions from the different samples into one file.
    Naivly merging: If a junction is present in one sample, it is present in the merged file.

    :param input_file_list:     List of files containing the junctions from all samples
    :param input_sample_list:   List of sample names
    :param output_file:         Output file for merged junctions
    :return:
    """

    # Column-names: Each line is a exon-exon junction
    common_cols = ["chrom", "exact_jct_start", "exact_jct_end", "strand", "add_info"]       # Shared column names
    column_types = {"chrom": "category", "exact_jct_start": "uint32", "exact_jct_end": "uint32", "strand": "category"}
    summary_df = pd.DataFrame(columns=common_cols)
    summary_df.astype(column_types)

    # Iterate over all samples and merge the junctions into one file
    for counter, input_file in enumerate(input_file_list):
        sample_name = input_sample_list[counter]

        current_df = pd.read_csv(input_file, low_memory=False, sep="\t")
        current_df.fillna(0, inplace=True)
        reduced_df = current_df.iloc[:, 0:6]
        reduced_df.columns = ["chrom", "exact_jct_start", "exact_jct_end", "sample_name", "score", "strand"]
        reduced_df["add_info"] = current_df.iloc[:, -2]
        reduced_df.rename(columns={"score": sample_name}, inplace=True)     # Rename the score column to the sample name
        reduced_df = reduced_df[["chrom", "exact_jct_start", "exact_jct_end", "strand", "add_info", sample_name]]
        # Set datatypes to save memory
        reduced_df.astype(column_types)
        reduced_df.astype({sample_name: "uint16"})

        # Merge the junctions
        summary_df = summary_df.merge(reduced_df, on=common_cols, how="outer")

    # Save the results to a file
    summary_df.fillna(0).to_csv(output_file, index=False, sep="\t")


if __name__ == "__main__":
    # Input files
    snakemake_input_files = snakemake.input.all_junc_files
    snakemake_input_samples = snakemake.params.sample_names
    # Output file
    snakemake_output_file = snakemake.output.output_file

    # Merge the junctions
    merge_junctions_naivly(snakemake_input_files, snakemake_input_samples, snakemake_output_file)

