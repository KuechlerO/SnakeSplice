import pandas as pd


if __name__ == '__main__':

    count_files = snakemake.input["count_files"]
    output_total_count_file = snakemake.output["total_counts_file"]

    # Iterate over all count files and merge them into one dataframe
    total_counts_df = pd.DataFrame()
    for i, count_file in enumerate(count_files):
        # Current sample ID
        current_sample_id = count_file.split("/")[-1].split(".")[0]
        # Read count file & add sample ID as column
        current_count_df = pd.read_csv(count_file, sep="\t", low_memory=False, skiprows=1)
        current_count_df = current_count_df.set_axis([*current_count_df.columns[:-1], current_sample_id], axis=1)

        if i == 0:  # First count file
            total_counts_df = current_count_df
        else:    # All other count files
            # Add only last column to the total count dataframe
            total_counts_df = total_counts_df.merge(current_count_df,  how="outer",
                                                    on=["Geneid", "Chr", "Start", "End", "Strand", "Length"])

    # Write total count dataframe to file
    total_counts_df.to_csv(output_total_count_file, sep="\t", index=False)
