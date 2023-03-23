import pandas as pd


def extract_actual_junctions_from_regtools_file(input_jct_file, output_jct_file, sample_id):
    """
    Extracts the actual junctions from the regtools file
    Meaning:
    ChromStart includes maximum overhang for the junction on the left side -> Add blockSizes[0] to get the actual start
    ChromEnd includes maximum overhang for the junction on the right side -> Subtract blockSizes[1] to get the actual end
    See docs here: https://regtools.readthedocs.io/en/latest/commands/junctions-extract/

    :param input_jct_file:      Input file: Contains all junctions for given sample
    :param output_jct_file:     Output file
    :params sample_id:          ID of current sample

    :return:
    """

    # Column-names: Each line is a exon-exon junction
    column_names = ["chrom", "chromStart", "chromEnd", "name", "score", "strand",
                    "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"]

    # Read the file
    junction_df = pd.read_csv(input_jct_file, names=column_names, sep="\t")

    # Extract the actual junctions
    junction_df["max_overhang_before_start"] = junction_df["blockSizes"].str.split(",").str[0].astype(int)
    junction_df["max_overhang_after_end"] = junction_df["blockSizes"].str.split(",").str[-1].astype(int)
    junction_df["exact_jct_start"] = junction_df["chromStart"].astype(int) + junction_df["max_overhang_before_start"]
    junction_df["exact_jct_end"] = junction_df["chromEnd"].astype(int) - junction_df["max_overhang_after_end"]

    # Save the results to a file, but without header...
    junction_df["sample_name"] = sample_id
    reduced_df = junction_df[["chrom", "exact_jct_start", "exact_jct_end", "sample_name", "score", "strand"]]
    reduced_df.to_csv(output_jct_file, index=False, header=False, sep="\t")


if __name__ == "__main__":
    # Input files
    snakemake_input_file = snakemake.input.regtools_junc_files
    # Output file
    snakemake_output_file = snakemake.output.output_file

    sample_id = snakemake.wildcards.sample_id

    # Extract the actual junctions
    extract_actual_junctions_from_regtools_file(snakemake_input_file, snakemake_output_file, sample_id)
