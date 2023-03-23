import pandas as pd


def clean_summary_csv_file(input, output):
    """
    Replaces the vector symbols in the last column of the summary file and replaces the commas with semicolons

    :param input:       Original summary file from DEXSeq (path)
    :param output:      Output file (path)
    :return:
    """
    with open(input, 'r') as f:
        filedata = f.read()

    # e.g.: c("ENST000000494424", "ENST00000373020") -> ENST000000494424;ENST00000373020
    # remove whitespaces that are included in c(...)
    filedata = filedata.replace('c(', '')
    filedata = filedata.replace(')', '')
    filedata = filedata.replace('", "', '";"')
    filedata = filedata.replace('"', '')
    filedata = filedata.replace(",\n", "")
    filedata = filedata.replace(", \n", "")
    filedata = filedata.replace(",  \n", "")

    with open(output, 'w') as file:
        file.write(filedata)


def extract_significant_results(input_file, output_file, p_value_cutoff=0.05, top_x_select=1000):
    """
    Extracts the significant results from the DEXSeq output file
    :param input_file:
    :param output_file:
    :param p_value_cutoff:  p-value cutoff value
    :return:
    """
    # Read in the results
    results = pd.read_csv(input_file, low_memory=False)

    # Remove count data (drop columns that start with "countData.")
    results = results.loc[:, ~results.columns.str.startswith('countData.')]

    # Remove entries, where the p-value is not significant
    # 1. Remove NA values
    results = results[results['padj'].notna()]      # "NA" values are not significant
    # 2. Remove entries with p-value > 0.05
    results = results[results['padj'] < p_value_cutoff]    # p-value > 0.05 are not significant

    # Sort the results by p-value
    results = results.sort_values(by=['padj'])

    # Select the top x results
    results = results.head(top_x_select)

    # Write the results to a file
    results.to_csv(output_file, index=False)


def integrate_gene_names(result_file, gene_mapping_file, output_file):
    """
    Replace first two columns (gene/group IDs) with gene names

    :param result_file:
    :param gene_mapping_file:
    :param output_file:
    :return:
    """
    # Read in the results
    results = pd.read_csv(result_file, low_memory=False)
    results['ensembl_gene_id'] = results['groupID'].str.split('+').str[0]

    # Get gene names from Ensembl IDs
    ensembl_mappings_df = pd.read_csv(gene_mapping_file, low_memory=False)

    # Merge the two dataframes
    results = pd.merge(results, ensembl_mappings_df, on='ensembl_gene_id')
    sorted_results = results.sort_values(by=['padj'])   # Sort the results by p-value -> ascending order

    # Replace the first two columns with the last two columns
    columns = sorted_results.columns.tolist()[-2:] + sorted_results.columns.tolist()[2:-2]
    sorted_results = sorted_results[columns]

    # Save the results to a file
    sorted_results.to_csv(output_file, index=False)


if __name__ == "__main__":
    # Summary file from DEXSeq & Gene mapping file
    snakemake_summary_file = snakemake.input.summary_csv_file
    snakemake_gene_mapping_file = snakemake.params.gene_mapping_file
    top_x_results = snakemake.params.top_x_results

    # Output file
    snakemake_output_file = snakemake.output.filtered_results_file

    # Clean the summary file
    clean_summary_csv_file(snakemake_summary_file, snakemake_output_file)
    # Extract the significant results
    extract_significant_results(snakemake_output_file, snakemake_output_file, top_x_select=top_x_results)
    # Integrate gene names
    integrate_gene_names(snakemake_output_file, snakemake_gene_mapping_file, snakemake_output_file)
