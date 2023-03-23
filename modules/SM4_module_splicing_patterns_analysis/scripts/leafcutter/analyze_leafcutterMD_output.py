import pandas as pd


# filters input_df according to given p-val threshold
# -> Returns filtered dataframe
def filter_dfs(input_df, p_val_threshold=0.01, columns=None):
    """
    Filters input_df according to given p-val threshold
    Sort & extract top 1000 rows
    :param input_df:
    :param p_val_threshold:
    :param columns:
    :return:
    """
    if columns is None or len(columns) == 0:
        columns = input_df.columns

    sub_df = input_df[columns]

    # filtering: only keep rows with p-val < p_val_threshold
    # .any(axis=1) -> keep rows with at least one True value
    output_df = sub_df[(sub_df <= p_val_threshold).any(axis=1)]

    # Sort by minimal p-value
    output_df['min_pvalue'] = output_df.min(axis=1)
    output_df = output_df.sort_values(by=['min_pvalue'])

    # Extract only top 1000 results
    output_df = output_df.head(n=1000)

    return output_df


if __name__ == '__main__':
    # Load files
    all_outlier_introns_pVals_file = snakemake.input["all_outlier_introns_pVals_file"]
    all_outlier_clusters_pVals_file = snakemake.input["all_outlier_clusters_pVals_file"]
    all_outlier_effSize_file = snakemake.input["all_outlier_effSize_file"]

    # Output files
    all_filtered_introns_file = snakemake.output["all_filtered_introns_file"]
    condition_filtered_introns_file = snakemake.output["condition_filtered_introns_file"]
    all_filtered_clusters_file = snakemake.output["all_filtered_clusters_file"]
    condition_filtered_clusters_file = snakemake.output["condition_filtered_clusters_file"]

    # Load sample names for affected samples
    sample_ids = snakemake.params.patient_sample_ids
    pvalue_threshold = snakemake.params.pvalue_threshold

    # Load dataframes of LeafcutterMD results
    introns_df = pd.read_csv(all_outlier_introns_pVals_file, sep='\t')
    clusters_df = pd.read_csv(all_outlier_clusters_pVals_file, sep='\t')

    # Intron assessment results
    filter_dfs(introns_df, p_val_threshold=pvalue_threshold).to_csv(all_filtered_introns_file, sep='\t')
    filter_dfs(introns_df, p_val_threshold=pvalue_threshold, columns=sample_ids).to_csv(condition_filtered_introns_file, sep='\t')

    # Cluster assessment results
    filter_dfs(clusters_df, p_val_threshold=pvalue_threshold).to_csv(all_filtered_clusters_file, sep='\t')
    filter_dfs(clusters_df, p_val_threshold=pvalue_threshold, columns=sample_ids).to_csv(condition_filtered_clusters_file, sep='\t')
