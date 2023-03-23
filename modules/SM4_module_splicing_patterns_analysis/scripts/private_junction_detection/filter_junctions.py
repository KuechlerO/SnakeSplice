import pandas as pd


def filter_junctions(junction_collection_file, control_samples, condition_samples,
                     only_ctr_junc_file, only_cond_junc_file,
                     max_contrast=0.0,
                     in_all_samples=True):
    """
    Merges the junctions from the different samples into one file
    :param junction_collection_file:    File containing the junctions from all samples
    :param control_samples:             List of samples from condition 1
    :param condition_samples:           List of samples from condition 2
    :param only_ctr_junc_file:          Output file for junctions only in condition 1
    :param only_cond_junc_file:         Output file for junctions only in condition 2
    :param max_contrast:                Part of total nr of reads of all samples
                                        over all conditions, that is allowed to appear.
                                        E.g. 0.2, where total number of reads is 100 -> then max a total of 20 (-> 20%)
                                         reads are allowed to be in the contra-condition
    :param in_all_samples:              If True, junctions have to be in all samples
    :return:
    """

    # Column-names: Each line is a exon-exon junction
    total_junction_df = pd.read_csv(junction_collection_file, low_memory=False, sep="\t")

    # Select first 5 info columns and respective sample columns
    print("Control samples: ", control_samples)
    print("Condition samples: ", condition_samples)
    selected_columns = total_junction_df.columns[:5].tolist() + control_samples + condition_samples
    total_junction_df_selected = total_junction_df[selected_columns]

    # Collect total sum & compute maximum contrast
    total_junction_df_selected["total_sum"] = (total_junction_df_selected[control_samples].sum(axis=1)
                                               + total_junction_df_selected[condition_samples].sum(axis=1))
    total_junction_df_selected["max_contrast"] = total_junction_df_selected["total_sum"]*max_contrast

    # ------------- 1. Filter only control junctions -------------
    only_control_junctions_df = apply_filtering(total_junction_df_selected, control_samples, condition_samples, in_all_samples)
    only_control_junctions_df = only_control_junctions_df.sort_values(by="total_sum", ascending=False)

    # ------------- 2. Filter only condition junctions -------------
    only_condition_junctions_df = apply_filtering(total_junction_df_selected, condition_samples, control_samples, in_all_samples)
    only_condition_junctions_df = only_condition_junctions_df.sort_values(by="total_sum", ascending=False)

    # Save the results in output dir
    only_control_junctions_df.to_csv(only_ctr_junc_file, index=False, sep="\t")
    only_condition_junctions_df.to_csv(only_cond_junc_file, index=False, sep="\t")


def apply_filtering(input_df, samples_1, samples_2, in_all_samples_bool):
    """
    Filters the junctions in the input_df
    1. Only rows where at least one read in all samples of samples_1 are collected
    2. If in_all_samples_bool is True, then only junctions in all samples_1 are collected
    3. Only junctions where the sum of all samples_1 is higher than the sum of all samples_2 (depending on max_contrast)

    :param input_df:            Input dataframe
    :param samples_1:           List of samples from condition 1
    :param samples_2:           List of samples from condition 2
    :param in_all_samples_bool:     If True, junctions have to be in all samples of condition 1
    :return:
    """

    # 1. Select only rows where at least one sample has a value > 0
    df_with_s1_jcts = input_df[input_df[samples_1].sum(axis=1) > 0]

    # 2. Every junction has to be in all s1 samples
    if in_all_samples_bool:
        for sample in samples_1:
            df_with_s1_jcts = df_with_s1_jcts[df_with_s1_jcts[sample] > 0]

    # 3. Filtering depending on the max_contrast column
    df_with_s1_jcts = df_with_s1_jcts[(df_with_s1_jcts[samples_1].sum(axis=1) > df_with_s1_jcts["max_contrast"])
                                      & (df_with_s1_jcts[samples_2].sum(axis=1) <= df_with_s1_jcts["max_contrast"])]

    return df_with_s1_jcts


if __name__ == "__main__":
    # Input files
    snakemake_junction_collection_file = snakemake.input.junction_collection_file
    # params
    snakemake_control_samples = snakemake.params.control_samples
    snakemake_condition_samples = snakemake.params.condition_samples
    # Output file
    snakemake_only_control_junctions_file = snakemake.output.only_control_junctions_file
    snakemake_only_condition_junctions_file = snakemake.output.only_condition_junctions_file

    # params
    snakemake_max_contrast = snakemake.params.max_contrast

    # Filter junctions
    filter_junctions(snakemake_junction_collection_file, snakemake_control_samples, snakemake_condition_samples,
                     snakemake_only_control_junctions_file, snakemake_only_condition_junctions_file,
                     snakemake_max_contrast)
