import pandas as pd


def load_simple_result_file(input_file_path, input_gene_name_col, input_adjusted_pval_col, tool_name):
    """
    Simply loads the given file (automatic separator detection) and
    returns a dataframe with the following columns:
    - gene_name
    - <tool_name> adjusted-p-value
    - <tool_name> ranking
    """

    df = None

    if input_file_path:
        # New column names
        ranking_col = tool_name + ": ranking"
        adjusted_pval_col = tool_name + ": adjusted-p-value"

        # Load file -> automatically detect separator
        df = pd.read_csv(str(input_file_path), sep=None)
        df[ranking_col] = df.index +1   # +1 because index starts at 0

        # Rename gene name column
        df = df.rename(columns={input_gene_name_col: "gene_name"})

        # Rename adjusted p-value column
        if input_adjusted_pval_col:
            df = df.rename(columns={input_adjusted_pval_col: adjusted_pval_col})
            df = df[["gene_name", adjusted_pval_col, ranking_col]]
        else:   # no adjusted p-value column given -> So do not use it
            df = df[["gene_name", ranking_col]]

    else:
        print("No input file for " + tool_name + " provided. Skipping...")

    return df


def load_result_files():
    """
    Loads result files from all different tools of the given workflow.
    Returns an array of the resulting dataframes.

    Attention: Makes use of snakemake.input and snakemake.params.
    :return:
    """
    output_result_dfs = []

    # 1. Load pjd results
    try:
        # PJD condition
        if snakemake.input.pjd_condition:
            pjd_condition_df = load_simple_result_file(snakemake.input.pjd_condition,
                                                       snakemake.params.pjd_gene_col_name,
                                                       None, "PJD-condition")
            output_result_dfs.append(pjd_condition_df)
        # PJD control
        if snakemake.input.pjd_control:
            pjd_control_df = load_simple_result_file(snakemake.input.pjd_control, snakemake.params.pjd_gene_col_name,
                                                     None, "PJD-control")
            output_result_dfs.append(pjd_control_df)

    except AttributeError as e:
        print("1. No input file for PJD provided. Skipping...")
        print(e)

    # 2. Load Leafcutter results
    try:
        if snakemake.input.leafcutter_results:
            leafcutter_df = load_simple_result_file(snakemake.input.leafcutter_results, snakemake.params.leafcutter_gene_col_name,
                                                    snakemake.params.leafcutter_adjusted_pval_col_name, "Leafcutter")
            output_result_dfs.append(leafcutter_df)
    except AttributeError as e:
        print("2. No input file for Leafcutter provided. Skipping...")
        print(e)

    # 3. Load fraser results
    try:
        if snakemake.input.fraser_results:
            fraser_df = pd.read_csv(str(snakemake.input.fraser_results), sep=",")

            # select only condition samples -> params.samples_with_condition
            fraser_df = fraser_df[fraser_df["sampleID"].isin(snakemake.params.samples_with_condition)]

            fraser_df["FRASER: ranking"] = fraser_df.index
            fraser_df = fraser_df.rename(columns={snakemake.params.fraser_gene_col_name: "gene_name"})
            fraser_df = fraser_df.rename(columns=
                                         {snakemake.params.fraser_adjusted_pval_col_name: "FRASER: adjusted-p-value"})
            fraser_df = fraser_df[["gene_name", "FRASER: adjusted-p-value", "FRASER: ranking"]]
            output_result_dfs.append(fraser_df)
    except AttributeError as e:
        print("3. No input file for FRASER provided. Skipping...")
        print(e)

    # 4. Load dexseq results
    try:
        if snakemake.input.dexseq_results:
            dexseq_df = load_simple_result_file(snakemake.input.dexseq_results,
                                                snakemake.params.dexseq_gene_col_name,
                                                snakemake.params.dexseq_adjusted_pval_col_name, "DEXSeq")
            output_result_dfs.append(dexseq_df)
    except AttributeError as e:
        print("2. No input file for DEXSeq provided. Skipping...")
        print(e)

    # 5. Load rMATS results
    try:
        # 1. Load rMATS results for A3SS
        if snakemake.input.rmats_results_a3ss_jcec:
            rmats_a3ss_jcec_df = load_simple_result_file(snakemake.input.rmats_results_a3ss_jcec,
                                                         snakemake.params.rmats_gene_col_name,
                                                         snakemake.params.rmats_adjusted_pval_col_name, "rMATS-A3SS")
            output_result_dfs.append(rmats_a3ss_jcec_df)

        # 2. Load rMATS results for A5SS
        if snakemake.input.rmats_results_a5ss_jcec:
            rmats_a5ss_jcec_df = load_simple_result_file(snakemake.input.rmats_results_a5ss_jcec,
                                                         snakemake.params.rmats_gene_col_name,
                                                         snakemake.params.rmats_adjusted_pval_col_name, "rMATS-A5SS")
            output_result_dfs.append(rmats_a5ss_jcec_df)

        # 3. Load rMATS results for MXE
        if snakemake.input.rmats_results_mxe_jcec:
            rmats_mxe_jcec_df = load_simple_result_file(snakemake.input.rmats_results_mxe_jcec,
                                                        snakemake.params.rmats_gene_col_name,
                                                        snakemake.params.rmats_adjusted_pval_col_name, "rMATS-MXE")
            output_result_dfs.append(rmats_mxe_jcec_df)

        # 4. Load rMATS results for RI
        if snakemake.input.rmats_results_ri_jcec:
            rmats_ri_jcec_df = load_simple_result_file(snakemake.input.rmats_results_ri_jcec,
                                                       snakemake.params.rmats_gene_col_name,
                                                       snakemake.params.rmats_adjusted_pval_col_name, "rMATS-RI")
            output_result_dfs.append(rmats_ri_jcec_df)

        # 5. Load rMATS results for SE
        if snakemake.input.rmats_results_se_jcec:
            rmats_se_jcec_df = load_simple_result_file(snakemake.input.rmats_results_se_jcec,
                                                       snakemake.params.rmats_gene_col_name,
                                                       snakemake.params.rmats_adjusted_pval_col_name, "rMATS-SE")
            output_result_dfs.append(rmats_se_jcec_df)
    except AttributeError as e:
        print("5. No input file for rMATS provided. Skipping...")
        print(e)

    # Finally return the list of dataframes
    return output_result_dfs


def merge_results(output_result_dfs):
    """"
    Load result files from different tools and merge them into one dataframe.
    ATTENTION: Remove empty gene_name rows, since otherwise a memory overload occurs during merging.
    """
    merged_df = output_result_dfs[0]
    for df in output_result_dfs[1:]:
        merged_df = merged_df.merge(df, on="gene_name", how="outer")

    # count non-empty cells per row for chosen columns
    # Select columns that contain "ranking" in the column name
    ranking_cols_without_rmats = [col for col in merged_df.columns if "ranking" in col and "rMATS" not in col]
    ranking_cols_with_rmats = [col for col in merged_df.columns if "ranking" in col]
    # A: Agreement score without rMATS
    if len(ranking_cols_without_rmats) == 0:
        merged_df["Agreement Sum without rMATS"] = 0
    else:
        merged_df["Agreement Sum without rMATS"] = merged_df[ranking_cols_without_rmats].notnull().sum(axis=1)
    # B: Agreement score with rMATS
    merged_df["Agreement Sum with rMATS"] = merged_df[ranking_cols_with_rmats].notnull().sum(axis=1)

    # sort by detection sum
    merged_df = merged_df.sort_values(by=["Agreement Sum without rMATS", "Agreement Sum with rMATS"],
                                      ascending=[False, False])

    # replace NaN with empty string
    merged_df = merged_df.fillna("None")
    # Place col "gene_name", "Agreement Sum without rMATS", "Agreement Sum with rMATS" at the beginning
    cols = ["gene_name", "Agreement Sum without rMATS", "Agreement Sum with rMATS"] + \
           list([col for col in merged_df.columns if col not in ["gene_name", "Agreement Sum without rMATS",
                                                                 "Agreement Sum with rMATS"]
           ])
    merged_df = merged_df[cols]

    return merged_df


if __name__ == "__main__":
    # input
    input = snakemake.input
    # output
    output = snakemake.output
    # params
    params = snakemake.params

    # load results
    print("Loading result files", flush=True)
    result_dfs = load_result_files()    # Makes use of snakemake.input and snakemake.params
    # assert that list is not empty
    assert result_dfs, "No results loaded. Check input files."

    # clean results_dfs
    reduced_result_dfs = []
    print("Cleaning result dataframes", flush=True)
    for df in result_dfs:
        # 1. Remove rows where gene_name is "None", or ".", or "NA"
        df = df[~df["gene_name"].isin(["None", ".", "NA"])]
        # 2. Remove rows where gene_name is NaN
        df = df.dropna(subset=["gene_name"])
        # 3. Remove duplicate entries where gene_name is not unique
        df = df.drop_duplicates(subset="gene_name", keep="first")

        reduced_result_dfs.append(df)

    # merge results
    print("Now merging results", flush=True)
    merged_df = merge_results(reduced_result_dfs)
    print("Merging done", flush=True)

    # write output
    merged_df.to_csv(output[0], sep="\t", index=False)

