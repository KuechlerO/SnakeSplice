# This script runs the differential expression analysis using the DEXSeq package

library("BiocParallel")   # For parallelization needed
library("DEXSeq")         # For differential expression analysis

library("dplyr")
source(snakemake@params[["load_subreadOutput_script"]])      # Source script to import Subread output


runDexseqAnalysis <- function(countFile, input_flattened_gtf_file, sampleTable,
                              dexseq_results_object_file, output_csv_file, BPPARAM) {
    # ------ A. Load count data into table -------
    # 3. argument: specifies the formula for creating a linear regression model -> interaction between condition and exon
    # Using this formula, we are interested in differences in exon usage due to the “condition” variable changes.
    dxd <- DEXSeqDataSetFromFeatureCounts(
        countFile,
        flattenedfile = input_flattened_gtf_file,
        sampleData = sampleTable
    )

#     dxd = DEXSeqDataSetFromHTSeq(
#        countFile,
#        sampleData=sampleTable,
#        design= ~ sample + exon + condition:exon,
#        flattenedfile=input_flattened_gtf_file
#     )

    # ------ B. Normalization -------
    dxd <- estimateSizeFactors(dxd)      # Normalization -> Uses same method as DESeq2

    # -------- 4.3 Dispersion estimation ---------
    # To test for differential exon usage, we need to estimate the variability of the data.
    # This is necessary to be able to distinguish technical and biological variation (noise) from real
    # effects on exon usage due to the different conditions.
    dxd <- estimateDispersions(dxd, BPPARAM=BPPARAM)

    # --------- C. Testing for differential exon usage ----------------
    # For each gene, DEXSeq fits a generalized linear model with the formula
    # ~sample + exon + condition:exon
    # and compares it to the smaller model (the null model)
    # ~ sample + exon.

    # exon: Factor with 2 levels: this and others
    # Explanation for linear models in R: For every coefficient to add into formula: use symbol "+"
    # Interactions are separated by colon -> condition:exon -> interpreted as multiplication term

    # Testing: The deviances of both fits are compared using a χ2-distribution, providing a p value.
    # Based on this p-value, we can decide whether the null model is sufficient to explain the data,
    # or whether it may be rejected in favour of the alternative model, which contains an interaction
    # coefficient for condition:exon. The latter means that the fraction of the gene’s reads that fall
    # onto the exon under the test differs significantly between the experimental conditions.
    dxd <- testForDEU(dxd, BPPARAM=BPPARAM)

    # --------- D. Compute exon fold change ----------------
    # Compute exon fold change numbers with formula:
    # count ~ condition + exon + condition:exon
    dxd <- estimateExonFoldChanges(dxd, fitExpToVar="condition", BPPARAM=BPPARAM)

    # ------- E. Results ------------
    # Summarize results, save R-object and write result summary to file
    dxr1 <- DEXSeqResults(dxd)

    print("Save DEXSeq results object in R-file")
    saveRDS(dxr1, file=dexseq_results_object_file)

    print("Save summary of results in CSV-file")
    write.csv(dxr1, file=output_csv_file, row.names=TRUE)
}

# -------------------- Main function --------------------
main <- function(){
    # ----------------- 1. Load snakemake variables -----------------
    # inputs
    input_flattened_gtf_file <- snakemake@input[["flattened_gtf_file"]]
    input_exon_counting_bin_file <- snakemake@input[["exon_counting_bin_file"]]

    # params
    input_sample_ids <- snakemake@params[["sample_ids"]]
    input_sample_conditions <- snakemake@params[["sample_conditions"]]
    threads <- snakemake@threads

    # outputs
    dexseq_results_object_file <- snakemake@output[["dexseq_results_object_file"]]
    output_csv_file <- snakemake@output[["result_summary_csv_file"]]

    # ----------------- 2. Prepare analysis -----------------
    # ------ 2.1 Set number of threads --------
    BPPARAM <- MulticoreParam(threads)
    # register(MulticoreParam(threads))

    # ------ 2.2 Load annotation data into table -------
    # Table: One row for each library (each sample)
    # Columns: For all relevant information -> covariates
    # If only one covariant, it has to be named "condition"!
    overallSampleTable <- data.frame(
        row.names = input_sample_ids,
        condition = input_sample_conditions
       )
    print("Sample table:")
    print(overallSampleTable)

    # --------------- 2.3 Extract only needed columns from Subread output ----------------
    # Extract only needed columns from Subread output
    # counts_file <- fread(input_exon_counting_bin_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    counts_file <- read.table(input_exon_counting_bin_file, header=TRUE, check.names=FALSE, sep="\t", stringsAsFactors=FALSE)
    keep_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length", input_sample_ids)
    # Extract only needed columns from Subread output and ensure order of columns
    subset_counts <- counts_file[names(counts_file) %in% keep_cols][keep_cols]
    # Write selected columns to temporary file
    tmp_subset_counts_file <- tempfile(fileext = ".tsv")
    write.table(subset_counts, file=tmp_subset_counts_file, sep="\t", quote=FALSE, row.names=FALSE)


    # ----------------- 3. Run analysis -----------------
    # Run dexseq analysis
    print("Run analysis")
    runDexseqAnalysis(tmp_subset_counts_file, input_flattened_gtf_file, overallSampleTable,
                  dexseq_results_object_file, output_csv_file, BPPARAM)
}

# Run main function
main()

