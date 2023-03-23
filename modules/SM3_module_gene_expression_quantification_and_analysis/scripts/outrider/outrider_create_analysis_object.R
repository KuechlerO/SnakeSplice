# For efficient data handling and visualization, we further load data.table, ggplot2, and ggpubr.
# ATTENTION: Parallelization of the analysis is possible, however
# on a HPC cluster, the parallelization is not working properly.
# A hotfix was published on GitHub: https://github.com/gagneurlab/OUTRIDER/issues/11
# TODO: Upgrade from serial to parallel processing

# ------- 1. Load libraries -------------
library("OUTRIDER")
library("annotables")
library("data.table")
library("ggplot2")
library("ggpubr")
library("plotly")

library("BiocParallel")   # For parallelization needed
# BPPARAM = MulticoreParam(snakemake@threads)


create_outrider_data_set <- function(ctsFile_path) {
    ### Create OUTRIDER data set
    # Input: ctsFile_path: File path leading to count file
    # Output: OUTRIDER data set

    ctsTable <- read.table(ctsFile_path, sep="\t", header=FALSE, stringsAsFactors=FALSE)

    countDataMatrix <- as.matrix(ctsTable[-1,-1])   # Extract counts & Ignore first column (geneID) and first row (sample names)
    mode(countDataMatrix) <- "integer"              # Convert to integer
    rownames(countDataMatrix) <- ctsTable[-1,1]     # Set rownames to geneIDs
    colnames(countDataMatrix) <- ctsTable[1,-1]     # Set colnames to sample names

    # Create OutriderDataSet
    ods <- OutriderDataSet(countData=countDataMatrix)
    print("Done creating OutriderDataSet")

    return(ods)
}


filter_outrider_data_set <- function(ods) {
    ### Filter OUTRIDER data set: Remove genes with low counts
    # Input: ods: OUTRIDER data set
    # Output: Filtered OUTRIDER data set

    # --------- 3. Filter out non expressed genes --------
    # filter out non expressed genes
    # minCounts: If True, only genes wit 0 counts in all samples are filtered out.
    # ALTERNATIVE: If one provides also GTF-annotation, then based on FPKM values filtering is applied
    ods <- filterExpression(ods, minCounts=TRUE, filterGenes=FALSE)
    print("Done filtering out non expressed genes")

    # -------- 3.1 Plotting of the filtered data ---------
    # TODO: Might be not applicable since we do not use FPKM values for filtering

    # Plot FPKM distribution across all sample/gene pairs
    # png_file_path <- file.path(plot_output_dir, paste("fpkm_distribution_across_all_samples_and_genes.png"))
    # png(png_file_path)
    # # TODO: This might not work, since we are not working with FPKM values
    # plotFPKM(ods) + theme(legend.position = 'bottom')
    # dev.off()

    # Apply filter
    ods <- ods[mcols(ods)[['passedFilter']]]
    print("Done applying filter")

    # -------- 3.2 Plotting of potential co-variation ---------
    # TODO: Might be not applicable since we do not use FPKM values
    # -> Requires also a sample annotation file

    # # Make heatmap figure bigger
    # options(repr.plot.width=6, repr.plot.height=5)
    # png_file_path = file.path(plot_output_dir, paste("covariation_heatmap.png"))
    # png(png_file_path)
    #
    # # use normalize=FALSE since the data is not yet corrected
    # # use columns from the annotation to add labels to the heatmap
    # plotCountCorHeatmap(ods, colGroups=c("adaptors_file"), rowGroups="condition", normalize=FALSE)
    # dev.off()

    return(ods)
}


# TODO: clean up naming of output files
# TODO: Do not give output dir, but directly the output files?!
main_function <- function() {
    # Counts file
    ctsFile_path <- snakemake@input[["counts_file"]]

    # With estimated size factors
    outrider_obj_with_estimated_size_factors_txt <- snakemake@output[["outrider_obj_with_estimated_size_factors_txt"]]
    outrider_obj_with_estimated_size_factors_rds <- snakemake@output[["outrider_obj_with_estimated_size_factors_rds"]]

    # Final Outrider object
    output_final_outrider_obj_file <- snakemake@output[["outrider_object_file"]]

    # ------- 1. Create OutriderDataSet ---------
    ods <- create_outrider_data_set(ctsFile_path)

    # ------- 2. Filter out non expressed genes ---------
    ods <- filter_outrider_data_set(ods)

    # -------- 3. Run full Outrider pipeline ------------
    # run full OUTRIDER pipeline (control, fit model, calculate P-values)

    # Crash in case where sample groups are not equally sized!
    #     -> Size Factors, which are accounting for sequencing depth,
    #     are differing strongly between sample groups -> Controlling for confounders leads to crash?!
    #     - first remote error: L-BFGS-B benötigt endliche Werte von 'fn' -> filter values?!
    #     - https://github.com/gagneurlab/OUTRIDER/issues/25

    # --------- For debugging purposes: Saves estimated size factors ---------
    print("Investigate estimated size factors")
    ods <- OUTRIDER::estimateSizeFactors(ods)
    saveRDS(ods, outrider_obj_with_estimated_size_factors_rds)
    sink(outrider_obj_with_estimated_size_factors_txt)
    print(OUTRIDER::sizeFactors(ods))
    sink()
    # ------------------------------------------

    print("Start running full OUTRIDER pipeline")
    ods <- OUTRIDER(ods, BPPARAM=SerialParam())
    # Save the OutriderDataSet -> can be used for further analysis via R: readRDS(file)
    saveRDS(ods, output_final_outrider_obj_file)
}

main_function()
