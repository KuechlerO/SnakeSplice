# Gene-level expression analysis using DESeq2
# http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#preparing-quantification-input-to-deseq2

# Gene-level expression analysis based on transcript abundance quantifiers
# 1. Corrects for potential changes in gene length across samples (e.g. differential isoform usage)
# 2. Methods are fast and require less memory and disk usage compared to alignment methods
# 3. Possible avoidance of discarding fragments that can align to multiple genes with homologous sequence

library("BiocParallel")

# Load data
library("tximeta")	# Import transcript quantification data from Salmon
library("tximport")	# Import transcript-level quantification data from Kaleidoscope, Sailfish, Salmon, Kallisto, RSEM, featureCounts, and HTSeq
library("rhdf5")
library("SummarizedExperiment")

library("magrittr")	# Pipe operator
library("DESeq2")		# Differential gene expression analysis

# Plotting libraries
library("pheatmap")
library("RColorBrewer")
library("ggplot2")

# Mixed
library("PoiClaClu")
library("glmpca")
library("apeglm")
library("ashr")
library("genefilter")
library("AnnotationDbi")
library("org.Hs.eg.db")

library("ReportingTools")	# For creating HTML reports



# ---------------- Helper functions ----------------
savely_create_deseq2_object <- function(dds) {
	### Function to create a DESeq2 object -> Handles errors that can appear due to parallelization
	### Input: dds object (DESeq dataset object)
	### Output: dds object (DESeq2 object)

	# ------------- 5. Run the differential expression analysis ---------------
	# The respective steps of this function are printed out
	# 1. Estimation of size factors: Controlling for differences
	# in the sequencing depth of the samples
	# 2. Estimation of dispersion values for each gene & fitting a generalized
	# linear model
	print("Creating DESeq2 object")

	# Try to create the DESeq2 object with results in Parallel
	create_obj_parallelized <- function(){
		print("Creating DESeq2 object in parallel")
		dds <- DESeq2::DESeq(dds, parallel=TRUE)
		return(list("dds"=dds, "run_parallel"=TRUE))
	}
	# Try to create the DESeq2 object with results in Serial
	create_obj_not_parallelized <- function(error){
		print("Error in parallelized DESeq2 object creation"); print(error)
		print("Creating DESeq2 object not in parallel")
		dds <- DESeq2::DESeq(dds, parallel=FALSE)
		return(list("dds"=dds, "run_parallel"=FALSE))
	}

	result_list <- tryCatch(create_obj_parallelized(), error=create_obj_not_parallelized)
	print("DESeq2 object created!")
	return(result_list)
}

rename_rownames_with_ensembl_id_matching <- function(dds_object, input_algorithm) {
  	"
	Extracts Ensembl-Gene-IDs from rownames of SummarizedExperiment object and
	renames rownames with Ensembl-Gene-IDs.
	"
  print("Gene annotations")
  if (input_algorithm == "salmon") {
    # Ensembl-Transcript-IDs at first place
    gene_ids_in_rows <- substr(rownames(dds_object), 1, 15)
  }
  else if (input_algorithm == "kallisto") {
    # Ensembl-Transcript-IDs at second place (delimeter: "|")
    gene_ids_in_rows <- sapply(rownames(dds_object), function(x) strsplit(x, '\\|')[[1]], USE.NAMES=FALSE)[2,]
    gene_ids_in_rows <- sapply(gene_ids_in_rows, function(x) substr(x, 1, 15), USE.NAMES=FALSE)
  }
  else {
    stop("Unknown algorithm used for quantification")
  }

  # Set new rownames
  rownames(dds_object) <- gene_ids_in_rows
  return(dds_object)
}


add_gene_symbol_and_entrez_id_to_results <- function(result_object, with_entrez_id=FALSE) {
	"
	Adds gene symbols and entrez-IDs to results object.
	"
	gene_ids_in_rows <- rownames(result_object)

	# Add gene symbols
	# Something breaks here when setting a new column name
	result_object$symbol <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
												 keys=gene_ids_in_rows,
												 column="SYMBOL",
												 keytype="ENSEMBL",
												 multiVals="first")
	if (with_entrez_id) {
		# Add ENTREZ-ID
		result_object$entrez <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
													 keys=gene_ids_in_rows,
													 column="ENTREZID",
													 keytype="ENSEMBL",
													 multiVals="first")
	}

	return(result_object)
}


# ---------------- DESeq2 analysis ----------------
explore_deseq2_results <- function(dds, false_discovery_rate, output_file_paths, run_parallel=FALSE,
                                   used_algorithm) {
	# Results: Metadata
	# 1. baseMean: Average/Mean of the normalized count values divided by size factors, taken over ALL samples
	# 2. log2FoldChange: Effect size estimate. Change of gene's expression
	# 3. lfcSE: Standard Error estimate for log2FoldChange
	# 4. Wald statistic results
	# 5. Wald test p-value ->  p value indicates the probability that a fold change as strong as the observed one, or even stronger, would be seen under the situation described by the null hypothesis.
	# 6. BH adjusted p-value
	print("Creating DESeq2 results object")
	results_obj <- results(dds, alpha=false_discovery_rate, parallel=run_parallel)
	capture.output(summary(results_obj), file=output_file_paths[1])


	# ------------------ 6. Plotting results --------------------
	# Contrast usage
	# TODO: Failed... -> Remove
	# print("Plotting results")
	# chosen_contrast <- tail(resultsNames(results_obj), n=1)	     # get the last contrast: Comparison of states
	# print("resultsNames(results_obj)"); print(resultsNames(results_obj))
	# print("chosen_contrast"); print(chosen_contrast)

	# ------------ 6.1 MA plot without shrinking --------------
	# - M: minus <=> ratio of log-values -> log-Fold-change on Y-axis
	# - A: average -> Mean of normalized counts on X-axis

	# res.noshr <- results(dds, contrast=chosen_contrast, parallel=run_parallel)
	res.no_shrink <- results(dds, parallel=run_parallel)
	jpeg(output_file_paths[2], width=800, height=800)
	DESeq2::plotMA(res.no_shrink, ylim = c(-5, 5), main="MA plot without shrinkage")
	dev.off()

	# ------------ 6.2 MA plot with apeGLM shrinking --------------
	# apeglm method for shrinking coefficients
	# -> which is good for shrinking the noisy LFC estimates while
	# giving low bias LFC estimates for true large differences

	# TODO: apeglm requires coefficients. However, resultsNames(results_obj) does not return any coefficients...
	# res <- lfcShrink(dds, coef=chosen_contrast, type="apeglm", parallel=run_parallel)      # Pass contrast and shrink results
	# Use ashr as shrinkage method
	res <- DESeq2::lfcShrink(dds, res=res.no_shrink, type="ashr", parallel=run_parallel)
	jpeg(output_file_paths[3], width=800, height=800)
	DESeq2::plotMA(res, ylim = c(-5, 5), main="MA plot with ashr shrinkage")
	dev.off()

	# ------------ 6.3 Plot distribution of p-values in histogram --------------
	jpeg(output_file_paths[4], width=800, height=800)
	hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white",
		main="Histogram of distribution of p-values (non-adjusted)", xlab="p-value", ylab="Frequency")
	dev.off()

	jpeg(output_file_paths[5], width=800, height=800)
	hist(res$padj[res$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white",
		main="Histogram of distribution of p-values (adjusted)", xlab="p-value", ylab="Frequency")
	dev.off()

	# ------------- 6.4 Gene clustering -----------------
	print("Plotting results: Gene clustering")
	# Gene clustering -> Heatmap of divergence of gene's expression in comparison to average over all samples
	# Transform count results to reduce noise for low expression genes
	transformed_dds <- DESeq2::rlog(dds, blind=FALSE)

	# Get top Genes -> with most variance in VSD-values/rlog-transformed counts
	topVarGenes <- head(order(genefilter::rowVars(SummarizedExperiment::assay(transformed_dds)), decreasing=TRUE), 20)
	mat  <- SummarizedExperiment::assay(transformed_dds)[ topVarGenes, ]
	mat  <- mat - rowMeans(mat)   # difference to mean expression
	# Transform row names to gene symbols
	rownames(mat) <- add_gene_symbol_and_entrez_id_to_results(mat)$symbol
	# Additional annotations
	anno <- as.data.frame(SummarizedExperiment::colData(transformed_dds)[, c("condition", "add_info")])
	# Create plot
	jpeg(output_file_paths[6], width=800, height=800)
	pheatmap::pheatmap(mat, annotation_col=anno,
			 main="Divergence in gene expression in comparison to average over all samples")
	dev.off()

	# ---------- 7. Gene annotations --------------
	res <- add_gene_symbol_and_entrez_id_to_results(res)
	resOrdered <- res[order(res$pvalue),]					# Sort results by p-value

	# Exporting results
	resOrderedDF <- as.data.frame(resOrdered)
	write.csv(resOrderedDF, file=output_file_paths[7])
}



# ----------------- Main function -----------------
main_function <- function(){
	threads <- snakemake@threads[[1]]
	register(MulticoreParam(workers=threads))

	# Snakemake variables
	deseq_dataset_obj <- snakemake@input[["deseq_dataset_r_obj"]]
	output_file_paths <- snakemake@params[["output_file_paths"]]
	# For gene-ID matching: Used in rename_rownames_with_ensembl_id_matching()
	used_algorithm <- snakemake@params["used_algorithm"]
	# Adjusted p-value threshold
	false_discovery_rate <- 0.05

    # Load deseq dataset object
    dds <- readRDS(deseq_dataset_obj)

	# Create DESeq2 results object
	print("Creating DESeq2 results object")
	result_list <- savely_create_deseq2_object(dds)
	deseq2_obj <- result_list$dds
	run_parallel <- result_list$run_parallel

	# Rename rows
	deseq2_obj <- rename_rownames_with_ensembl_id_matching(deseq2_obj, used_algorithm)

	# Run statistical analysis
	print("Running statistical analysis")
	explore_deseq2_results(deseq2_obj, false_discovery_rate, output_file_paths, run_parallel=run_parallel,
	                       used_algorithm=used_algorithm)
}


# ----------------- Run main function -----------------
main_function()
