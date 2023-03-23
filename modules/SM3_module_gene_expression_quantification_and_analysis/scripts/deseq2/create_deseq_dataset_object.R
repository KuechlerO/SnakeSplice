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
library("genefilter")
library("AnnotationDbi")
library("org.Hs.eg.db")


# ---------------- Loading read/fragment quantification data from Salmon output ----------------
load_in_salmon_generated_counts <- function(annotation_table_file) {
	# ----------------- 1. Load annotation -----------------
	# Columns: 1. names, 2. files, 3. condition, 4. additional information
	annotation_table <- read.csv(file=annotation_table_file, sep="\t")

	annotation_data <- data.frame(
									names=annotation_table[,"sample_name"],
									files=file.path(annotation_table[,"salmon_results_file"]),
									condition=annotation_table[,"condition"],
									add_info=annotation_table[,"additional_comment"]
									)
	# Replace None in condition column with "Control"
	annotation_data$condition[annotation_data$condition=="None"] <- "Control"
	annotation_data$condition[annotation_data$condition==""] <- "Control"

	# ----------------- 2. Load into Bioconductor experiment objects -----------------
	# Summarized experiment: Imports quantifications & metadata from all samples -> Each row is a transcript
	se <- tximeta(annotation_data)

	# Summarize transcript-level quantifications to the gene level -> reduces row number: Each row is a gene
	# Includes 3 matrices:
	# 1. counts: Estimated fragment counts per gene & sample
	# 2. abundance: Estimated transcript abundance in TPM
	# 3. length: Effective Length of each gene (including biases as well as transcript usage)
	gse <- summarizeToGene(se)

	# ----------------- 3. Load experiments into DESeq2 object -----------------
	# SummarizedExperiment
	# assayNames(gse)   		# Get all assays -> counts, abundance, length, ...
	# head(assay(gse), 3)     	# Get count results for first 3 genes
	# colSums(assay(gse))     	# Compute sums of mapped fragments
	# rowRanges(gse)          	# Print rowRanges: Ranges of individual genes
	# seqinfo(rowRanges(gse))   # Metadata of sequences (chromosomes in our case)

	gse$condition <- as.factor(gse$condition)
	gse$add_info <- as.factor(gse$add_info)

	# Use relevel to make sure untreated is listed first
	gse$condition %<>% relevel("Control")   # Concise way of saying: gse$condition <- relevel(gse$condition, "Control")

	# Construct DESeqDataSet from gse
	if (gse$add_info %>% unique %>% length >1) {
		# Add info column with more than 1 unique value
		print("More than 1 unique value in add_info column")
		# TODO Need to make sure to avoid:
		# the model matrix is not full rank, so the model cannot be fit as specified.
		#   One or more variables or interaction terms in the design formula are linear
		#   combinations of the others and must be removed.
		# dds <- DESeqDataSet(gse, design = ~condition + add_info)
		print("However, simple DESeq2 analysis will be performed without add_info column")
		dds <- DESeqDataSet(gse, design = ~condition)

	} else {
		print("Only 1 unique value in add_info column")
		dds <- DESeqDataSet(gse, design = ~condition)
	}

	return(dds)
}


# ---------------- Loading read/fragment quantification data from RSEM output ----------------
load_in_rsem_generated_counts <- function(annotation_table_file) {
	# ----------------- 1. Load annotation -----------------
	annotation_table <- read.csv(file=annotation_table_file, sep="\t")
	files <- file.path(annotation_table[,"rsem_results_file"])
	# For sample.genes.results: txIn= FALSE & txOut= FALSE
	# For sample.isoforms.results: txIn= TRUE & txOut= TRUE
	# Check: https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html
	txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

	annotation_data <- data.frame(condition=factor(annotation_table[,"condition"]),
									add_info=factor(annotation_table[,"additional_comment"])
						)
	rownames(annotation_data) <- annotation_table[,"sample_name"]

	# Construct DESeqDataSet from tximport
	if (annotation_data$add_info %>% unique %>% length >1) {
		# Add info column with more than 1 unique value
		# dds <- DESeqDataSetFromTximport(txi.rsem, annotation_data, ~condition + add_info)
		dds <- DESeqDataSetFromTximport(txi.rsem, annotation_data, ~condition)
	} else {
		dds <- DESeqDataSetFromTximport(txi.rsem, annotation_data, ~condition)
	}
	return(dds)
}


load_in_kallisto_generated_counts <- function(annotation_table_file) {
	# ----------------- 1. Load annotation -----------------
	annotation_table <- read.csv(file=annotation_table_file, sep="\t")

	files <- file.path(annotation_table[,"kallisto_results_file"])
	txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)

	annotation_data <- data.frame(condition=factor(annotation_table[,"condition"]),
									add_info=factor(annotation_table[,"additional_comment"])
						)
	rownames(annotation_data) <- annotation_table[,"sample_name"]

	# Construct DESeqDataSet from tximport
	if (annotation_data$add_info %>% unique %>% length >1) {
		# Add info column with more than 1 unique value
		# dds <- DESeqDataSetFromTximport(txi.kallisto, annotation_data, ~condition + add_info)
		dds <- DESeqDataSetFromTximport(txi.kallisto, annotation_data, ~condition)
	} else {
		dds <- DESeqDataSetFromTximport(txi.kallisto, annotation_data, ~condition)
	}
	return(dds)
}


# ----------------- Main function -----------------
main_function <- function(){
	threads <- snakemake@threads[[1]]
	register(MulticoreParam(workers=threads))

	# Snakemake variables
	annotation_table_file <- snakemake@input[["annotation_table_file"]]
	output_file <- snakemake@output[["deseq_dataset_r_obj"]]
	count_algorithm <- snakemake@params[["count_algorithm"]]

	# Load annotation table & Salmon data into a DESeq2 object
	if (count_algorithm == "salmon") {
		dds <- load_in_salmon_generated_counts(annotation_table_file)
	} else if (count_algorithm == "kallisto") {
		dds <- load_in_kallisto_generated_counts(annotation_table_file)
	} else if (count_algorithm == "rsem") {
		dds <- load_in_rsem_generated_counts(annotation_table_file)
	} else {
		stop("Count algorithm not supported!")
	}

	# Remove rows that have no or nearly no information about the amount of gene expression
	print(paste(c("Number of rows before filtering out counts with values <1", nrow(dds))))
	keep <- rowSums(counts(dds)) > 1 # Counts have to be greater than 1
	dds <- dds[keep,]
	print(paste(c("Number of rows after filtering out counts with values <1", nrow(dds))))

    # Save deseq dataset object
    saveRDS(dds, output_file)
}


# ----------------- Run main function -----------------
main_function()
