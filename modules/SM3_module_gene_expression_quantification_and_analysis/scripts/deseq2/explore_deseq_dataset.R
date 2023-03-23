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



# ---------------- DESeq2 explorative analysis ----------------
run_deseq2_explorative_analysis <- function(dds, output_files) {

	# ----------- 4.2 Variance stabilizing transformation and the rlog -------------
	# Problem: PCA depends mostly on points with highest variance
	# -> For gene-counts: Genes with high expression values, and therefore high variance
	# are the ones the PCA is mostly depending on
	# Solution: Apply stabilizing transformation to variance
	# -> transform data, so it becomes more homoskedastic (expected amount of variance the same across different means)
	# 1. Variance stabilizing transformation: VST-function -> fast for large datasets (> 30n)
	# 2. Regularized-logarithm transformation or rlog -> Works well on small datasets (< 30n)

	# The transformed values are no longer counts, and are stored in the assay slot.
	# 1. VST
	# transformed_dds <- vst(dds, blind = FALSE)
	# head(assay(transformed_dds), 3)

	# 2. rlog
	transformed_dds <- rlog(dds, blind=FALSE)

	# ----------- A. Sample distances -------------
	# ----------- A.1 Euclidian distances -------------
	# Sample distances -> Assess overall similarity between samples
	# dist: takes samples as rows and genes as columns -> we need to transpose
	sampleDists <- dist(t(assay(transformed_dds)))

	# Heatmap of sample-to-sample distances using the transformed values
	# Uses euclidian distance between samples
	sampleDistMatrix <- as.matrix(sampleDists)
	rownames(sampleDistMatrix) <- paste(transformed_dds$names, transformed_dds$condition, sep = " - " )
	colnames(sampleDistMatrix) <- paste(transformed_dds$names, transformed_dds$condition, sep = " - " )
	colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
	jpeg(output_files[1], width=800, height=800)
	pheatmap(sampleDistMatrix,
			 clustering_distance_rows = sampleDists,
			 clustering_distance_cols = sampleDists,
			 col = colors,
			 main = "Heatmap of sample-to-sample distances (Euclidian) after normalization")
	dev.off()


	# ----------- A.2 Poisson distances -------------
	# Use Poisson distance
	# -> takes the inherent variance structure of counts into consideration
	# The PoissonDistance function takes the original count matrix (not normalized) with samples as rows instead of
	# columns -> so we need to transpose the counts in dds.
	poisd <- PoissonDistance(t(counts(dds)))

	# heatmap
	samplePoisDistMatrix <- as.matrix(poisd$dd)
	rownames(samplePoisDistMatrix) <- paste(transformed_dds$names, transformed_dds$condition, sep=" - ")
	colnames(samplePoisDistMatrix) <- paste(transformed_dds$names, transformed_dds$condition, sep=" - ")
	jpeg(output_files[2], width=800, height=800)
	pheatmap(samplePoisDistMatrix,
			 clustering_distance_rows = poisd$dd,
			 clustering_distance_cols = poisd$dd,
			 col = colors,
			 main = "Heatmap of sample-to-sample distances (Poisson) without normalization")
	dev.off()

	# ------------ 4.4 PCA plot -------------------

	# ----------------- 4.4.1 Custom PCA plot --------------
	# Build own plot with ggplot -> to distinguish subgroups more clearly
	# Each unique combination of treatment and cell-line has unique color
	# Use function that is provided with DeSeq2
	pcaData <- plotPCA(transformed_dds, intgroup = c("condition", "add_info"), returnData=TRUE)
	percentVar <- round(100 * attr(pcaData, "percentVar"))

	print("Creating custom PCA plot")
	jpeg(output_files[3], width=800, height=800)
	customPCAPlot <- ggplot(pcaData, aes(x=PC1, y=PC2, color=condition, shape=add_info, label=name)) +
		geom_point(size =3) +
		geom_text(check_overlap=TRUE, hjust=0, vjust=1) +
		xlab(paste0("PC1: ", percentVar[1], "% variance")) +
		ylab(paste0("PC2: ", percentVar[2], "% variance")) +
		coord_fixed() +
		ggtitle("PCA on transformed (rlog) data with subgroups (see shapes)")
	print(customPCAPlot)
	dev.off()

	# ----------------- 4.4.2 Generalized PCA plot --------------
	# Generalized PCA: Operates on raw counts, avoiding pitfalls of normalization
	print("Creating generalized PCA plot")
	gpca <- glmpca(counts(dds), L=2)
	gpca.dat <- gpca$factors
	gpca.dat$condition <- dds$condition
	gpca.dat$add_info <- dds$add_info

	jpeg(output_files[4], width=800, height=800)
	generalizedPCAPlot <- ggplot(gpca.dat, aes(x=dim1, y=dim2, color=condition, shape=add_info,
											   label=rownames(gpca.dat))) +
	  	geom_point(size=2) +
		geom_text(check_overlap=TRUE, hjust=0.5,vjust=1) +
		coord_fixed() +
		ggtitle("glmpca - Generalized PCA of samples")
	print(generalizedPCAPlot)
	dev.off()
}


# ----------------- Main function -----------------
main_function <- function(){
	threads <- snakemake@threads[[1]]
	register(MulticoreParam(workers=threads))

	# Snakemake variables
	deseq_dataset_obj <- snakemake@input[["deseq_dataset_r_obj"]]
	output_file_paths <- snakemake@params[["output_file_paths"]]

    # Load deseq dataset object
    dds <- readRDS(deseq_dataset_obj)

	# Run explorative analysis
	run_deseq2_explorative_analysis(dds, output_file_paths)
}


# ----------------- Run main function -----------------
main_function()
