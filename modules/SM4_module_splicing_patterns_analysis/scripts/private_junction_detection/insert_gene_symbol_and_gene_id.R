library("AnnotationDbi")
library("org.Hs.eg.db")


extract_gene_id_from_info_col <- function(data_frame_obj, info_col, gene_id_col="gene_ensembl_id") {
	"
	Extracts gene ID from info column.
	"
	# Extract gene-IDs from info_col
	# Each entry in info_col looks like this:
	# gene_id "ENSG00000186092"; transcript_id "ENST00000335137"; exon_number "1"; gene_name "OR4F5"; gene_biotype "protein_coding"; transcript_name "OR4F5-201"; exon_id "ENSE00002234944";
	# Extract the first part of the string, i.e. the gene_id
	gene_ids <- lapply(data_frame_obj[info_col], FUN=function(x) {
		gene_id <- gsub(pattern=".*gene_id \"", replacement="", x=x)
		gene_id <- gsub(pattern="\";.*", replacement="", x=gene_id)
		return(gene_id)
		}
	)
	data_frame_obj[gene_id_col] <- gene_ids

	return(data_frame_obj)
}


add_gene_symbol_and_entrez_id_to_results <- function(data_frame_obj,
														gene_ensembl_id_col="gene_ensembl_id",
														gene_name_col="gene_name") {
	"
	Adds gene symbols and entrez-IDs to results object.
	"
	gene_ids_vector <- as.vector(t(data_frame_obj[gene_ensembl_id_col]))

	# If empty gene_ids_vector, then return fill with NA
	if (length(gene_ids_vector) == 0) {
		data_frame_obj[gene_name_col] <- character(0)
	}

	else {
		# Add gene symbols
		# Something breaks here when setting a new column name
		data_frame_obj[gene_name_col] <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
															  keys=gene_ids_vector,
															  column="SYMBOL",
															  keytype="ENSEMBL",
															  multiVals="first")
	}
	return(data_frame_obj)
}




# Main function
main <- function() {
	# Input
	input_table_files <- snakemake@input
	# Output
	output_files <- snakemake@output

	# info_col
	info_col_name <- snakemake@params[["info_col_name"]]
	gene_ensembl_id_col_name = snakemake@params[["gene_ensembl_id_col_name"]]
	gene_name_col_name = snakemake@params[["gene_name_col_name"]]


	# Loop over input files
	for (i in seq_along(input_table_files)) {
		# Read input table
		df <- read.table(toString(input_table_files[i]), sep="\t", header=TRUE, stringsAsFactors=FALSE)

		# Extract gene ID from info column
		df <- extract_gene_id_from_info_col(df, info_col=info_col_name, gene_id_col=gene_ensembl_id_col_name)

		# Add gene symbols and entrez-IDs
		df <- add_gene_symbol_and_entrez_id_to_results(df,
			gene_ensembl_id_col=gene_ensembl_id_col_name, gene_name_col=gene_name_col_name)


		# Put gene_ensembl_id_col and gene_name_col to the front
		input_table <- df[, c(gene_ensembl_id_col_name, gene_name_col_name,
			setdiff(colnames(df), c(gene_ensembl_id_col_name, gene_name_col_name)))]

		# Write output table
		write.table(input_table, file=toString(output_files[i]), sep="\t", quote=FALSE, row.names=FALSE)
	}
}


# Run main function
main()
