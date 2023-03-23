library("DESeq2")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("clusterProfiler")
library("ggplot2")
library("enrichplot")


create_gsea_plots <- function(gsea_obj, dotplot_file_path, gsea_plot_file_path_1, gsea_plot_file_path_2, method) {
  # Create plots for GSEA results

  # -------- Plottings ---------
  # 1. Dotplot
  jpeg(dotplot_file_path,width=800, height=800)
  print(clusterProfiler::dotplot(gsea_obj, showCategory=30)
          + ggplot2::ggtitle(paste0("DotPlot for GSE-analysis (top 30 results) with method: ", method)))
  dev.off()

  # 2. GSEA-Plot for top 10 results
  jpeg(gsea_plot_file_path_1, width=800, height=800)
  print(enrichplot::gseaplot2(gsea_obj, geneSetID = 1:5, pvalue_table=FALSE,
                             title = paste0("GSEA-Plot for top 1-5 results with method: ", method)))
  dev.off()

  jpeg(gsea_plot_file_path_2, width=800, height=800)
  print(enrichplot::gseaplot2(gsea_obj, geneSetID = 6:10, pvalue_table=FALSE,
                             title = paste0("GSEA-Plot for top 6-10 results with method: ", method)))
  dev.off()

  # for (count in c(1:10)) {
  #   jpeg(paste0(plot_output_dir, "/gsea_plot_", method, "_", count, ".jpg"), width=800, height=800)
  #   print(clusterProfiler::gseaplot(gsea_obj, geneSetID=1, pvalue_table=TRUE))
  #   dev.off()
  # }
}

explore_gsea_go <- function(ordered_gene_list, summary_file_path, gsea_obj_file_path) {
  # Gene Set Enrichtment Analyses
  # params:
  #   ordered_gene_list: Ordered (i.e. deseq2 stat) list of genes

  # ----- 1. GO Enrichment --------
  # GO comprises three orthogonal ontologies, i.e. molecular function (MF),
  # biological process (BP), and cellular component (CC)
  go_gsea <- clusterProfiler::gseGO(ordered_gene_list,
               ont = "BP",
               keyType = "ENSEMBL",
               OrgDb = "org.Hs.eg.db",
               verbose = TRUE)

  df_go_gsea <- as.data.frame(go_gsea)
  df_go_gsea <- df_go_gsea[order(df_go_gsea$p.adjust),]
  write.csv(df_go_gsea, file=summary_file_path)

  # Save GSEA object
  saveRDS(go_gsea, file=gsea_obj_file_path)

  return(go_gsea)
}


explore_gsea_kegg <- function(ordered_gene_list, summary_file_path, gsea_obj_file_path) {
  # KEGG pathway enrichment analysis
  # params:
  #   ordered_gene_list: Ordered (i.e. deseq2 stat) list of genes

  names(ordered_gene_list) <- mapIds(org.Hs.eg.db, keys=names(ordered_gene_list), column="ENTREZID",
                                                      keytype="ENSEMBL", multiVals="first")
  # res$symbol <- mapIds(org.Hs.eg.db, keys=row.names(res), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
  # res$entrez <- mapIds(org.Hs.eg.db, keys=row.names(res), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
  # res$name <- mapIds(org.Hs.eg.db, keys=row.names(res), column="GENENAME", keytype="ENSEMBL", multiVals="first")

  # ----- 1. GO Enrichment --------
  # GO comprises three orthogonal ontologies, i.e. molecular function (MF),
  # biological process (BP), and cellular component (CC)
  kegg_gsea <- clusterProfiler::gseKEGG(geneList=ordered_gene_list,
                                        organism='hsa',
                                        verbose=TRUE)

  df_kegg_gsea <- as.data.frame(kegg_gsea)
  df_kegg_gsea <- df_kegg_gsea[order(df_kegg_gsea$p.adjust),]
  write.csv(df_kegg_gsea, file=summary_file_path)

  # Save GSEA object
  saveRDS(kegg_gsea, file=gsea_obj_file_path)
  return(kegg_gsea)
}


explore_gsea_wp <- function(ordered_gene_list,
                            summary_file_path, gsea_obj_file_path) {
  # WikiPathway
  # params:
  #   ordered_gene_list: Ordered (i.e. deseq2 stat) list of genes

  names(ordered_gene_list) <- mapIds(org.Hs.eg.db, keys=names(ordered_gene_list), column="ENTREZID",
                                                      keytype="ENSEMBL", multiVals="first")

  # ----- 1. GO Enrichment --------
  # GO comprises three orthogonal ontologies, i.e. molecular function (MF),
  # biological process (BP), and cellular component (CC)
  wp_gsea <- clusterProfiler::gseWP(
    geneList=ordered_gene_list,
    organism="Homo sapiens",
    verbose=TRUE)

  df_wp_gsea <- as.data.frame(wp_gsea)
  df_wp_gsea <- df_wp_gsea[order(df_wp_gsea$p.adjust),]
  write.csv(df_wp_gsea, file=summary_file_path)

  # Save GSEA object
  saveRDS(wp_gsea, file=gsea_obj_file_path)
  return(wp_gsea)
}


main <- function() {
  # Input
  input_dseq_dataset_obj <- snakemake@input$deseq_dataset_obj

  # Outputs
  gsea_go_obj_file_path <- snakemake@output$gsea_go_obj_file_path
  gsea_go_summary_file_path <- snakemake@output$gsea_go_summary_file_path
  gsea_kegg_obj_file_path <- snakemake@output$gsea_kegg_obj_file_path
  gsea_kegg_summary_file_path <- snakemake@output$gsea_kegg_summary_file_path
  gsea_wp_obj_file_path <- snakemake@output$gsea_wp_obj_file_path
  gsea_wp_summary_file_path <- snakemake@output$gsea_wp_summary_file_path

  # DotPlots
  dotplot_gsea_go_file_path <- snakemake@output$dotplot_gsea_go_file_path
  dotplot_gsea_kegg_file_path <- snakemake@output$dotplot_gsea_kegg_file_path
  dotplot_gsea_wp_file_path <- snakemake@output$dotplot_gsea_wp_file_path
  # GSEA Plots
  gsea_go_top10_plot_file_path_1 <- snakemake@output$gsea_go_top10_plot_file_path_1
  gsea_kegg_top10_plot_file_path_1 <- snakemake@output$gsea_kegg_top10_plot_file_path_1
  gsea_wp_top10_plot_file_path_1 <- snakemake@output$gsea_wp_top10_plot_file_path_1
  gsea_go_top10_plot_file_path_2 <- snakemake@output$gsea_go_top10_plot_file_path_2
  gsea_kegg_top10_plot_file_path_2 <- snakemake@output$gsea_kegg_top10_plot_file_path_2
  gsea_wp_top10_plot_file_path_2 <- snakemake@output$gsea_wp_top10_plot_file_path_2

  # Params
  input_algorithm <- snakemake@params$input_algorithm

  # Load DataSet
  dds <- readRDS(input_dseq_dataset_obj)
  # Create DESeq2 object
  dds <- DESeq(dds)
  res <- DESeq2::results(dds)

  # Filtering
  res <- na.omit(res)
  res <- res[res$baseMean >50,] # Filter out genes with low expression

  # Order output -> We choose stat, which takes log-Fold as well as SE into account
  # Alternative: lfc * -log10(P-value)
  # order descending so use minus sign
  res <- res[order(-res$stat),]

  # --------- Create input gene list ---------------
  # Extract stat values
  gene_list <- res$stat
  # Add rownames
  if (input_algorithm == "salmon") {
    # Ensembl-Transcript-IDs at first place
    names(gene_list) <- substr(rownames(res), 1, 15)
  }
  else if (input_algorithm == "kallisto") {
    # Ensembl-Transcript-IDs at second place (delimeter: "|")
    gene_ids_in_rows <- sapply(rownames(res), function(x) strsplit(x, '\\|')[[1]], USE.NAMES=FALSE)[2,]
    gene_ids_in_rows <- sapply(gene_ids_in_rows, function(x) substr(x, 1, 15), USE.NAMES=FALSE)
    names(gene_list) <- gene_ids_in_rows
  }
  else {
    stop("Unknown algorithm used for quantification")
  }

  # =========== Run  GSEA ===========
  # ----- 1. GO Enrichment --------
  go_gsea_obj <- explore_gsea_go(gene_list, gsea_go_summary_file_path, gsea_go_obj_file_path)
  create_gsea_plots(go_gsea_obj, dotplot_gsea_go_file_path,
                    gsea_go_top10_plot_file_path_1, gsea_go_top10_plot_file_path_2, "go")

  # ----- 2. KEGG Enrichment --------
  kegg_gsea_obj <- explore_gsea_kegg(gene_list, gsea_kegg_summary_file_path, gsea_kegg_obj_file_path)
  create_gsea_plots(kegg_gsea_obj, dotplot_gsea_kegg_file_path,
                    gsea_kegg_top10_plot_file_path_1, gsea_kegg_top10_plot_file_path_2, "kegg")

  # ----- 3. WikiPathway Enrichment --------
  wp_gsea_obj <- explore_gsea_wp(gene_list, gsea_wp_summary_file_path, gsea_wp_obj_file_path)
  create_gsea_plots(wp_gsea_obj, dotplot_gsea_wp_file_path,
                    gsea_wp_top10_plot_file_path_1, gsea_wp_top10_plot_file_path_2, "wp")
}

# Run main function
main()