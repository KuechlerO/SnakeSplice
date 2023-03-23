# Cummerbund manual: 
# http://bioconductor.org/packages/2.11/bioc/vignettes/cummeRbund/inst/doc/cummeRbund-manual.pdf

# load package
library("cummeRbund")


global_statistics_and_qc <- function(cuff_obj, output_dir) {
  # ---- Global Statistics and Quality Control -----------

  # Dispersion explained:
  # https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8
  # https://www.biostars.org/p/167688/
  # https://support.bioconductor.org/p/75260/

  # ------------- 1. Dispersion ----------------
  # -> visualizes the estimated overdispersion for each sample
  # uses cufflinks emitted data (mean counts, variance, & dispersion)
  # -> http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/
  genes.disp<-dispersionPlot(genes(cuff_obj))
  pdf(file=file.path(output_dir, "cummerbund_figures/dispersion_genes.pdf"))
  plot(genes.disp)				# Plot is displayed
  dev.off()

  cuff_obj.disp<-dispersionPlot(cuff_obj)
  pdf(file=file.path(output_dir, "dispersion_cuff.pdf"))
  plot(cuff_obj.disp)				# Plot is displayed
  dev.off()

  # ------ 2. Distributions of FPKM scores across samples ----------
  # 2.1.) csDensity plots
  dens<-csDensity(genes(cuff_obj))
  pdf(file=file.path(output_dir, "fpkm_density_genes.pdf"))
  plot(dens)
  dev.off()

  dens<-csDensity(isoforms(cuff_obj))
  pdf(file=file.path(output_dir, "fpkm_density_isoforms.pdf"))
  plot(dens)
  dev.off()

  # 2.2.) Boxplots
  b<-csBoxplot(genes(cuff_obj))
  pdf(file=file.path(output_dir, "fpkm_density_boxplot_genes.pdf"))
  plot(b)
  dev.off()

  b<-csBoxplot(isoforms(cuff_obj))
  pdf(file=file.path(output_dir, "fpkm_density_boxplot_isoforms.pdf"))
  plot(b)
  dev.off()

  # 2.3.) Matrix of pairwise scatterplots
  s<-csScatterMatrix(genes(cuff_obj))
  pdf(file=file.path(output_dir, "fpkm_density_matrix_genes.pdf"))
  plot(s)
  dev.off()

  s<-csScatterMatrix(isoforms(cuff_obj))
  pdf(file=file.path(output_dir, "fpkm_density_matrix_isoforms.pdf"))
  plot(s)
  dev.off()

  # 2.4.) Volcano plots -> Explore relationship between fold-change and significance
  v<-csVolcanoMatrix(genes(cuff_obj))
  pdf(file=file.path(output_dir, "cs_volcano_matrix_genes.pdf"))
  plot(v)
  dev.off()

  v<-csVolcanoMatrix(isoforms(cuff_obj))
  pdf(file=file.path(output_dir, "cs_volcano_matrix_isoforms.pdf"))
  plot(v)
  dev.off()
}


analyse_differential_expression <- function(cuff_obj, output_dir) {
  # ------ Differential expression -------------
  # 1.) genes
  # all
  gene.diff <- diffData(genes(cuff_obj))
  write.csv(gene.diff, file.path(output_dir, "gene_diff.csv"), row.names=F)

  # only significant -> with gene names
  sig_gene_ids <- getSig(cuff_obj,level="genes",alpha=0.05)
  if (NROW(sig_gene_ids) > 0) {
    sigFeatures <- getFeatures(cuff_obj,sig_gene_ids,level="genes")
    sigData <- diffData(sigFeatures)
    sigData <- subset(sigData, (significant == 'yes'))
    names <- featureNames(sigFeatures)
    sigOutput <- merge(names, sigData, by.x="tracking_id", by.y="gene_id")

    # Patch the merged table to have the original name for the ID column.
    # This is always the first column for the examples we've seen.
    colnames(sigOutput)[1] <- "gene_id"
    write.table(sigOutput, file.path(output_dir, "gene_diff_only_significant.tsv"), sep='\t', row.names = F,
                col.names = T, quote = F)
  } else {
    sink(file.path(output_dir, "gene_diff_only_significant.tsv"))
    cat("No significantly differently expressed genes detected")
    sink()
  }

  # 2.) isoforms
  # all
  isoform.diff <- diffData(isoforms(cuff_obj))
  write.csv(isoform.diff, file.path(output_dir, "isoform_diff.csv"), row.names=F)

  # only significant -> with gene names
  sig_isoforms_ids <- getSig(cuff_obj, level="isoforms", alpha=0.05)
  if (NROW(sig_isoforms_ids) > 0) {
    sigFeatures <- getFeatures(cuff_obj, sig_isoforms_ids, level="isoforms")
    sigData <- diffData(sigFeatures)
    sigData <- subset(sigData, (significant == 'yes'))
    names <- featureNames(sigFeatures)
    sigOutput <- merge(names, sigData, by.x="tracking_id", by.y="isoform_id")

    # Patch the merged table to have the original name for the ID column.
    # This is always the first column for the examples we've seen.
    colnames(sigOutput)[1] <- "gene_id"
    write.table(sigOutput, file.path(output_dir, "isoform_diff_only_significant.tsv"), sep='\t', row.names = F,
                col.names = T, quote = F)
  } else {
    sink(file.path(output_dir, "isoform_diff_only_significant.tsv"))
    cat("No significantly differently expressed isoforms detected")
    sink()
  }
}


create_count_matrices <- function(cuff_obj, output_dir) {
  " Create count matrices for genes and isoforms "

  # ------------- CSV files ----------------
  # access feature lvl data
  gene.features <- annotation(genes(cuff_obj))
  write.csv(gene.features, file.path(output_dir, "gene_features.csv"), row.names = FALSE)

  gene.fpkm <- fpkm(genes(cuff_obj))
  write.csv(gene.fpkm, file.path(output_dir, "gene_fpkm.csv"), row.names = FALSE)

  # raw and normalized (on sequencing depth?) fragment counts
  gene.counts <- count(genes(cuff_obj))
  write.csv(gene.counts, file.path(output_dir, "gene_counts.csv"), row.names = FALSE)

  # -- isoforms
  isoform.features <- annotation(isoforms(cuff_obj))
  write.csv(isoform.features, file.path(output_dir, "isoform_features.csv"), row.names = FALSE)

  isoform.fpkm <- fpkm(isoforms(cuff_obj))
  write.csv(isoform.fpkm, file.path(output_dir, "isoform_fpkm.csv"), row.names = FALSE)

  isoform.counts <- count(isoforms(cuff_obj))
  write.csv(isoform.counts, file.path(output_dir, "isoform_counts.csv"), row.names = FALSE)


  # ----------- create PDFs -----------
  # FPKM matrices
  gene.fpkm.matrix<-fpkmMatrix(genes(cuff_obj))
  pdf(file=file.path(output_dir, "fpkm_matrix_genes.pdf"))
  plot(gene.fpkm.matrix)
  dev.off()

  isoform.fpkm.matrix<-fpkmMatrix(isoforms(cuff_obj))
  pdf(file=file.path(output_dir, "fpkm_matrix_isoforms.pdf"))
  plot(isoform.fpkm.matrix)
  dev.off()

  # Count matrices
  gene.count.matrix<-countMatrix(genes(cuff_obj))
  pdf(file=file.path(output_dir, "count_matrix_genes.pdf"))
  plot(gene.count.matrix)
  dev.off()

  isoforms.count.matrix<-countMatrix(isoforms(cuff_obj))
  pdf(file=file.path(output_dir, "count_matrix_isoforms.pdf"))
  plot(isoforms.count.matrix)
  dev.off()
}


single_gene_analysis <- function(cuff_obj, gene_of_interest_id, output_dir) {
  " Detailed analysis of a single gene of interest"
  # IV) ---- Single gene -----------
  myGeneId <- gene_of_interest_id
  myGene<-getGene(cuff_obj,myGeneId)
  header <- paste("Single Gene", arg_gene, ":", sep=" ")
  capture.output(myGene, file=file.path(output_dir, "analysis_output.txt"), append = TRUE)

  output_file_genes <- file.path(output_dir, paste(arg_gene, "_gene_fpkm.csv", sep=""))
  write.csv(fpkm(myGene), output_file_genes, row.names = FALSE)

  output_file_isoforms <- file.path(output_dir, paste(arg_gene, "_isoforms_fpkm.csv", sep=""))
  write.csv(fpkm(isoforms(myGene)), output_file_isoforms, row.names = FALSE)

  # Plots
  gl <- expressionPlot(myGene)
  output_file_1 <- file.path(output_dir, paste("expressionPlot_singleGene_", arg_gene, ".pdf", sep=""))
  pdf(file=output_file_1)
  plot(gl)
  dev.off()

  # Expression plot of all isoforms of a single gene with FPKMs exposed
  gl.iso.rep <- expressionPlot(isoforms(myGene))
  output_file_2 <- file.path(output_dir, paste("expressionPlot_isoforms_singleGene_", arg_gene, ".pdf", sep=""))
  pdf(file=output_file_2)
  plot(gl.iso.rep)
  dev.off()

  # Expression plot of all CDS for a single gene with FPKMS exposed
  gl.cds.rep<-expressionPlot(CDS(myGene))
  output_file_3 <- file.path(output_dir, paste("expressionPlot_cds_singleGene_", arg_gene, ".pdf", sep=""))
  pdf(file=output_file_3)
  plot(gl.cds.rep)
  dev.off()

  # Detailed feature graph
  trackList<-list()
  myStart<-min(features(myGene)$start)
  myEnd<-max(features(myGene)$end)
  myChr<-unique(features(myGene)$seqnames)
  genome<-arg_genome
  ideoTrack <- IdeogramTrack(genome = genome, chromosome = myChr)
  trackList<-c(trackList,ideoTrack)   # appending ideoTrack -> chromosome
  axtrack<-GenomeAxisTrack()
  trackList<-c(trackList,axtrack)     # appending axtrack -> genome-Axis
  genetrack<-makeGeneRegionTrack(myGene)

  trackList<-c(trackList,genetrack)   # appending genetrack -> the mapping results
  biomTrack<-BiomartGeneRegionTrack(genome=genome,chromosome=as.character(myChr), start=myStart,end=myEnd,name="ENSEMBL",showId=T)
  trackList<-c(trackList,biomTrack)   # Biomart transcripts
  conservation <- UcscTrack(genome = genome, chromosome = myChr, track = "Conservation", table = "multiz100way",from = myStart-2000, to = myEnd+2000, trackType = "DataTrack",start = "start", end = "end", data = "score",type = "hist", window = "auto", col.histogram = "darkblue",fill.histogram = "darkblue", ylim = c(-3.7, 4),name = "Conservation")
  trackList<-c(trackList,conservation)    # conservation
  # Plot detailed graph: Chromosome on top...
  pdf(file.path(output_dir, "detailed_track.pdf"))
  plotTracks(trackList,from=myStart-2000,to=myEnd+2000)
  dev.off()
}

# ----------- Main -----------
main <- function() {
  # Inputs
  cufflinks_merged_transcriptome_assemblies_gtf <- snakemake@input["merged_cufflinks_transcriptome_assemblies_gtf"]

  # Params
  cufflinks_output_files_dir <- snakemake@params["cufflinks_output_files_dir"]
  genome_build <- snakemake@params["original_genome_build"]
  chosen_genes_of_interest <- snakemake@params[["chosen_genes_of_interest"]]

  # Output directory
  cummerbund_output_dir <- snakemkae@output[["cummerbund_output_dir"]]
  cummerbund_summary_results <- snakemake@output[["cummerbund_summary_results"]]

  # read results -> create SQLite db
  # Rebuild is important: Create always new database
  cuff_obj <- readCufflinks(dir=cufflinks_output_files_dir, gtfFile=cufflinks_merged_transcriptome_assemblies_gtf,
                            genome=genome_build, rebuild=T)
  capture.output(cuff_obj, file=file.path(cummerbund_output_dir, cummerbund_summary_results), append = FALSE)

  # Global statistics and qc
  global_statistics_and_qc(cuff_obj, cummerbund_output_dir)

  # Count matrices
  create_count_matrices(cuff_obj, cummerbund_output_dir)

  # Single gene analysis
  for (gene_of_interest in chosen_genes_of_interest) {
      single_gene_analysis(cuff_obj, gene_of_interest, cummerbund_output_dir)
  }
}

# Execute main
main()
