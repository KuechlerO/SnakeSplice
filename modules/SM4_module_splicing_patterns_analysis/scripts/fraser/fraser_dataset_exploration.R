# load FRASER library
library("FRASER")

library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("org.Hs.eg.db")

# Requierements: 1. Sample annotation,
# 2. Two count matrices are needed: one containing counts for the splice junctions, i.e. the
# split read counts, and one containing the splice site counts, i.e. the counts of non
# split reads overlapping with the splice sites present in the splice junctions.


set_up_fraser_dataset_object <- function(sample_annotation_file_path) {
  #' Function to set up a FRASER object
  #'
  #' @param sample_annotation_file_path     Path to sample annotation file
  #' @param output_dir_path     Path to output directory
  #'
  #' @return FRASER object

  # Load annotation file
  annotationTable <- fread(sample_annotation_file_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  annotationTable$bamFile <- file.path(annotationTable$bamFile)   # Required for FRASER

  # --------------- Creating a FRASER object ----------------
  # create FRASER object
  settings <- FraserDataSet(colData=annotationTable, name="Fraser Dataset")

  # Via count reads
  fds <- countRNAData(settings)

  # Via raw counts
  # junctionCts <- fread(additional_junction_counts_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  # spliceSiteCts <- fread(additional_splice_site_counts_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  # fds <- FraserDataSet(colData=annotationTable, junctions=junctionCts, spliceSites=spliceSiteCts, workingDir="FRASER_output")

  return(fds)
}


run_filtering <- function(fraser_object,
                          plot_filter_expression_file, plot_cor_psi5_heatmap_file,
                          plot_cor_psi3_heatmap_file, plot_cor_theta_heatmap_file) {
  #' Function to run filtering
  #'
  #' @param fraser_object     FRASER object
  #' @param output_dir_path     Path to output directory
  #'
  #' @return FRASER object


  # --------------- Filtering ----------------
  # Compute main splicing metric -> The PSI-value
  fds <- calculatePSIValues(fraser_object)
  # Run filters on junctions: At least one sample has 20 reads, and at least 5% of the samples have at least 1 reads
  # Filter=FALSE, since we first plot and subsequently apply subsetting
  fds <- filterExpressionAndVariability(fds,
                                        minExpressionInOneSample=20,
                                        minDeltaPsi=0.0,  # Only junctions with a PSI-value difference of at least x% between two samples are considered
                                        filter=FALSE       # If TRUE, a subsetted fds containing only the introns that passed all filters is returned.
                                        )

  # Plot filtering results
  jpeg(plot_filter_expression_file, width=800, height=800)
  print(plotFilterExpression(fds, bins=100))
  dev.off()

  # Finally apply filter results
  fds_filtered <- fds[mcols(fds, type="j")[,"passed"],]

  # ---------------- Heatmaps of correlations ----------------
  # 1. Correlation of PSI5
  tryCatch(
    expr = {
      # Heatmap of the sample correlation
      jpeg(plot_cor_psi5_heatmap_file, width=800, height=800)
      plotCountCorHeatmap(fds_filtered, type="psi5", logit=TRUE, normalized=FALSE)
      dev.off()
    },
    error = function(e) {
        print("Error in creating Heatmap of the sample correlation")
        print(e)
    }
  )
  # tryCatch(
  #   expr = {
  #     # Heatmap of the intron/sample expression
  #     jpeg(plot_cor_psi5_top100_heatmap_file, width=800, height=800)
  #     plotCountCorHeatmap(fds_filtered, type="psi5", logit=TRUE, normalized=FALSE,
  #                     plotType="junctionSample", topJ=100, minDeltaPsi = 0.01)
  #     dev.off()
  #   },
  #   error = function(e) {
  #       print("Error in creating Heatmap of the intron/sample expression")
  #       print(e)
  #   }
  # )

  # 2. Correlation of PSI3
  tryCatch(
      expr = {
      # Heatmap of the sample correlation
      jpeg(plot_cor_psi3_heatmap_file, width=800, height=800)
      plotCountCorHeatmap(fds_filtered, type="psi3", logit=TRUE, normalized=FALSE)
      dev.off()
      },
      error = function(e) {
          print("Error in creating Heatmap of the sample correlation")
          print(e)
      }
  )
  # tryCatch(
  #   expr = {
  #     # Heatmap of the intron/sample expression
  #     jpeg(plot_cor_psi3_top100_heatmap_file, width=800, height=800)
  #     plotCountCorHeatmap(fds_filtered, type="psi3", logit=TRUE, normalized=FALSE,
  #                     plotType="junctionSample", topJ=100, minDeltaPsi = 0.01)
  #     dev.off()
  #   },
  #   error = function(e) {
  #       print("Error in creating Heatmap of the intron/sample expression")
  #       print(e)
  #   }
  # )

  # 3. Correlation of Theta
  tryCatch(
      expr = {
      # Heatmap of the sample correlation
      jpeg(plot_cor_theta_heatmap_file, width=800, height=800)
      plotCountCorHeatmap(fds_filtered, type="theta", logit=TRUE, normalized=FALSE)
      dev.off()
      },
      error = function(e) {
          print("Error in creating Heatmap of the sample correlation")
          print(e)
      }
  )
  # tryCatch(
  #   expr = {
  #     # Heatmap of the intron/sample expression
  #     jpeg(plot_cor_theta_top100_heatmap_file, width=800, height=800)
  #     plotCountCorHeatmap(fds_filtered, type="theta", logit=TRUE, normalized=FALSE,
  #                     plotType="junctionSample", topJ=100, minDeltaPsi = 0.01)
  #     dev.off()
  #   },
  #   error = function(e) {
  #       print("Error in creating Heatmap of the intron/sample expression")
  #       print(e)
  #   }
  # )

  return(fds_filtered)
}


detect_dif_splice <- function(fraser_object, output_fraser_analysis_set_object_file,
                              plot_normalized_cor_psi5_heatmap_file,
                              plot_normalized_cor_psi3_heatmap_file,
                              plot_normalized_cor_theta_heatmap_file) {
  #' Function to detect differential splicing
  #'
  #' @param fraser_object     FRASER object
  #' @param output_dir_path     Path to output directory
  #' @param summary_table_file     Path to summary table file
  #'
  #' @return FRASER object


  # ----------------- Detection of differential splicing -----------------
  # 1. Fitting the splicing model:
  # Normalizing data and correct for confounding effects by using a denoising autoencoder
  # This is computational heavy on real size datasets and can take awhile

  # q: The encoding dimension to be used during the fitting procedure. Can be fitted with optimHyperParams
  # see: https://rdrr.io/bioc/FRASER/man/optimHyperParams.html
  fds <- FRASER(fraser_object, q=c(psi5=3, psi3=5, theta=2))

  # Plot 1: PSI5
  tryCatch(
    expr = {
      # Check results in heatmap
      jpeg(plot_normalized_cor_psi5_heatmap_file, width=800, height=800)
      plotCountCorHeatmap(fds, type="psi5", normalized=TRUE, logit=TRUE)
      dev.off()
    },
    error = function(e) {
        print("Error in creating Heatmap of the sample correlation")
        print(e)
    }
  )

  # Plot 2: PSI3
  tryCatch(
      expr = {
      # Check results in heatmap
      jpeg(plot_normalized_cor_psi3_heatmap_file, width=800, height=800)
      plotCountCorHeatmap(fds, type="psi3", normalized=TRUE, logit=TRUE)
      dev.off()
      },
      error = function(e) {
          print("Error in creating Heatmap of the sample correlation")
          print(e)
      }
  )

  # Plot 3: Theta
  tryCatch(
      expr = {
      # Check results in heatmap
      jpeg(plot_normalized_cor_theta_heatmap_file, width=800, height=800)
      plotCountCorHeatmap(fds, type="theta", normalized=TRUE, logit=TRUE)
      dev.off()
      },
      error = function(e) {
          print("Error in creating Heatmap of the sample correlation")
          print(e)
      }
  )


  # 2. Differential splicing analysis
  # 2.1 annotate introns with the HGNC symbols of the corresponding gene
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  orgDb <- org.Hs.eg.db
  fds <- annotateRangesWithTxDb(fds, txdb=txdb, orgDb=orgDb)

  # 2.2 retrieve results with default and recommended cutoffs (padj <= 0.05 and |deltaPsi| >= 0.3)
  print("Saving FraserAnalysisDataSetTest results")
  # Saves RDS-file into savedObjects folder
  saveFraserDataSet(fds, dir=dirname(dirname(output_fraser_analysis_set_object_file)),
                    name=basename(output_fraser_analysis_set_object_file))


  # ----------------- Finding splicing candidates in patients -----------------
  # -> Plotting the results
  # tryCatch(
  #   expr = {
      # -------- Sample specific plots --------
      # jpeg(file.path(output_dir_path, "psi5_volcano_plot_sample1.jpg"), width=800, height=800)
      # plotVolcano(fds, type="psi5", annotationTable$sampleID[1])
      # dev.off()

      # jpeg(file.path(output_dir_path, "psi5_expression_sample1.jpg"), width=800, height=800)
      # plotExpression(fds, type="psi5", result=sampleRes[1])
      # dev.off()

      # jpeg(file.path(output_dir_path, "expected_vs_observed_psi_sample1.jpg"), width=800, height=800)
      # plotExpectedVsObservedPsi(fds, result=sampleRes[1])
      # dev.off()
  #   },
  #   error = function(e) {
  #       print("Error in creating plots")
  #       print(e)
  #   }
  # )

  return(fds)
  }


main_function <- function() {
  in_sample_annotation_file <- snakemake@input[["sample_annotation_file"]]

  # Output: Plot files - After filtering, no normalization
  plot_filter_expression_file <- snakemake@output[["plot_filter_expression_file"]]
  plot_cor_psi5_heatmap_file <- snakemake@output[["plot_cor_psi5_heatmap_file"]]
  plot_cor_psi3_heatmap_file <- snakemake@output[["plot_cor_psi3_heatmap_file"]]
  plot_cor_theta_heatmap_file <- snakemake@output[["plot_cor_theta_heatmap_file"]]

  # ToDO: Set plotType to "sampleCorrelation", however this plots are not helpful and can be ignored...
  # plot_cor_psi5_top100_heatmap_file <- snakemake@output[["plot_cor_psi5_top100_heatmap_file"]]
  # plot_cor_psi3_top100_heatmap_file <- snakemake@output[["plot_cor_psi3_top100_heatmap_file"]]
  # plot_cor_theta_top100_heatmap_file <- snakemake@output[["plot_cor_theta_top100_heatmap_file"]]

  # Output: Plot files - After filtering, normalization
  plot_normalized_cor_psi5_heatmap_file <- snakemake@output[["plot_normalized_cor_psi5_heatmap_file"]]
  plot_normalized_cor_psi3_heatmap_file <- snakemake@output[["plot_normalized_cor_psi3_heatmap_file"]]
  plot_normalized_cor_theta_heatmap_file <- snakemake@output[["plot_normalized_cor_theta_heatmap_file"]]

  # Output: Differential splicing analysis
  output_fraser_dataset_object_file <- snakemake@output[["fraser_data_set_object_file"]]


  # TODO: Integrate additional count files from external resources -> Failed...
  # additional_junction_counts_file <- snakemake@params[["additional_junction_counts_file"]]
  # additional_splice_site_counts_file <- snakemake@params[["additional_splice_site_counts_file"]]

  threads <- snakemake@threads
  register(MulticoreParam(workers=threads))

  # 1. Create FRASER object
  fraser_obj <- set_up_fraser_dataset_object(in_sample_annotation_file)
  print("FRASER: FRASER dataset object created")

  # 2. Run filtering
  filtered_fraser_obj <- run_filtering(fraser_obj,
                                       plot_filter_expression_file,
                                       plot_cor_psi5_heatmap_file,
                                       plot_cor_psi3_heatmap_file,
                                       plot_cor_theta_heatmap_file)
  print("FRASER: Filtering done")

  # 3. Detect differential splicing
  detect_dif_splice(filtered_fraser_obj, output_fraser_dataset_object_file,
                    plot_normalized_cor_psi5_heatmap_file,
                    plot_normalized_cor_psi3_heatmap_file,
                    plot_normalized_cor_theta_heatmap_file
                    )
  print("FRASER: Differential splicing analysis done")
}

main_function()
