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


explore_outrider_results <- function(ods, output_dir,
                                     sample_ids, genes_of_interest,
                                     p005_file, p010_file,
                                     ab_genes_per_sample_file, ab_samples_per_gene_file) {
    ### Extract OUTRIDER results
    # Input: ods: OUTRIDER data set
    #        sample_ids: Sample IDs
    #        genes_of_interest: Genes of interest
    #        plot_output_dir: Output directory for plots
    #        outrider_obj_file: File path to save OUTRIDER object
    #        p005_file: File path to save p-value 5% results
    #        p010_file: File path to save p-value 10% results
    #        ab_genes_per_sample_file: File path to save aberrant genes per sample
    #        ab_samples_per_gene_file: File path to save aberrant samples per gene
    # Output: OUTRIDER object, p005 results, p010 results, aberrant genes per sample, aberrant samples per gene

    # -------- 5. Output Results ------------
    # -------- 5.1 significant results ----------
    # get results (default only significant, padj < 0.05)
    res <- results(ods)
    res <- res[order(res$padjust),]
    write.table(res, file=p005_file, sep="\t", quote=FALSE, row.names=FALSE)

    # get results (default only significant, padj < 0.10)
    res <- results(ods, padjCutoff=0.1)
    res <- res[order(res$padjust),]
    write.table(res, file=p010_file, sep="\t", quote=FALSE, row.names=FALSE)

    # -------- 5.2 Aberrant expression ----------
    # number of aberrant genes per sample
    nr_aberrant_genes_per_sample <- sort(aberrant(ods, by="sample"))
    # Use sink to include sample-names in output file
    sink(ab_genes_per_sample_file); print(nr_aberrant_genes_per_sample); sink()

    # number of aberrant samples per gene
    nr_aberrant_samples_per_gene <- sort(aberrant(ods, by="gene"))
    sink(ab_samples_per_gene_file); print(nr_aberrant_samples_per_gene); sink()


    # -------- 5.3 Volcano Plots for p-values ----------
    print("Given sample IDs:")
    print(sample_ids)

    # Convert sample IDs
    # 1. "-" are converted into "."
    # 2. Samples which start with numbers get the prefix: "X"
    for (i in (1:length(sample_ids))) {
        sample_ids[i] <- gsub("-", ".", sample_ids[i])
        # if ( grepl("^[1-9].*$", sample_ids[i]) ) {
        #     sample_ids[i] <- paste0("X", sample_ids[i])  # Add X to sample IDs starting with a number -> Comply with Outrider conversion
        # }
    }

    tryCatch(
        expr = {
            for (sample_id in sample_ids) {
                html_file_path <- file.path(output_dir, paste(sample_id, "_pvalues_volcano.html", sep=""))
                png_file_path <- file.path(output_dir, paste(sample_id, "_pvalues_volcano.png", sep=""))

                # A. Create interactive Plotly plot
                interactive_plot <- plotVolcano(ods, sample_id, basePlot=FALSE)
                htmlwidgets::saveWidget(as_widget(interactive_plot), html_file_path)

                # B. Create static plot
                png(png_file_path)
                print(plotVolcano(ods, sample_id, basePlot=TRUE))
                dev.off()
            }
        },
        error = function(e) {
            print("Error in creating volcano plots")
            print(e)
        },
        finally = {
            print("Yes, Volcano plots are done!")
        }
    )

    # -------- 5.4 Gene level plots ----------
    # 5.4.1 Expression Rank
    tryCatch(
        expr = {
            for (gene_name in genes_of_interest) {
                html_file_path <- file.path(output_dir, paste(gene_name, "_expressionRank.html"))
                jpeg_file_path <- file.path(output_dir, paste(gene_name, "_expressionRank.jpg"))

                # A. Create interactive Plotly plot
                interactive_plot <- plotExpressionRank(ods, gene_name, basePlot=FALSE)
                htmlwidgets::saveWidget(as_widget(interactive_plot), html_file_path)

                # B. Create static plot
                jpeg(file=jpeg_file_path)
                print(plotExpressionRank(ods, gene_name, basePlot=TRUE))
                dev.off()
            }
        },
        error = function(e) {
            print("Error in plotExpressionRank")
            print(e)
        },
        finally = {
            print("Yes, Expression Rank plots are done!")
        }
    )

    # 5.4.2 Quantile-Quantile-Plots (Q-Q plots)
    tryCatch(
        expr = {
            for (gene_name in genes_of_interest) {
                jpeg_file_path <- file.path(output_dir, paste(gene_name, "_qqPlot.jpg", sep=""))

                # B. Create static plot
                jpeg(file=jpeg_file_path)
                print(plotQQ(ods, gene_name))
                dev.off()
            }
        },
        error = function(e) {
            print("Error in plotQQ")
            print(e)
        },
        finally = {
            print("Yes, Q-Q plots are done!")
        }
    )

    # 5.4.3 Observed versus expected Expression
    tryCatch(
        expr = {
            for (gene_name in genes_of_interest) {
                html_file_path <- file.path(output_dir, paste(gene_name, "_expectedVsObservedExpression.html", sep=""))
                jpeg_file_path <- file.path(output_dir, paste(gene_name, "_expectedVsObservedExpression.jpg", sep=""))

                # A. Create interactive Plotly plot
                interactive_plot <- plotExpectedVsObservedCounts(ods, gene_name, basePlot=FALSE)
                htmlwidgets::saveWidget(as_widget(interactive_plot), html_file_path)

                # B. Create static plot
                jpeg(file=jpeg_file_path)
                print(plotExpectedVsObservedCounts(ods, gene_name, basePlot=TRUE))
                dev.off()
            }
        },
        error = function(e) {
            print("Error in plotExpectedVsObservedCounts")
            print(e)
        },
        finally = {
            print("Yes, Observed versus expected Expression plots are done!")
        }
    )
}


# TODO: clean up naming of output files
# TODO: Do not give output dir, but directly the output files?!
main_function <- function() {
  # Input
  input_final_outrider_obj_file <- snakemake@input[["outrider_object_file"]]

  # Output directory -> For saving plots
  output_dir <- snakemake@output[[1]]

  # Outputs
  significant_results_p005_output_file <- snakemake@output[["significant_results_p005_file"]]
  significant_results_p010_output_file <- snakemake@output[["significant_results_p010_file"]]
  nr_aberrant_genes_per_sample_output_file <- snakemake@output[["nr_aberrant_genes_per_sample"]]
  nr_aberrant_samples_per_gene_output_file <- snakemake@output[["nr_aberrant_samples_per_gene"]]

  # Params
  sample_ids <- snakemake@params[["sample_ids"]]
  genes_of_interest <- snakemake@params[["genes_of_interest"]]

  # ------- 1. Create output directory ---------
  dir.create(file.path(output_dir), recursive=TRUE, showWarnings=FALSE)     # Create plot-directory

  ods <- readRDS(input_final_outrider_obj_file)

  # -------- 5. Output Results ------------
  explore_outrider_results(ods, output_dir,
                           sample_ids, genes_of_interest,
                           significant_results_p005_output_file, significant_results_p010_output_file,
                           nr_aberrant_genes_per_sample_output_file, nr_aberrant_samples_per_gene_output_file)
}

main_function()
