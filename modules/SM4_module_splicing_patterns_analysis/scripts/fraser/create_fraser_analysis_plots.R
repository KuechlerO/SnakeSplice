# load FRASER library
library("FRASER")


main_function <- function() {
  input_fraser_analysis_set_object_file <- snakemake@input[["fraser_analysis_set_object_file"]]

  # Output: Differential splicing analysis - Plots
  output_summary_table_file <- toString(snakemake@output[["csv_summary_table_file"]])
  plot_aberrant_events_per_sample_file <- toString(snakemake@output["plot_aberrant_events_per_sample_file"][1])
  plot_qq_plot_file <- toString(snakemake@output["plot_qq_plot_file"][1])

  # 1. Create FRASER object
  dir_name <- dirname(dirname(input_fraser_analysis_set_object_file))
  file_name <- basename(input_fraser_analysis_set_object_file)
  fds <- FRASER::loadFraserDataSet(dir=dir_name, name=file_name)
  print("FRASER: FRASER dataset object loaded")


  # 2. Collect results and save them in a data frame
  res <- as.data.table(results(fds))
  resOrdered <- res[order(res$pValue),]					# Sort results by p-value
  # Exporting results
  resOrderedDF <- as.data.frame(resOrdered)
  write.csv(resOrderedDF, file=output_summary_table_file)


  # 3. Create Plots
  # 3.1 Plot the number of aberrant events per sample
  tryCatch(
    expr = {
      # Plot number of aberrant events per sample based on the given cutoff values
      print("Plotting number of aberrant events per sample")
      print(plot_aberrant_events_per_sample_file)
      png(filename=plot_aberrant_events_per_sample_file, width=800, height=800)
      print(FRASER::plotAberrantPerSample(fds))
      dev.off()
    },
    error = function(e) {
      print("Error in creating aberrant events per sample plot")
      print(e)
    }
  )

  # 3.2 Plot the qq-plot
  tryCatch(
    expr = {
      # Global qq-plot (on gene level since aggregate=TRUE)
      print("Plotting qq-plot")
      print(plot_qq_plot_file)
      jpeg(filename=plot_qq_plot_file, width=800, height=800)
      print(FRASER::plotQQ(fds, aggregate=TRUE, global=TRUE))
      dev.off()
    },
    error = function(e) {
        print("Error in creating qq-plot")
        print(e)
    }
  )
}

main_function()
