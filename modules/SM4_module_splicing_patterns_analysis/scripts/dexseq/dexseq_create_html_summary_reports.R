library(BiocParallel)   # For parallelization needed
library(DEXSeq)         # For differential expression analysis

library(dplyr)


# -------------------- Main function --------------------
main_function <- function() {
  # ----------------- 1. Load snakemake variables -----------------
  # inputs
  input_dxr_object_file_list <- snakemake@input[["dexseq_results_object_file_list"]]

  # outputs
  output_html_report_file_array <- snakemake@output[["result_html_summary_report_file_list"]]

  # params
  summary_report_fdr <- snakemake@params[["summary_report_fdr"]]
  threads <- snakemake@threads

  # ----------------- 2. Prepare analysis -----------------
  # ------ 2.1 Set number of threads --------
  BPPARAM <- MulticoreParam(threads)


  # ----------------- 3. Run analysis -----------------
  # ------ 3.1 Create html summary report for condition group --------
  for (i in c(1:length(input_dxr_object_file_list))) {
    # ------ 3.1 Load DEXSeq results object --------
    current_dexseq_result_obj <- readRDS(file=input_dxr_object_file_list[i])

    # Output directory & file
    current_output_html_report_file <- output_html_report_file_array[i]
    current_output_html_report_dir <- dirname(current_output_html_report_file)

    # HTML Summary with linkouts
    tryCatch(
      expr = {
        print("Creating HTML report")
        print("1. Create dir")
        dir.create(current_output_html_report_dir, showWarnings = FALSE)

        print("2. Create HTML report")
        DEXSeqHTML(current_dexseq_result_obj, FDR=summary_report_fdr, path=current_output_html_report_dir,
                   file=basename(current_output_html_report_file), BPPARAM=BPPARAM)
      },
      error = function(e) {
        print("Error in DEXSeqHTML")
        print(e)
        print("Saving Error instead of HTML report")
        sink(file=current_output_html_report_file); print(e); sink()
      },
      finally = {
          print("Finished DEXSeqHTML")
      }
    )
  }
}

# -------------------- Run main function --------------------
main_function()
