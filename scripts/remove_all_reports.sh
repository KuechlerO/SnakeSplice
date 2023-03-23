#!/bin/bash

# Module 0 reports
rm output/module0_report_generation/output/input_file_overview.html

# Module 1 reports
# Strandedness Check
rm -r output/module1_qc_preproc_and_alignment/output/check_strandedness/html_reports
# Samtools Statistics
rm -r output/module1_qc_preproc_and_alignment/output/bamstats/*/html_reports
# MultiQC
rm -r output/module1_qc_preproc_and_alignment/output/multiqc/multiqc.html

# Module 2 reports
# Arriba
rm -r output/module2_gene_fusion_detection/output/arriba/html_reports

# Module 3 reports
# Nothing yet...

# Module 4 reports
# Outrider
rm -r output/module3_gene_expression_quantification_and_analysis/output/outrider/html_reports
# Salmon
rm -r output/module3_gene_expression_quantification_and_analysis/output/salmon/*/html_reports
rm -r output/module3_gene_expression_quantification_and_analysis/output/salmon/*/gsea_analysis/html_reports/
# Kallisto
rm -r output/module3_gene_expression_quantification_and_analysis/output/kallisto/*/html_reports
rm -r output/module3_gene_expression_quantification_and_analysis/output/kallisto/*/gsea_analysis/html_reports/

# Module 5 reports
# DexSeq
rm -r output/module4_splicing_patterns_analysis/output/dexseq/html_reports
# FRASER
rm -r output/module4_splicing_patterns_analysis/output/fraser/html_reports
# Leafcutter
rm -r output/module4_splicing_patterns_analysis/output/leafcutter/html_reports/leafcutter_splicing_analysis
# LeafcutterMD
rm -r output/module4_splicing_patterns_analysis/output/leafcutter/html_reports/leafcutterMD_outlier_analysis
# PJD
rm -r output/module4_splicing_patterns_analysis/output/private_junction_detection/html_reports
# rmats
rm -r output/module4_splicing_patterns_analysis/output/rmats/*/html_reports
# summary
rm -r output/module4_splicing_patterns_analysis/output/overall_summary/html_reports



