# ===== Sample annotations: Made portable by using PEP ===========
pepfile: "pep/pep_config.yaml"
pepschema: "pep/pep_schema_config.yaml"

# ====== Import Workflow settings in configuration files ============
configfile: "config_files/config_main.yaml"                         # The general config-file
configfile: "config_files/config_module1_qc_preproc_alignment.yaml"  # The config-file for the module1: qc_preproc_alignment
configfile: "config_files/config_module2_gene_fusion.yaml"           # The config-file for the module2: gene_fusion
configfile: "config_files/config_module3_gene_expression_quantification_and_analysis_settings.yaml"     # The config-file for the module4: detection_quantification_transcripts
configfile: "config_files/config_module4_splicing_patterns_analysis.yaml"               # The config-file for the module6: splicing_patterns_analysis


# Create merged dict of main configuration settings
bundled_main_config = config["module_output_dir_names"]

# ====== 1: Import Snakemake modules ============
# 1.0 Module for report generation -> Just small helper rules
module module0_report_generation:
    snakefile:  "modules/SM0_report_generation/Snakefile"
    config:     bundled_main_config                                          # merged dicts
    prefix:     config["module_output_dir_names"]["output_dir_module0_report_generation"]

use rule * from module0_report_generation as module0_Report_Generation_*


# 1.1: Module Quality-Control, Preprocessing, and Alignment
module module1_qc_preprocessing_alignment:
    snakefile:  "modules/SM1_module_qc_preproc_and_alignment/Snakefile"
    config:     bundled_main_config | config["module1_qc_preprocessing_and_alignment_settings"]   # merged dicts
    prefix:     config["module_output_dir_names"]["output_dir_module1_qc_preprocessing_and_alignment"]

use rule * from module1_qc_preprocessing_alignment as module1_QcPreprocAlign_*


# 1.2: Module Detection of gene fusion events
module module2_gene_fusion_detection:
    snakefile: "modules/SM2_module_genefusion_detection/Snakefile"
    config:     bundled_main_config | config["module2_gene_fusion_detection_settings"] # merged dicts
    prefix:     config["module_output_dir_names"]["output_dir_module2_gene_fusion_detection"]

use rule * from module2_gene_fusion_detection as module2_GeneFusionDetection_*


# 1.4: Module Gene Expression Quantification and Analysis
module module3_detection_and_quantification_of_transcripts:
    snakefile: "modules/SM3_module_gene_expression_quantification_and_analysis/Snakefile"
    config:     bundled_main_config | config["module3_gene_expression_quantification_and_analysis_settings"] | config["module_output_dir_names"] # merged dicts
    prefix:     config["module_output_dir_names"]["output_dir_module3_gene_expression_quantification_and_analysis"]

use rule * from module3_detection_and_quantification_of_transcripts as module3_GeneExpressionAndAnalysis_*


# 1.6: Module Splicing Pattern Analysis
module module4_splicing_patterns_analysis:
    snakefile: "modules/SM4_module_splicing_patterns_analysis/Snakefile"
    config:     bundled_main_config | config["module4_splicing_patterns_analysis_settings"] # merged dicts
    prefix:     config["module_output_dir_names"]["output_dir_module4_splicing_patterns_analysis"]

use rule * from module4_splicing_patterns_analysis as module4_SplicingPatternsAnalysis_*


# 1.6: Module Analysis of non-coding RNA
# TODO


# ============== 2: Output depending on user choices in config-file ==============
# ------------ 2.1. Check if chosen submodules are compatible ----------------
# Nothing...

# ------------ 2.2. Collect output-files by referencing to submodules "all/complete_output"-rules  ----------------
conditional_output_files = []

# Module 1: Quality-Control, Preprocessing, and Alignment
if config["module_switches"]["run_module1_qc_preprocessing_and_alignment"]:
    conditional_output_files.append(rules.module1_QcPreprocAlign_complete_output.input) # Main rule of module qc_preprocessing_alignment
# Module 2: Detection of gene fusion events
if config["module_switches"]["run_module2_gene_fusion_detection"]:
    conditional_output_files.append(rules.module2_GeneFusionDetection_complete_output.input)    # Main rule of module gene_fusion_detection
# Module 3: Gene Expression Quantification and Analysis
if config["module_switches"]["run_module3_gene_expression_quantification_and_analysis"]:
    conditional_output_files.append(rules.module3_GeneExpressionAndAnalysis_complete_output.input)    # Main rule of module detection_quantification_transcripts
# Module 4: Splicing Pattern Analysis
if config["module_switches"]["run_module4_splicing_pattern_analysis"]:
    conditional_output_files.append(rules.module4_SplicingPatternsAnalysis_complete_output.input)    # Main rule of module splicing_patterns_analysis



# ======== Main rule ===========
localrules: all         # Defines trivial tasks, that do not need to run on cluster nodes

rule all:
    input:
        # Report generation module -> help outputs
        rules.module0_Report_Generation_complete_output.input,      # Main rule of module report_generation

        # Output of remaining modules
        conditional_output_files,
    default_target: True    # Makes this rule the default rule

# ========= Workflow settings =========
report: "report_source/workflow_overview_text.rst"      # Description shown on the start page of workflow report
