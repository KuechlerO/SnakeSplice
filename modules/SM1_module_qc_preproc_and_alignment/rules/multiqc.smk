# ========== Submodule: MultiQC - Summarize Quality-Control Steps ==========
import os

# Output directory
multiqc_output_path = os.path.join('output',config["output_directories"]["multiqc_output_dir"])

# --------- 0: Depending on user's choice, identify reports that should be included in MultiQC report. ---------
multiqc_input_files = []
if config["switch_variables"]["run_kraken2"]:
    kraken2_output_path = os.path.join('output',config["output_directories"]["kraken2_output_dir"])
    multiqc_input_files.append(expand(os.path.join(kraken2_output_path, "{sample_id}_report_summary.report"),
                                      sample_id=pep.sample_table["sample_name"]))

if config["switch_variables"]["run_fastqc_before_trimming"]:
    fastqc_output_path = os.path.join('output',config["output_directories"]["fastqc_output_dir"], "non_trimmed")
    multiqc_input_files.append(expand(os.path.join(fastqc_output_path, "{sample_id}.html"),
                                      sample_id=pep.sample_table["read1"]))
    multiqc_input_files.append(expand(os.path.join(fastqc_output_path,"{sample_id}.html"),
                                      sample_id=pep.sample_table["read2"]))

if config["switch_variables"]["run_fastqc_after_trimming"]:
    fastqc_output_path = os.path.join('output',config["output_directories"]["fastqc_output_dir"], "trimmed")
    multiqc_input_files.append(expand(os.path.join(fastqc_output_path, "{sample_id}_R1.html"),
                                      sample_id=pep.sample_table["sample_name"]))
    multiqc_input_files.append(expand(os.path.join(fastqc_output_path, "{sample_id}_R2.html"),
                                      sample_id=pep.sample_table["sample_name"]))

if config["switch_variables"]["run_qualimap"] and config["switch_variables"]["run_alignment"]["use_star"]:
    qualimap_output_path = os.path.join('output',config["output_directories"]["qualimap_output_dir"], "star")
    multiqc_input_files.append(expand(os.path.join(qualimap_output_path, "{sample_id}/qualimapReport.html"),
                                      sample_id=pep.sample_table["sample_name"]))

if config["switch_variables"]["run_qualimap"] and config["switch_variables"]["run_alignment"]["use_olego"]:
    qualimap_output_path = os.path.join('output',config["output_directories"]["qualimap_output_dir"], "olego")
    multiqc_input_files.append(expand(os.path.join(qualimap_output_path,"{sample_id}/qualimapReport.html"),
                                      sample_id=pep.sample_table["sample_name"]))


# 1: Create MultiQC report.
# Bind all fastqc & Qualimap statistics & Kraken together
# Use default resources: 1 cpu & 2 GB RAM
rule run_multiqc:
    input:
        multiqc_input_files
    output:
        report(
            os.path.join(multiqc_output_path,"multiqc.html"),
            # Path relative to Snakefile
            caption="../../../report_source/module1_qc_preproc_alignment/multiqc/caption.rst",
            category=config["output_dir_module1_qc_preprocessing_and_alignment"].replace("output/",""),
            subcategory="MultiQC",
        )
    params:
        "-d"  	# -d: Uses urls as sample-names -> Needed to include both: fastqc before and after trimming
    log:
        "logs/multiqc.log"
    wrapper:
        "v0.86.0/bio/multiqc"

