import os

# Output directory for all generated files of this Tool
pjd_output_path = os.path.join('output', config["output_directories"]["private_junction_detection_output_dir"])

def pjd_extract_sample_ids_with_specific_condition(wildcards):
    """
    Depending on the condition-wildcard, this function returns a list of sample ids
    :param wildcards:
    :return:
    """
    condition = wildcards.condition
    return pep.sample_table[pep.sample_table["condition"] == condition]["sample_name"].tolist()


def pjd_get_regtools_strandedness_parameter(sample_id):
    """
    Returns the strandedness parameter for regtools junctions extract
    :param sample_id:
    :return:
    """
    strandedness = pep.sample_table[pep.sample_table["sample_name"] == sample_id]["stranded"][0]
    if strandedness == "no":
        return "0"
    elif strandedness == "reverse":
        return "1"
    elif strandedness == "yes":
        return "2"
    else:
        raise ValueError("Strandedness value not recognized: " + strandedness)


rule pjd_regtools_bam_to_junc:
    """
    Extracts junctions from BAM-files
    """
    input:
        # Return URLs of BAM-files with {sample_id}-wildcard
        bam_file=main_helper_get_bam_file_for_sample(config["input_dir_of_bam_files"],
            filename_extension=config["bam_files_attributes"]["filename_extension"]),
        bam_index=main_helper_get_bam_file_for_sample(config["input_dir_of_bam_files"],
            filename_extension=config["bam_files_attributes"]["filename_extension"]) + ".bai"
    output:
        os.path.join(pjd_output_path, "regtools/{sample_id}.bam.junc")
    log:
        "logs/private_junction_detection/regtools/{sample_id}.log"
    params:
        a=config["private_junction_detection_settings"]["regtools_junctions_anchor_length"],
        m=config["private_junction_detection_settings"]["regtools_junctions_minimum_intron_length"],
        M=config["private_junction_detection_settings"]["regtools_junctions_maximum_intron_length"],
        strandedness=lambda wildcards: get_regtools_strandedness_parameter(wildcards.sample_id)
    threads:
        1
    conda:
        "../envs/private_junction_detection_env.yaml"
    shell:
        "regtools junctions extract "
        "-a {params.a} "
        "-m {params.m} "
        "-M {params.M} "
        "-s {params.strandedness} "
        "-o {output} "
        "{input.bam_file} 2>{log}"


rule pjd_extract_actual_junctions:
    """
    Based on: https://regtools.readthedocs.io/en/latest/commands/junctions-extract/
    Extracts the actual junctions from the .bam.junc file
    """
    input:
        # First line causes bug of DAG computation taking too long...
        #regtools_junc_files=os.path.join(pjd_output_path, "regtools/{sample_id}.bam.junc")
        regtools_junc_files=rules.pjd_regtools_bam_to_junc.output[0]
    output:
        output_file=os.path.join(pjd_output_path, "{sample_id}.bam.junc")
    log:
        "logs/private_junction_detection/extract_junctions_{sample_id}.log"
    threads:
        1
    resources:
        mem_mb=4096*4,
        cpus=1,
        time_min=60
    conda:
        "../envs/private_junction_detection_env.yaml"
    script:
        "../scripts/private_junction_detection/extract_actual_junctions.py"


rule pjd_sort_junc_file:
    """
    Sorts the .junc file, which is needed for Bedtools closest
    """
    input:
        in_file=rules.pjd_extract_actual_junctions.output.output_file
    output:
        output_file=os.path.join(pjd_output_path,"{sample_id}.bam.extracted.sorted.junc")
    log:
        "logs/private_junction_detection/{sample_id}_sorting_juncfile.log"
    threads:
        1
    resources:
        mem_mb=4096*4,
        cpus=1,
        time_min=60
    conda:
        "../envs/private_junction_detection_env.yaml"
    shell:
        "sort-bed {input.in_file} >{output.output_file} 2>{log}"


rule pjd_transform_gtf_to_bed_file:
    """
    Transform annotation file from GTF to BED format
    -> Only genes are included
    """
    input:
        gtf_file=os.path.join(config["private_junction_detection_settings"]["ensemble_sourced_gtf_file"]),
        chrom_transform_script=workflow.source_path("../../../scripts/transform_chro_names_of_gtf_file.py")
    output:
        bed_file=temp(os.path.join(pjd_output_path,"annotation.bed")),
        bed_file_transformed=os.path.join(pjd_output_path, "transformed_annotation.bed")
    log:
        "logs/private_junction_detection/transform_gtf_to_bed_file.log"
    threads:
        1
    resources:
        mem_mb=4096*4,
        cpus=1,
        time_min=60
    conda:
        "../envs/private_junction_detection_env.yaml"
    shell:
        "gtf2bed <{input.gtf_file} | grep -w gene | sort-bed - >{output.bed_file} 2>{log};"
        "python3 {input.chrom_transform_script} --input_gtf_file {output.bed_file} "
        "--output_gtf_file {output.bed_file_transformed} --log_file {log}"


rule pjd_map_genes_to_junctions:
    """
    Map gene names to coordinates that are given in BED-file
    -> adds closest match as additional columns in BED-file
    """
    input:
        junc_file=rules.pjd_sort_junc_file.output.output_file,
        annotation_file=rules.pjd_transform_gtf_to_bed_file.output.bed_file_transformed
    output:
        os.path.join(pjd_output_path, "{sample_id}.bam.extracted.sorted.named.junc")
    log:
        "logs/private_junction_detection/{sample_id}_mapping_genes_to_junc_files.log"
    params:
        # -d:           Reports also distance in extra column -> 0 for overlapping features
        # -t first:     When distance tie, then report only first match
        extra="-d -t first"          # When distance tie, then report only first match
    threads:
        1
    resources:
        mem_mb=4096*4,
        cpus=1,
        time_min=60
    conda:
        "../envs/private_junction_detection_env.yaml"
    shell:
        "bedtools closest {params.extra} -a {input.junc_file} -b {input.annotation_file} > {output} 2>{log}"


rule pjd_collect_all_junctions:
    input:
        all_junc_files=expand(os.path.join(pjd_output_path,"{sample_id}.bam.extracted.sorted.named.junc"),
            sample_id=pep.sample_table["sample_name"])
    output:
        output_file=os.path.join(pjd_output_path, "all_junctions.junc")
    params:
        sample_names=pep.sample_table["sample_name"]
    log:
        "logs/private_junction_detection/collect_all_junctions.log"
    conda:
        "../envs/private_junction_detection_env.yaml"
    threads:
        1
    resources:
        mem_mb=4096*4,
        cpus=1,
        time_min=60
    script:
        "../scripts/private_junction_detection/junction_collector.py"


rule pjd_filter_junctions_0:
    input:
        junction_collection_file=rules.pjd_collect_all_junctions.output.output_file
    output:
        only_control_junctions_file=os.path.join(pjd_output_path, "{condition}/filtered_junctions_0/only_control_junctions.junc"),
        only_condition_junctions_file=os.path.join(pjd_output_path, "{condition}/filtered_junctions_0/only_condition_junctions.junc")
    log:
        "logs/private_junction_detection/{condition}/filter_junctions_0.log"
    params:
        control_samples=pep.sample_table[(pep.sample_table["control"] == "true")]["sample_name"].tolist(),
        condition_samples=pjd_extract_sample_ids_with_specific_condition,
        max_contrast=0.0
    threads:
        1
    resources:
        mem_mb=4096*4,
        cpus=1,
        time_min=60
    conda:
        "../envs/private_junction_detection_env.yaml"
    script:
        "../scripts/private_junction_detection/filter_junctions.py"


rule pjd_extract_x_entries_filter_junctions_0:
    """
    Extract 250 entries from each junction file
    """
    input:
        rules.pjd_filter_junctions_0.output.only_control_junctions_file,
        rules.pjd_filter_junctions_0.output.only_condition_junctions_file
    output:
        control_junctions_file=os.path.join(pjd_output_path, "{condition}/filtered_junctions_0/only_control_junctions_250.junc"),
        condition_junctions_file=os.path.join(pjd_output_path, "{condition}/filtered_junctions_0/only_condition_junctions_250.junc")
    log:
        "logs/private_junction_detection/{condition}/extract_x_entries.log"
    threads:
        1
    resources:
        mem_mb=4096*4,
        cpus=1,
        time_min=20
    run:
        # Extract 150 entries from each file and write them to output
        for input_file, output_file in zip(input, output):
            with open(input_file, "r") as in_file:
                with open(output_file, "w") as out_file:
                    for i, line in enumerate(in_file):
                        if i < 251:     # 150 entries + header
                            out_file.write(line)
                        else:
                            break


rule pjd_insert_gene_symbol_and_gene_id:
    input:
        rules.pjd_extract_x_entries_filter_junctions_0.output.control_junctions_file,
        rules.pjd_extract_x_entries_filter_junctions_0.output.condition_junctions_file,
    output:
        control_junctions_file=os.path.join(pjd_output_path, "{condition}/filtered_junctions_0/named_only_control_junctions_250.gene_symbol.junc"),
        condition_junctions_file=os.path.join(pjd_output_path, "{condition}/filtered_junctions_0/named_only_condition_junctions_250.gene_symbol.junc")
    params:
        info_col_name="add_info",
        gene_ensembl_id_col_name="gene_id",
        gene_name_col_name="gene_name",
    threads:
        1
    resources:
        mem_mb=4096 * 4,
        cpus=1,
        time_min=20
    conda:
        "../envs/pjd_annotate_genes_env.yaml"
    script:
        "../scripts/private_junction_detection/insert_gene_symbol_and_gene_id.R"


rule pjd_create_html_reports:
    input:
        rules.pjd_insert_gene_symbol_and_gene_id.output.control_junctions_file,
        rules.pjd_insert_gene_symbol_and_gene_id.output.condition_junctions_file,
        # rules.pjd_filter_junctions_1.output,
    output:
        report(
            directory(os.path.join(pjd_output_path, "html_reports", "{condition}")),
            patterns=["{source}-junctions-{contrast}.junc.report.html"],
            # Path relative to Snakefile
            caption="../../../report_source/module4_splicing_patterns_analysis/pjd/caption.rst",
            category=config["output_dir_module4_splicing_patterns_analysis"].replace("output/",""),
            subcategory="Private Junction Detection",
            labels={  # Set labels manually
                "condition": "{condition}",
                "source": "{source}",
                "File:": "{source}-junctions-{contrast}.junc"
            }
        ),
        os.path.join(pjd_output_path, "html_reports", "{condition}", "only_control-junctions-000.junc.report.html"),
        os.path.join(pjd_output_path, "html_reports", "{condition}", "only_condition-junctions-000.junc.report.html"),
        #os.path.join(pjd_output_path, "html_reports", "{condition}", "only_control_junctions-010.junc.report.html"),
        #os.path.join(pjd_output_path, "html_reports", "{condition}", "only_condition_junctions-010.junc.report.html")
    params:
        input_files=lambda wildcards, input: input,
        data_separators=["\t" for i in range(2)],
        data_titles=lambda wildcards,output:
            [
                f"Private Junction Detection: Control group compared to the condition '{wildcards.condition}' group",
                f"Private Junction Detection: Junctions only appearing in the condition '{wildcards.condition}' group"
            ],
        # Cannot use lambda function with workflow.source_path -> ToDo Buggy
        info_texts=[
            workflow.source_path("../../../report_source/module4_splicing_patterns_analysis/pjd/info.html")
            for i in range(2)
        ],
        html_output_dir=lambda wildcards, output: output[0],
        html_output_file_basenames=lambda wildcards, output: [os.path.basename(file) for file in output[1:]]
    conda:
        "../../../report_source/report_env.yaml"
    script:
        "../../../scripts/create_report_html_files.R"



# ======================= Not needed anymore =======================
rule pjd_filter_junctions_1:
    input:
        junction_collection_file=rules.pjd_collect_all_junctions.output.output_file
    output:
        only_control_junctions_file=os.path.join(pjd_output_path, "{condition}/filtered_junctions_1/only_control_junctions.junc"),
        only_condition_junctions_file=os.path.join(pjd_output_path, "{condition}/filtered_junctions_1/only_condition_junctions.junc")
    log:
        "logs/private_junction_detection/{condition}/filter_junctions_1.log"
    params:
        control_samples=pep.sample_table[(pep.sample_table["control"] == "true")]["sample_name"].tolist(),
        condition_samples=pjd_extract_sample_ids_with_specific_condition,
        max_contrast=0.1
    threads:
        1
    resources:
        mem_mb=4096*4,
        cpus=1,
        time_min=60
    conda:
        "../envs/private_junction_detection_env.yaml"
    script:
        "../scripts/private_junction_detection/filter_junctions.py"

