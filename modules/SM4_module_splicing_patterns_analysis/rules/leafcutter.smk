import os

# Output directory for all generated files of this Tool
leafcutter_output_path = os.path.join('output', config["output_directories"]["leafcutter_output_dir"])


# -------- 0. Installation of Leafcutter ---------
rule install_leafcutter:
    output:
        os.path.join(config["leafcutter_settings"]["leafcutter_installation_dir"],
            "leafcutter/scripts/leafcutter_ds.R")
    params:
        installation_base_dir=config["leafcutter_settings"]["leafcutter_installation_dir"],
        installation_leafcutter_dir=os.path.join(config["leafcutter_settings"]["leafcutter_installation_dir"],
            "leafcutter"),
        git_url="https://github.com/davidaknowles/leafcutter"
    conda:
        "../envs/leafcutter_env.yaml"
    shell:
        # Check first whether the installation directory already exists -> if yes, then remove
        "if [ -d {params.installation_leafcutter_dir} ]; "
        "   then echo 'Leafcutter directory already exists, but will be removed and reinstalled.' && "
        " rm -rf {params.installation_leafcutter_dir};"
        "fi; "
        # Check whether the base installation directory exists
        "if [ -d {params.installation_base_dir} ]; "
        "   then echo 'Base installation directory already exists, so leafcutter will be installed there.'; "
        "   else echo 'Base installation directory does not exist, so it will be created.' && "
        " mkdir -p {params.installation_base_dir}; "
        "fi; "
        # Clone the leafcutter git repository and install leafcutter
        "cd {params.installation_base_dir}; "
        "git clone {params.git_url}; "
        "R -e 'devtools::install_github(\"davidaknowles/leafcutter/leafcutter\")'; "


# ----------------- 1. Preparation ------------------

def get_regtools_strandedness_parameter(sample_id):
    """
    Returns the library type for the given sample id.
    """
    strandedness = pep.sample_table[pep.sample_table["sample_name"] == sample_id]["stranded"][0]
    if strandedness == "no":
        return "0"
    elif strandedness == "reverse":
        return "1"
    elif strandedness == "yes":
        return "2"
    else:
        raise ValueError("Strandedness parameter not recognized: {}".format(strandedness))


# -------- 1.1 Converting BAMs to JUNCs --------
# Documentation of output: https://regtools.readthedocs.io/en/latest/commands/junctions-extract/
# Important: Olego- and Regtools JUNC-files are different formatted -> Regtool's files include Junction Anchors!
# -> But they have the same results!
rule leafcutter_regtools_bam_to_junc:
    """ Per sample/BAM-file """
    input:
        sorted_bam=main_helper_get_bam_file_for_sample(config["input_dir_of_bam_files"],
            filename_extension=config["bam_files_attributes"]["filename_extension"]),
        sorted_bai=main_helper_get_bam_file_for_sample(config["input_dir_of_bam_files"],
            filename_extension=config["bam_files_attributes"]["filename_extension"]) + ".bai"
    output:
        os.path.join(leafcutter_output_path, "junctions", "{sample_id}.junc")
    log:
        "logs/regtools/{sample_id}.log"
    params:
        a=config["leafcutter_settings"]["regtools_settings"]["regtools_junctions_anchor_length"],
        m=config["leafcutter_settings"]["regtools_settings"]["regtools_junctions_minimum_intron_length"],
        M=config["leafcutter_settings"]["regtools_settings"]["regtools_junctions_maximum_intron_length"],
        strandedness=lambda wildcards: get_regtools_strandedness_parameter(wildcards.sample_id),
    threads:
        1
    conda:
        "../envs/leafcutter_env.yaml"
    shell:
        "regtools junctions extract "
        "-a {params.a} "
        "-m {params.m} "
        "-M {params.M} "
        "-s {params.strandedness} "
        "-o {output} "
        "{input.sorted_bam} 2>{log}"


# -------- 1.2 Assert that for all samples a junc-file was created --------
# Use default resources
rule leafcutter_assert_production_of_junc_files:
    """
    Simply writes down the list of input files into a given output-txt-file
    """
    input:
        bam_files=expand(os.path.join(leafcutter_output_path, "junctions", "{sample_id}.junc"),
            sample_id=pep.sample_table["sample_name"])
    output:
        os.path.join(leafcutter_output_path, config["leafcutter_settings"]["leafcutter_project_output_dir"],
                     "juncfiles.txt")
    threads:
        1
    run:
        with open(output[0], "w") as f:        # Creates file / Rewrites the given file
            for file in input.bam_files:
                f.write(file + "\n")


# ======================== 2. Leafcutter ========================
# -------- 2.1 Intron clustering -> Uses all samples from input-file --------
# -> Filtering according to groups is applied later...
# Use default resources
rule leafcutter_clustering:
    """
    1. Cluster introns according to overlaps.
    Uses all samples from input-file. 
    The filtering according to groups is later applied.
    2. Unzip the output file
    3. Transform Leafcutter clustering results to BED file 
    (# IMPORTANT: Always includes all introns from all input-files!!! (-> Not group-specific))
    """
    input:
        rules.install_leafcutter.output,
        juncfile_list=rules.leafcutter_assert_production_of_junc_files.output
    output:
        output_dir=directory(os.path.join(leafcutter_output_path,
            config["leafcutter_settings"]["leafcutter_project_output_dir"], "clustering")),
        count_file=os.path.join(leafcutter_output_path,
            config["leafcutter_settings"]["leafcutter_project_output_dir"],
            "clustering", config["leafcutter_settings"]["output_prefix"] + "_perind_numers.counts.gz"),
        unzipped_count_file=os.path.join(leafcutter_output_path,
            config["leafcutter_settings"]["leafcutter_project_output_dir"], "clustering", "all_perind_numers.counts"),
        bed_file=os.path.join(leafcutter_output_path, config["leafcutter_settings"]["leafcutter_project_output_dir"],
            "clustering", "leafcutter_introns.bed")
    log:
        "logs/leafcutter/clustering.log"
    params:
        leafcutter_clustering_script_path=os.path.join(config["leafcutter_settings"]["leafcutter_installation_dir"],
            "leafcutter", "clustering", "leafcutter_cluster_regtools.py"),
        m=config["leafcutter_settings"]["leafcutter_clustering_supporting_split_reads"],   # m split reads needed so that a cluster is considered
        l=config["leafcutter_settings"]["leafcutter_clustering_maximum_intron_length"],    # max intron length
        o="all",                                                    # output prefix
        # r=            # output dir prefix is given in output
        bed_file_create_script=workflow.source_path("../scripts/leafcutter/transfer_leafcutter_introns_to_bed_file.py")  # Script to create BED file
    threads:
        1
    conda:
        "../envs/py2_env.yaml"
    shell:
        "mkdir -p {output.output_dir}; "
        "python {params.leafcutter_clustering_script_path} "
        "-j {input.juncfile_list} "    # Juncfile list
        "-m {params.m} "    # m split reads needed so that a cluster is considered
        "-l {params.l} "    # max intron length
        "-o {params.o} "    # output prefix
        "-r {output.output_dir} "    # output dir
        "2>{log} && "
        "gunzip --keep --force {output.count_file} 2>{log} && "
        "python {params.bed_file_create_script} "
        "--leafcutter_intron_file {output.unzipped_count_file} --bed_file {output.bed_file} 2>{log}"


# -------- 2.2 Leafcutter & Leafviz: Annotate reference sequence --------
rule leafcutter_annotate_reference_seq:
    """
    Annotate reference sequence with Leafcutter.
    The output file can then used by Leafcutter to label the respective clusters.
    """
    input:
        rules.install_leafcutter.output,
        reference_genome_annotation_file=config["leafcutter_settings"]["reference_genome_annotation_file"]
    output:
        output_dir=directory(os.path.join(leafcutter_output_path, "leafviz")),
        exon_file=os.path.join(leafcutter_output_path, "leafviz/annotated_genecode_all_exons.txt.gz")
    log:
        "logs/leafcutter/annotate_reference_seq.log"
    params:
        leafcutter_installation_dir=os.path.join(config["leafcutter_settings"]["leafcutter_installation_dir"],
            "leafcutter"),
        output_prefix="annotated_genecode"
    conda:
        "../envs/leafviz_env.yaml"
    threads:
        1
    shell:
        "{params.leafcutter_installation_dir}/leafviz/gtf2leafcutter.pl -o {output.output_dir}/{params.output_prefix} "
        "{input.reference_genome_annotation_file} 2>{log}"


# -------- 2.3 Leafcutter: Create for each condition group a separate group-file --------
# Conditions to query
all_conditions = pep.sample_table[(pep.sample_table["control"] != "true")]["condition"].unique()

rule leafcutter_create_group_files:
    """
    Create for each condition group a separate group-file.
    2 columns (no column names): sample_name, condition (control, or patient)
    """
    output:
        expand(os.path.join(leafcutter_output_path, config["leafcutter_settings"]["leafcutter_project_output_dir"],
               "{condition}_group_file.txt"), condition=all_conditions)
    log:
        "logs/leafcutter/annotate_reference_seq.log"
    params:
        output_dir=lambda wildcards, output: os.path.dirname(output[0]),
        control_samples=pep.sample_table[(pep.sample_table["condition"] == "None") | (pep.sample_table["condition"] == "")],
        condition_samples_array=[
            pep.sample_table[(pep.sample_table["condition"] == condition)]
            for condition in all_conditions
        ]
    threads:
        1
    script:
        "../scripts/leafcutter/create_leafcutter_group_files.py"



# ==================== 3. Differential splicing analysis =====================
# Use 32 CPUs & 131 GB
# -------- 3.1 Differential intron excision analysis --------
# Code: https://github.com/davidaknowles/leafcutter/blob/master/scripts/leafcutter_ds.R
rule leafcutter_differential_splicing_analysis:
    input:
        leafcutter_installed=rules.install_leafcutter.output,
        # This also includes the clustering output
        # Needed so that leafcutter can annotate pooled_intervals according to the corresponding gene
        ref_seq_annotation=rules.leafcutter_annotate_reference_seq.output.exon_file,
        clustering_output=rules.leafcutter_clustering.output.count_file,            # Counts file
        leafcutter_group_file=os.path.join(leafcutter_output_path,
            config["leafcutter_settings"]["leafcutter_project_output_dir"], "{condition}_group_file.txt")
    output:
        os.path.join(leafcutter_output_path, config["leafcutter_settings"]["leafcutter_project_output_dir"],
                     "diff_splicing_analysis/{condition}/sample_analysis_cluster_significance.txt"),
        os.path.join(leafcutter_output_path, config["leafcutter_settings"]["leafcutter_project_output_dir"],
                     "diff_splicing_analysis/{condition}/sample_analysis_effect_sizes.txt"),
    log:
        "logs/leafcutter/{condition}_differential_splicing_analysis.log"
    params:
        leafcutter_installation_dir=os.path.join(config["leafcutter_settings"]["leafcutter_installation_dir"],
            "leafcutter"),
        max_cluster_size=config["leafcutter_settings"]["leafcutter_max_cluster_size"],                 # Don't test clusters with more than this many introns
        min_samples_per_intron=config["leafcutter_settings"]["leafcutter_min_samples_per_intron"],     # Ignore introns used in fewer than n samples
        min_samples_per_group=config["leafcutter_settings"]["leafcutter_min_samples_per_group"],                           # Require this many samples in each group to have at least min_coverage reads
        min_coverage=config["leafcutter_settings"]["leafcutter_min_coverage"],                         # Require min_samples_per_group samples in each group to have at least this many reads
        output_prefix=lambda wildcards, output:
            os.path.join(os.path.dirname(output[0]), "sample_analysis")  # Output-prefix
    resources:
        mem_mb=32*4096,          # total: assign 4096MB per CPU = 131072
        cpus=32,                # uses 32 cpus
        time_min=120            # give max 2 hrs
    threads:
        32
    conda:
        "../envs/leafcutter_env.yaml"
    shell:
        # "R -e 'devtools::install_github(\"stan-dev/rstantools\")'; "
        "R -e 'devtools::install_github(\"davidaknowles/leafcutter/leafcutter\")'; "
        "{params.leafcutter_installation_dir}/scripts/leafcutter_ds.R "
        "--num_threads {threads} "
        "--max_cluster_size {params.max_cluster_size} "             # Don't test clusters with more than this many introns
        "--min_samples_per_intron {params.min_samples_per_intron} " # Ignore introns used if fewer than n samples (i.e. at least one supporting read)
        "--min_samples_per_group {params.min_samples_per_group} "   # Require this many samples in each group to have at least min_coverage reads
        "--min_coverage {params.min_coverage} "                     # Require min_samples_per_group samples in each group to have at least this many reads
        "--output_prefix {params.output_prefix} "           # Output-prefix: Prefix for the 2 output files
        "--exon_file {input.ref_seq_annotation}"                    # File defining known exons, example in data/gencode19_exons.txt.gz. Columns should be chr, start, end, strand, gene_name. Optional, only just to label the clusters
        " {input.clustering_output} {input.leafcutter_group_file} 2>{log}"


rule leafcutter_differential_splicing_analysis_extract_significant_clusters:
    """
    Extract significant clusters from the leafcutter differential splicing analysis.
    Also apply sorting.
    """
    input:
        os.path.join(leafcutter_output_path,config["leafcutter_settings"]["leafcutter_project_output_dir"],
            "diff_splicing_analysis/{condition}/sample_analysis_cluster_significance.txt"),
    output:
        os.path.join(leafcutter_output_path,config["leafcutter_settings"]["leafcutter_project_output_dir"],
            "diff_splicing_analysis/{condition}/sample_analysis_cluster_significance_extracted_significant_p010.tsv"),
    resources:
        mem_mb=4*4096,          # total: assign 4096MB per CPU = 131072
        cpus=1,                # uses 32 cpus
        time_min=30            # give max 2 hrs
    threads:
        1
    run:
        # Extract significant clusters (column: p.adjust)
        # p < 0.1
        import pandas as pd
        df = pd.read_csv(input[0], sep="\t")
        df = df[df["p.adjust"] < 0.1]
        df = df.sort_values(by="p.adjust")
        df.to_csv(output[0], sep="\t", index=False)


rule leafcutter_differential_splicing_analysis_create_report_html:
    input:
        rules.leafcutter_differential_splicing_analysis_extract_significant_clusters.output,
    output:
        report(
            directory(os.path.join(leafcutter_output_path,"html_reports", "leafcutter_splicing_analysis", "{condition}")),
            patterns=["{filename}.report.html"],
            # Path relative to Snakefile
            caption="../../../report_source/module4_splicing_patterns_analysis/leafcutter/caption.rst",
            category=config["output_dir_module4_splicing_patterns_analysis"].replace("output/", ""),
                subcategory="Leafcutter",
                labels={  # Set labels manually
                    "Tool": "Leafcutter",
                    "Condition": "{condition}",
                    "File": "{filename}"
                }
        ),
        os.path.join(leafcutter_output_path, "html_reports", "leafcutter_splicing_analysis", "{condition}",
            "sample_analysis_cluster_significance_extracted_significant_p010.tsv.report.html")
    params:
        input_files=lambda wildcards, input: input,
        data_separators=["\t" for i in range(1)],
        data_titles=["Leafcutter: Clusters with significant differential splicing for condition {condition}"],
        info_texts=[
            workflow.source_path("../../../report_source/module4_splicing_patterns_analysis/leafcutter/info.html")
        ],
        html_output_dir=lambda wildcards, output: output[0],
        html_output_file_basenames=lambda wildcards, output: os.path.basename(output[1])
    conda:
        "../../../report_source/report_env.yaml"
    script:
        "../../../scripts/create_report_html_files.R"


# Annotates each intron in each cluster at a given false discovery rate
rule leafcutter_leafviz_file_preparation:
    input:
        rules.leafcutter_differential_splicing_analysis.output,
        group_file=os.path.join(leafcutter_output_path, config["leafcutter_settings"]["leafcutter_project_output_dir"],
               "{condition}_group_file.txt"),
        count_file=rules.leafcutter_clustering.output.count_file,
        sample_analysis_cluster_significance_file=rules.leafcutter_differential_splicing_analysis.output[0],
        sample_analysis_effect_sizes_file=rules.leafcutter_differential_splicing_analysis.output[1],

        # Annotated gene code file -> Exon file as example...
        annotated_genecode_exon_file=os.path.join(leafcutter_output_path, "leafviz/annotated_genecode_all_exons.txt.gz")
    output:
        os.path.join(leafcutter_output_path,config["leafcutter_settings"]["leafcutter_project_output_dir"],
            "diff_splicing_analysis/{condition}/leafviz_obj.RData"),
    params:
        leafcutter_installation_dir=os.path.join(
            config["leafcutter_settings"]["leafcutter_installation_dir"],"leafcutter"),
        annotation_files_prefix=lambda wildcards, input:
                os.path.join(os.path.dirname(input.annotated_genecode_exon_file),"annotated_genecode"),
        fdr=0.10
    conda:
        "../envs/leafcutter_env.yaml"
    threads:
        1
    resources:
        mem_mb=4 * 4096,    # total: assign 4096MB per CPU = 131072
        cpus=1,             # uses 1 cpus
        time_min=30         # give max 30 min
    shell:
        "{params.leafcutter_installation_dir}/leafviz/prepare_results.R "
        "-m {input.group_file} "
        "-f {params.fdr} "
        "{input.count_file} "
        "{input.sample_analysis_cluster_significance_file} {input.sample_analysis_effect_sizes_file} "
        "{params.annotation_files_prefix} "
        "-o {output}"


# ------------------------------------------------------------------------------------------------------------------------

# ========== LeafcutterMD ==============
# ------------ 5. Outlier intron excision analysis ----------------
# ATTENTION: This is calculated for all given samples, since clusters are build with all samples!!!
# ATTENTION: This returns non-adjusted p-values!
rule leafcutterMD_outlier_intron_excision_analysis:
    input:
        rules.leafcutter_clustering.output.count_file                            # Counts file
    output:
        # Attention: Output tables include index column
        # all_outlier_introns_pVals_file: <prefix>_pVals.txt (containing p-value for each intron)
        os.path.join(leafcutter_output_path, config["leafcutter_settings"]["leafcutter_project_output_dir"],
            "diff_splicing_analysis_outliers/all_outliers_pVals.txt"),
        # all_outlier_clusters_pVals_file: <prefix>_clusterPvals.txt (containing p-value for each cluster)
        os.path.join(leafcutter_output_path, config["leafcutter_settings"]["leafcutter_project_output_dir"],
            "diff_splicing_analysis_outliers/all_outliers_clusterPvals.txt"),
        # all_outlier_effSize_file: <prefix>_effSize.txt (containing the effect sizes for each intron)
        os.path.join(leafcutter_output_path, config["leafcutter_settings"]["leafcutter_project_output_dir"],
            "diff_splicing_analysis_outliers/all_outliers_effSize.txt"),
        # Output directory
        output_dir=directory(os.path.join(leafcutter_output_path, config["leafcutter_settings"]["leafcutter_project_output_dir"],
            "diff_splicing_analysis_outliers"))
    log:
        "logs/leafcutter/outlier_analysis.log"
    params:
        leafcutter_installation_dir=os.path.join(config["leafcutter_settings"]["leafcutter_installation_dir"],
                                                 "leafcutter"),
        output_prefix="all_outliers"
    resources:
        mem_mb=32*4096,     # total: assign 4096MB per CPU = 131072
        cpus=32,            # uses 32 cpus
        time_min=120        # give max 2 hrs
    threads:
        32
    conda:
        "../envs/leafcutter_env.yaml"
    shell:
        "R -e 'devtools::install_github(\"davidaknowles/leafcutter/leafcutter\")'; "
        "{params.leafcutter_installation_dir}/scripts/leafcutterMD.R --num_threads {threads} "
        "--output_prefix {output.output_dir}/{params.output_prefix} {input} 2>{log}"


# -------- 6. Analyze LeafcutterMD output --------
rule leafcutterMD_output_analysis:
    """
    Analyze LeafcutterMD output (for affected samples)
    Currently included:
        - Extract only non-control samples
        - Extract outliers with p-value < 0.01
        - Extract only top 1000 outliers
    """
    input:
        # Attention: Input tables include index column
        all_outlier_introns_pVals_file = rules.leafcutterMD_outlier_intron_excision_analysis.output[0],
        all_outlier_clusters_pVals_file = rules.leafcutterMD_outlier_intron_excision_analysis.output[1],
        all_outlier_effSize_file = rules.leafcutterMD_outlier_intron_excision_analysis.output[2]
    output:
        all_filtered_introns_file=os.path.join(leafcutter_output_path,
            config["leafcutter_settings"]["leafcutter_project_output_dir"], "diff_splicing_analysis_outliers/filtered_all_introns_pvals.tsv"),
        condition_filtered_introns_file=os.path.join(leafcutter_output_path,
            config["leafcutter_settings"]["leafcutter_project_output_dir"], "diff_splicing_analysis_outliers/filtered_condition_introns_pvals.tsv"),
        all_filtered_clusters_file=os.path.join(leafcutter_output_path,
            config["leafcutter_settings"]["leafcutter_project_output_dir"], "diff_splicing_analysis_outliers/filtered_all_clusters_pvals.tsv"),
        condition_filtered_clusters_file=os.path.join(leafcutter_output_path,
            config["leafcutter_settings"]["leafcutter_project_output_dir"], "diff_splicing_analysis_outliers/filtered_condition_clusters_pvals.tsv"),
    log:
        "logs/leafcutter/outlier_analysis.log"
    params:
        patient_sample_ids=(pep.sample_table[pep.sample_table["control"] == "false"]["sample_name"]).tolist(),
        pvalue_threshold=0.01
    threads:
        1
    conda:
        "../envs/pandas_env.yaml"
    script:
        "../scripts/leafcutter/analyze_leafcutterMD_output.py"


rule leafcutter_create_html_reports_for_leafcutterMD:
    input:
        condition_filtered_introns_file = rules.leafcutterMD_output_analysis.output.condition_filtered_introns_file,
        condition_filtered_clusters_file = rules.leafcutterMD_output_analysis.output.condition_filtered_clusters_file,
    output:
        report(
            directory(os.path.join(leafcutter_output_path,"html_reports", "leafcutterMD_outlier_analysis")),
            patterns=["{filename}.tsv.report.html"],
            # Path relative to Snakefile
            caption="../../../report_source/module4_splicing_patterns_analysis/leafcutter/caption.rst",
            category=config["output_dir_module4_splicing_patterns_analysis"].replace("output/",""),
            subcategory="LeafcutterMD",
            labels={  # Set labels manually
                "File:": "{filename}"
            }
        ),
        # os.path.join(leafcutter_output_path, "html_reports", "leafcutterMD_outlier_analysis",
        #     "filtered_all_introns_pvals.tsv.report.html"),
        os.path.join(leafcutter_output_path,"html_reports","leafcutterMD_outlier_analysis",
            "filtered_condition_introns_pvals.tsv.report.html"),
        # os.path.join(leafcutter_output_path,"html_reports","leafcutterMD_outlier_analysis",
        #     "filtered_all_clusters_pvals.tsv.report.html"),
        os.path.join(leafcutter_output_path,"html_reports","leafcutterMD_outlier_analysis",
            "filtered_condition_clusters_pvals.tsv.report.html"),
    params:
        input_files=lambda wildcards, input: input,
        data_separators=["\t" for i in range(2)],
        data_titles=["LeafcutterMD results" for i in range(2)],
        info_texts=["Only top 1000 results are displayed!" for i in range(2)],
        html_output_dir=lambda wildcards, output: output[0],
        html_output_file_basenames=lambda wildcards, input: [
            (os.path.basename(input_file) + ".report.html") for input_file in input
        ]
    conda:
        "../../../report_source/report_env.yaml"
    script:
        "../../../scripts/create_report_html_files.R"
