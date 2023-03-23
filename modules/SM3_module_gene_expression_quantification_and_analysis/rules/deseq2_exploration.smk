import os


rule deseq2_exploration_create_deseq_dataset_obj:
    """
    Creates a DESeq2 object from the quantification output
    """
    input:
        annotation_table_file=os.path.join("output", "{quant_tool}", "{condition}", "{quant_tool}_output_annotation_table.tsv")
    output:
        deseq_dataset_r_obj=os.path.join("output", "{quant_tool}", "{condition}", "{quant_tool}_deseq_dataset_object.rds")
    params:
        count_algorithm=lambda wildcards: wildcards.quant_tool,
    conda:
        "../envs/deseq2_gene_level_env.yaml"
    threads:
        12
    resources:
        mem_mb=12 * 4096,# total: assign 4096MB per CPU
        cpus=12,	# uses 12 cpus
        time_min=360		    # use 6 hrs, to make sure
    script:
        "../scripts/deseq2/create_deseq_dataset_object.R"



# ============= Due to not supported wildcards in reports kallisto and salmon rules are defined separately =============
rule deseq2_exploration_salmon_explore_deseq_dataset_obj:
    """
    Explores the deseq dataset object
    -> Creates 2 Heatmaps & 2 PCA plots:
    1. Euclidian distance
    2. Poisson distance
    3. PCA on normalized counts
    4. Generalized PCA using glmPCA
    """
    input:
        deseq_dataset_r_obj=os.path.join("output", "salmon", "{condition}", "salmon_deseq_dataset_object.rds")
    output:
        report(
            directory(os.path.join("output", "salmon", "{condition}", "explore_deseq_dataset_obj")),
            patterns=["{filename}.jpg"],
            # Path relative to Snakefile
            caption="../../../report_source/module3_gene_expression_quantification_and_analysis/salmon/caption.rst",
            category=config["output_dir_module3_gene_expression_quantification_and_analysis"].replace("output/", ""),
            subcategory="Salmon: Explore quantification results",
            labels={  # Set labels manually
                "Tool": "Salmon-DESeq2",
                "Condition": "{condition}",
                "Modus": "Exploration",
                "Type": "Plot",
                "File:": "{filename}"
            }
        ),
        os.path.join("output", "salmon", "{condition}", "explore_deseq_dataset_obj", "sampleEuclidianDistMatrix.jpg"),
        os.path.join("output", "salmon", "{condition}", "explore_deseq_dataset_obj", "samplePoissonDistMatrix.jpg"),
        os.path.join("output", "salmon", "{condition}", "explore_deseq_dataset_obj", "PCA_rlog_transformed.jpg"),
        os.path.join("output", "salmon", "{condition}", "explore_deseq_dataset_obj", "PCA_glmPCA.jpg")
    log:
        "logs/salmon/exploration_{condition}.log",
    params:
        output_file_paths = lambda wildcards, output: output[1:]
    conda:
        "../envs/deseq2_gene_level_env.yaml"
    threads:
        12
    resources:
        mem_mb=12*4096,         # total: assign 4096MB per CPU
        cpus=12,	            # uses 12 cpus
        time_min=360		    # use 6 hrs, to make sure
    script:
        "../scripts/deseq2/explore_deseq_dataset.R"


rule deseq2_exploration_kallisto_explore_deseq_dataset_obj:
    """
    Explores the deseq dataset object
    -> Creates 2 Heatmaps & 2 PCA plots:
    1. Euclidian distance
    2. Poisson distance
    3. PCA on normalized counts
    4. Generalized PCA using glmPCA
    """
    input:
        deseq_dataset_r_obj=os.path.join("output", "kallisto", "{condition}", "kallisto_deseq_dataset_object.rds")
    output:
        report(
            directory(os.path.join("output", "kallisto", "{condition}", "explore_deseq_dataset_obj")),
            patterns=["{filename}.jpg"],
            # Path relative to Snakefile
            caption="../../../report_source/module3_gene_expression_quantification_and_analysis/kallisto/caption.rst",
            category=config["output_dir_module3_gene_expression_quantification_and_analysis"].replace("output/", ""),
            subcategory="Kallisto: Explore quantification results",
            labels={  # Set labels manually
                "Tool": "Kallisto-DESeq2",
                "Condition": "{condition}",
                "Modus": "Exploration",
                "Type": "Plot",
                "File:": "{filename}"
            }
        ),
        os.path.join("output", "kallisto", "{condition}", "explore_deseq_dataset_obj", "sampleEuclidianDistMatrix.jpg"),
        os.path.join("output", "kallisto", "{condition}", "explore_deseq_dataset_obj", "samplePoissonDistMatrix.jpg"),
        os.path.join("output", "kallisto", "{condition}", "explore_deseq_dataset_obj", "PCA_rlog_transformed.jpg"),
        os.path.join("output", "kallisto", "{condition}", "explore_deseq_dataset_obj", "PCA_glmPCA.jpg")
    log:
        "logs/kallisto/exploration_{condition}.log",
    params:
        output_file_paths = lambda wildcards, output: output[1:]
    conda:
        "../envs/deseq2_gene_level_env.yaml"
    threads:
        12
    resources:
        mem_mb=12*4096,         # total: assign 4096MB per CPU
        cpus=12,	            # uses 12 cpus
        time_min=360		    # use 6 hrs, to make sure
    script:
        "../scripts/deseq2/explore_deseq_dataset.R"

