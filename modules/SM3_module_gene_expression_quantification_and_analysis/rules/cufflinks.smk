# https://www.biostars.org/p/160808/   <- Cufflinks workflow advices
# http://compbio.mit.edu/cummeRbund/manual_2_0.html <- Use CummeRbund for result visualizations

# ToDo: Check this: http://homer.ucsd.edu/homer/basicTutorial/rnaseqCufflinks.html -> To infer isoforms: cufflinks -o ...

import os

# Output directory
cufflinks_output_path = os.path.join('output', config["output_directories"]["cufflinks_output_dir"])


rule cufflinks_prepare_bam_files:
    """
    Cufflinks is not able to handle soft-clipped reads, that span over the ende of contigs/chromosomes.
    Hence we need to remove them from the bam files.
    Reference: https://groups.google.com/g/rna-star/c/Ta1Z2u4bPfc
    """
    input:
        # Returns wildcard: '{sample_id}'
        main_helper_get_bam_file_for_sample(config["input_dir_of_bam_files"],
            filename_extension=config["bam_files_attributes"]["filename_extension"])
    output:
        temp(os.path.join(cufflinks_output_path, "cufflinks_formated_sam_files", '{sample_id}.noS.sam'))
    conda:
        "../envs/samtools_env.yaml"
    shell:
        "samtools view -h {input} | "
        "awk 'BEGIN {{OFS=\"\t\"}} {{"
            "split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); "
                "if (C[2]==\"S\") {{$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}}; "
                "if (C[length(C)]==\"S\") {{L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }}; "
            "gsub(/[0-9]*S/,\"\",$6); print}}' - >{output}"


def get_cufflinks_library_type(sample_id):
    """
    Returns the library type for the given sample id.
    See details here:
    http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/#library-types
    """
    strandedness = pep.sample_table[pep.sample_table["sample_name"] == sample_id]["stranded"][0]
    if strandedness == "no":
        return "fr-unstranded"
    elif strandedness == "yes":
        return "fr-secondstrand"
    elif strandedness == "reverse":
        return "fr-firststrand"
    else:
        raise ValueError("Unknown strandedness: {}".format(strandedness))


# TODO gets often stuck...
rule cufflinks_transcriptome_assembly:
    """
    Assemble transcriptome per sample -> per BAM-file
    """
    input:
        rules.cufflinks_prepare_bam_files.output
    output:
        os.path.join(cufflinks_output_path, "transcriptome_assembly", "{sample_id}", "transcripts.gtf")
    params:
        output_dir=lambda wildcards, output: os.path.dirname(output[0]),
        gtf_file=config["cufflinks_settings"]["reference_genome_annotation_file"],
        ref_seq=config["cufflinks_settings"]["reference_genome_fasta_file"],
        min_isoform_fraction=config["cufflinks_settings"]["min_isoform_fraction"],
        min_frags_per_transfrag=config["cufflinks_settings"]["min_frags_per_transfrag"],
        library_type=lambda wildcards: get_cufflinks_library_type(wildcards.sample_id),
        extra_options=config["cufflinks_settings"]["extra_options_cufflinks"]
    threads: 32
    resources:
        mem_mb=32*4096,		# total: assign 4096MB per CPU = 131072 MB, 131GB
        cpus=32,			# uses 32 cpus
        time_min=24*60		# use 24 hrs, to make sure -> normally takes about 1hr 30 min [PolyA], but can take longer for RiboD (>4hrs...)
    conda:
        "../envs/cufflinks.yaml"
    shell:
        'cufflinks '
        '--num-threads {threads}'
        ' --library-type {params.library_type}'
        ' --GTF-guide {params.gtf_file}'                            # use reference transcript annotation to guide assembly, but also includes novel transcripts
        ' --frag-bias-correct {params.ref_seq} '                    # use bias correction - reference fasta required 
        ' --min-isoform-fraction {params.min_isoform_fraction}'     # suppress transcripts below this abundance level (compared with major isoform of the gene)
        ' --min-frags-per-transfrag {params.min_frags_per_transfrag}'   # assembled transfrags supported by fewer than this many aligned RNA-Seq fragments are ignored
        ' --output-dir {params.output_dir} '
        '{params.extra_options}'                                    # additional options
        '{input}'


rule cufflinks_compose_merge:
    """
    Create a text-file listing all GTF-files (created by cufflinks) that should be merged together by Cufflinks
    """
    input:
        # cufflinks_transcriptome_assembly
        expand(os.path.join(cufflinks_output_path, "transcriptome_assembly", "{sample_id}", "transcripts.gtf"),
            sample_id=pep.sample_table["sample_name"].tolist())
    output:
        all_transcriptome_assemblies_file=os.path.join(cufflinks_output_path, "transcriptome_assembly",
            "all_samples_transcriptome_assemblies.txt")
    threads: 1
    run:
        with open(output.all_transcriptome_assemblies_file, 'w') as out:
            print(*input, sep="\n", file=out)


rule cufflinks_merge_assemblies:
    """
    Merge all assembled transcriptome assemblies
    """
    input:
        rules.cufflinks_compose_merge.output
    output:
        os.path.join(cufflinks_output_path, "transcriptome_assembly", 'merged.gtf')
    params:
        output_dir=lambda wildcards, output: os.path.dirname(output[0]),
        gtf_file=config["cufflinks_settings"]["reference_genome_annotation_file"],
        ref_seq=config["cufflinks_settings"]["reference_genome_fasta_file"]
    threads: 32
    resources:
        mem_mb=131072,		# total: assign 4096MB per CPU = 131072 MB, 131GB
        cpus=32,			# uses 32 cpus
        time_min=2*60		# use 2 hrs -> takes 1 hr with 1 CPU
    conda:
        "../envs/cufflinks.yaml"
    shell:
        'cuffmerge -o {params.output_dir}'
        ' -g {params.gtf_file} '
        ' -s {params.ref_seq} '
        ' -p {threads} '
        '{input}'


def get_samples_for_given_condition(wildcards):
    """
    Returns for a given condition a list of samples (control & condition), that should be
    analyzed with cufflinks
    """
    condition = wildcards.condition
    # Get samples
    control_samples_ids = pep.sample_table[pep.sample_table["control"] == "true"]["sample_name"].tolist()
    condition_samples_ids = pep.sample_table[pep.sample_table["condition"] == condition]["sample_name"].tolist()

    # Get paths to BAM-files
    input_dir_choice = config["input_dir_of_bam_files"]
    filename_extension = config["bam_files_attributes"]["filename_extension"]
    control_samples_bam_files = main_helper_get_all_bam_file_paths(control_samples_ids,
        input_dir_choice, filename_extension=filename_extension)
    condition_samples_bam_files = main_helper_get_all_bam_file_paths(condition_samples_ids,
        input_dir_choice, filename_extension=filename_extension)

    return {
        "control_samples_bam_files": control_samples_bam_files,
        "condition_samples_bam_files": condition_samples_bam_files
    }


rule cufflinks_complete_analysis:
    """
    Differential expression analysis
    """
    # TODO: Samples are grouped by control and patient, but then summarized as replicates
    #   One could split it more up, but then cuffdiff would perform a pairwise comparison for each sample combination
    # ToDo: Process seem to stuck for a long time, before finally finishing
    input:
        unpack(get_samples_for_given_condition),
        merged_cufflinks_transcriptomes_gtf=rules.cufflinks_merge_assemblies.output,
    output:
        # 1. FPKM tracking files
        expand(os.path.join(cufflinks_output_path,"{{condition}}",'{level}.fpkm_tracking'),
                level=["isoforms", "genes", "cds", "tss_groups"]),
        # 2. Count tracking files
        expand(os.path.join(cufflinks_output_path,"{{condition}}",'{level}.count_tracking'),
                level=["isoforms", "genes", "cds", "tss_groups"]),
        # 3. Read group tracking files
        expand(os.path.join(cufflinks_output_path,"{{condition}}",'{level}.read_group_tracking'),
                level=["isoforms", "genes", "cds", "tss_groups"]),
        # 4. Differential expression tests
        expand(os.path.join(cufflinks_output_path,"{{condition}}",'{level}_exp.diff'),
                level=["isoform", "gene", "cds", "tss_groups"]),
        # 5. Differential splicing tests
        os.path.join(cufflinks_output_path,"{{condition}}",'splicing.diff'),
        # 6. Differential coding output (Only genes producing two or more distinct CDS -> multi-protein genes)
        # are listed here
        os.path.join(cufflinks_output_path,"{{condition}}",'cds.diff'),
        # 7. Differential promoter usage
        os.path.join(cufflinks_output_path,"{{condition}}",'promoters.diff'),
        # 8. Read group info
        os.path.join(cufflinks_output_path,"{{condition}}",'read_groups.info'),

    params:
        output_dir=lambda wildcards, output: os.path.dirname(output[0]),
        labels=lambda wildcards: ",".join(["control", wildcards.condition]),
        ref_seq=config["cufflinks_settings"]["reference_genome_fasta_file"],
        extra_options=config["cufflinks_settings"]["extra_options_cuffdiff"],

        # Entries connected by commas are treated as replicates
        ctrl_replicates=lambda wildcards, input:
            ",".join(input.control_samples_bam_files),
        con_replicates=lambda wildcards, input:
            ",".join(input.condition_samples_bam_files),
    threads:
        100
    resources:
        mem_mb=500000,		# total: assign 4096MB per CPU = ...
        cpus=100,			# uses 100 cpus
        time_min=60*24*3	# use 24 hrs (with these settings normally only 1:20 hrs:min are needed)
    conda:
        "../envs/cufflinks.yaml"
    shell:
        'cuffdiff '
        '--num-threads {threads} '
        '--output-dir {params.output_dir} '
        '--labels {params.labels} '
        '--frag-bias-correct {params.ref_seq} '
        '{params.extra_options} '
        '{input.merged_cufflinks_transcriptomes_gtf} '
        '{params.ctrl_replicates} '
        '{params.con_replicates}'	# samples are separated by space, replicates are separated by comma


rule cufflinks_create_html_reports:
    input:
        # 4. Differential expression tests
        expand(os.path.join(cufflinks_output_path,"{{condition}}",'{level}_exp.diff'),
            level=["isoform", "gene", "cds", "tss_groups"]),
        # 5. Differential splicing tests
        os.path.join(cufflinks_output_path,"{{condition}}",'splicing.diff'),
        # 6. Differential coding output (Only genes producing two or more distinct CDS -> multi-protein genes)
        # are listed here
        os.path.join(cufflinks_output_path,"{{condition}}",'cds.diff'),
        # 7. Differential promoter usage
        os.path.join(cufflinks_output_path,"{{condition}}",'promoters.diff'),
    output:
        report(
            directory(os.path.join(cufflinks_output_path, "html_reports", "{condition}")),
            patterns=["{filename}.report.html"],
            # Path relative to Snakefile
            caption="../../../report_source/module3_gene_expression_quantification_and_analysis/cufflinks/caption.rst",
            category=config["output_dir_module3_gene_expression_quantification_and_analysis"].replace("output/",""),
            subcategory="Cufflinks",
            labels={  # Set labels manually
                "Indication": "{condition}",
                "File:": "{filename}"
            }
        ),
        # 4. Differential expression tests
        expand(os.path.join(cufflinks_output_path, "html_reports", "{{condition}}",'{level}_exp.diff.report.html'),
            level=["isoform", "gene", "cds", "tss_groups"]),
        # 5. Differential splicing tests
        os.path.join(cufflinks_output_path, "html_reports", "{{condition}}",'splicing.diff.report.html'),
        # 6. Differential coding output (Only genes producing two or more distinct CDS -> multi-protein genes)
        # are listed here
        os.path.join(cufflinks_output_path, "html_reports", "{{condition}}",'cds.diff.report.html'),
        # 7. Differential promoter usage
        os.path.join(cufflinks_output_path, "html_reports", "{{condition}}",'promoters.diff.report.html'),
    params:
        input_files=lambda wildcards,input: input,
        data_separators=["\t" for i in range(4)],
        data_titles=["Title needs to be found" for i in range(4)],
        info_texts=["Info-text needs to be found" for i in range(4)],
        html_output_dir=lambda wildcards, output: output[0],
        html_output_file_basenames=lambda wildcards, input: [
            (os.path.basename(input_file) + ".report.html") for input_file in input
        ]
    conda:
        "../../../report_source/report_env.yaml"
    script:
        "../../../scripts/create_report_html_files.R"

