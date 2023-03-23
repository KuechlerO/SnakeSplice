import os

rmats_output_path = os.path.join('output',config["output_directories"]["rmats_output_dir"])


rule rmats_create_control_group_file:
    """
    Create a control group file for a given control group.
    """
    output:
        os.path.join(rmats_output_path,"control_group.txt")
    log:
        "logs/rmats_create_control_group_file.log"
    params:
        # # Get all BAM-file paths for the given controls. -> Comma-separated string.
        control_samples_bam_paths=
        ",".join(
            main_helper_get_all_bam_file_paths(
                pep.sample_table[pep.sample_table["control"] == "true"]["sample_name"],
                config["input_dir_of_bam_files"],
                config["bam_files_attributes"]["filename_extension"]),
        )
    threads:
        1
    shell:
        # Print all BAM-file paths to a file.
        """
        echo {params.control_samples_bam_paths} > {output} 2> {log}
        """


rule rmats_create_patient_group_file:
    """
    Create a patient group file for a given control group.
    """
    output:
        os.path.join(rmats_output_path,"{condition}_group.txt")
    log:
        "logs/rmats_create_{condition}_group_file.log"
    params:
        # Get all BAM-file paths for the given condition. -> Comma-separated string.
        patient_samples_bam_paths=lambda wildcards:
        ",".join(
            main_helper_get_all_bam_file_paths(
                pep.sample_table[pep.sample_table["condition"] == wildcards.condition]["sample_name"],
                config["input_dir_of_bam_files"],
                config["bam_files_attributes"]["filename_extension"],
            )
        )
    threads:
        1
    shell:
        # Print all BAM-file paths to a file.
        """
        echo {params.patient_samples_bam_paths} > {output} 2> {log}
        """


# Error -> Folder creation error if multiple instances are run at the same time...
rule rmats_run_prep:
    """
    Run the preparation step of rMATS.
    The input files are processed and a summary is saved to .rmats files in the --tmp directory.
    The .rmats files track info from each BAM separately according to the full path of the BAM 
    specified in the input .txt file. 
    """
    input:
        # BAM file -> returns with {sample_id}-wildcard
        bam_file=main_helper_get_bam_file_for_sample(
            config["input_dir_of_bam_files"],
            filename_extension=config["bam_files_attributes"]["filename_extension"]),
        gtf_file=config["rmats_settings"]["reference_genome_annotation_file"],
    output:
        rmats_tmp_file=os.path.join(rmats_output_path,"tmp","{sample_id}_done.txt")
    log:
        "logs/rmats_run_prep_{sample_id}.log"
    params:
        tmp_dir=lambda wildcards, output: os.path.dirname(output.rmats_tmp_file),
        output_dir_main=lambda wildcards, output: os.path.dirname(os.path.dirname(output.rmats_tmp_file)),

        paired_flag="-t paired" if config["bam_files_attributes"]["paired"] else "",
        libType=config["rmats_settings"]["library_type"],
        average_read_length=config["rmats_settings"]["average_read_length"],
        variable_read_length_flag="--variable-read-length" if config["rmats_settings"][
                                                                  "consider_variable_read_length"] == True else "",
        consider_novel_splice_sites_flag="--novelSS" if config["rmats_settings"]["consider_novel_splice_sites"] else "",
    threads:
        64
    resources:
        mem_mb=64 * 4096,# total: assign 4096MB per CPU = 131072
        cpus=64,# uses 8 cpus
        time_min=12*60  # give max 12 hrs
    conda:
        "../envs/rmats_env.yaml"
    shell:
        # Save path to BAM file in a text file
        "echo {input.bam_file} > {output.rmats_tmp_file} && "

        # Run rMATS pre
        "rmats.py "
        "--b1 {output.rmats_tmp_file} "  # input file
        "--gtf {input.gtf_file} "  # GTF file
        "{params.paired_flag} "  # Paired-end or single-end reads
        "--libType {params.libType} "  # Library type
        "--readLength {params.average_read_length} "  # Average read length
        "{params.variable_read_length_flag} "  # Variable read length
        "{params.consider_novel_splice_sites_flag} "  # Consider novel splice sites
        "--nthread {threads} "
        "--od {params.output_dir_main} "
        "--tmp {params.tmp_dir} "
        "--task prep "  # IMPORTANT: Prep Task
        "2> {log}"


rule rmats_run_post:
    """
    In the post step, .rmats files are read and the final output files are created.
    -> Statistical analysis
    """
    input:
        # Requires all samples to have gone through the prep step.
        rmats_tmp_files=expand(os.path.join(rmats_output_path,"tmp","{sample_id}_done.txt"),
            sample_id=pep.sample_table["sample_name"]),
        # Control group file
        group_file_control=rules.rmats_create_control_group_file.output,
        # Patient group file
        group_file_patient=rules.rmats_create_patient_group_file.output,
        gtf_file=config["rmats_settings"]["reference_genome_annotation_file"],
    output:
        os.path.join(rmats_output_path,"{condition}","summary.txt"),
        # Final output including only reads that span junctions defined by rmats (Junction Counts)
        expand(os.path.join(rmats_output_path,"{{condition}}", "{AS_event}.MATS.JC.txt"),
            AS_event=["A3SS", "A5SS", "SE", "MXE", "RI"]),
        # Final output including both reads that span junctions defined by rmats (Junction Counts) and
        # reads that do not cross an exon boundary (Exon Counts)
        expand(os.path.join(rmats_output_path,"{{condition}}","{AS_event}.MATS.JCEC.txt"),
               AS_event=["A3SS", "A5SS", "SE", "MXE", "RI"]),
        # All identified alternative splicing (AS) events derived from GTF and RNA
        expand(os.path.join(rmats_output_path,"{{condition}}","fromGTF.{AS_event}.txt"),
               AS_event=["A3SS", "A5SS", "SE", "MXE", "RI"]),
        # Alternative splicing (AS) events which were identified only after considering the
        # RNA (as opposed to analyzing the GTF in isolation). This does not include events
        # with an unannotated splice site.
        expand(os.path.join(rmats_output_path,"{{condition}}","fromGTF.novelJunction.{AS_event}.txt"),
           AS_event=["A3SS", "A5SS", "SE", "MXE", "RI"]),
        # This file contains only those events which include an unannotated splice site.
        expand(os.path.join(rmats_output_path,"{{condition}}","fromGTF.novelSpliceSite.{AS_event}.txt"),
           AS_event=["A3SS", "A5SS", "SE", "MXE", "RI"]),

        # ------ Count Files ---------
        # Event counts including only reads that span junctions defined by rmats (Junction Counts)
        expand(os.path.join(rmats_output_path,"{{condition}}","JC.raw.input.{AS_event}.txt"),
               AS_event=["A3SS", "A5SS", "SE", "MXE", "RI"]),
        # Event counts including both reads that span junctions defined by rmats (Junction Counts) and reads that do not cross an exon boundary (Exon Counts)
        expand(os.path.join(rmats_output_path,"{{condition}}","JCEC.raw.input.{AS_event}.txt"),
               AS_event=["A3SS", "A5SS", "SE", "MXE", "RI"]),
    log:
        "logs/rmats_run_post_{condition}.log"
    params:
        out_dir=lambda wildcards, output: os.path.dirname(output[0]),
        tmp_dir=lambda wildcards, input: os.path.dirname(input["rmats_tmp_files"][0]),

        paired_flag="-t paired" if config["bam_files_attributes"]["paired"] else "",
        libType=config["rmats_settings"]["library_type"],
        average_read_length=config["rmats_settings"]["average_read_length"],
        variable_read_length_flag="--variable-read-length" if config["rmats_settings"][
                                                                  "consider_variable_read_length"] == True else "",
        consider_novel_splice_sites_flag="--novelSS" if config["rmats_settings"]["consider_novel_splice_sites"] else "",
    threads:
        64
    resources:
        mem_mb=64 * 4096,# total: assign 4096MB per CPU = 131072
        cpus=64,# uses 8 cpus
        time_min=12 *60  # give max 12 hrs
    conda:
        "../envs/rmats_env.yaml"
    shell:
        # Run rMATS post
        "rmats.py "
        "--b1 {input.group_file_control} "  # input file
        "--b2 {input.group_file_patient} "  # input file
        "--gtf {input.gtf_file} "  # GTF file
        "{params.paired_flag} "  # Paired-end or single-end reads
        "--libType {params.libType} "  # Library type
        "--readLength {params.average_read_length} "  # Average read length
        "{params.variable_read_length_flag} "  # Variable read length
        "{params.consider_novel_splice_sites_flag} "  # Consider novel splice sites
        "--nthread {threads} "
        "--od {params.out_dir} "
        "--tmp {params.tmp_dir} "
        "--task post "  # IMPORTANT: Post Task
        # "--statoff "                                  # Do not create statistical analysis, only create count results
        "2> {log}"


rule rmats_filter_output:
    """
    Filter output files to only include significant events.
    """
    input:
        expand(os.path.join(rmats_output_path,"{{condition}}","{AS_event}.MATS.JC.txt"),
            AS_event=["A3SS", "A5SS", "SE", "MXE", "RI"]),
        # Final output including both reads that span junctions defined by rmats (Junction Counts) and
        # reads that do not cross an exon boundary (Exon Counts)
        expand(os.path.join(rmats_output_path,"{{condition}}","{AS_event}.MATS.JCEC.txt"),
            AS_event=["A3SS", "A5SS", "SE", "MXE", "RI"]),
    output:
        expand(os.path.join(rmats_output_path,"{{condition}}","filtered_{AS_event}.MATS.JC.txt"),
            AS_event=["A3SS", "A5SS", "SE", "MXE", "RI"]),
        expand(os.path.join(rmats_output_path,"{{condition}}","filtered_{AS_event}.MATS.JCEC.txt"),
            AS_event=["A3SS", "A5SS", "SE", "MXE", "RI"]),
    params:
        top_x_results=1000
    threads:
        1
    run:
        import pandas as pd

        for input_file, output_file in zip(input, output):
            df = pd.read_csv(input_file, sep="\t")
            # filter for significant results (adjusted p-value < 0.10)
            df = df[(df["PValue"] < 0.05) & (df["FDR"] < 0.10)]
            df = df.sort_values(by=["PValue"])

            # Extract top x results
            df = df.head(params.top_x_results)

            # write to file
            df.to_csv(output_file, sep="\t",index=False)


rule rmats_create_html_reports:
    input:
        # ------- Take filtered files as input -------
        os.path.join(rmats_output_path,"{condition}","summary.txt"),
        # Final output including only reads that span junctions defined by rmats (Junction Counts)
        # Exclude to save report space...
        # expand(os.path.join(rmats_output_path,"{{condition}}","filtered_{AS_event}.MATS.JC.txt"),
        #     AS_event=["A3SS", "A5SS", "SE", "MXE", "RI"]),
        # Final output including both reads that span junctions defined by rmats (Junction Counts) and
        # reads that do not cross an exon boundary (Exon Counts)
        expand(os.path.join(rmats_output_path,"{{condition}}","filtered_{AS_event}.MATS.JCEC.txt"),
            AS_event=["A3SS", "A5SS", "SE", "MXE", "RI"]),
    output:
        report(
            directory(os.path.join(rmats_output_path, "{condition}", "html_reports")),
            patterns=["filtered_{splice_event}.MATS.{counting}.txt.report.html"],
            # Path relative to Snakefile
            caption="../../../report_source/module4_splicing_patterns_analysis/rmats/caption.rst",
            category=config["output_dir_module4_splicing_patterns_analysis"].replace("output/",""),
            subcategory="rMATS",
            labels={  # Set labels manually
                "Condition": "{condition}",
                "Splice event": "{splice_event}",
                "Counting": "{counting}",
                "File:": "filtered_{splice_event}.MATS.{counting}"
            }
        ),
        # Just two output examples... -> Total of 26 files will be generated
        os.path.join(rmats_output_path, "{condition}", "html_reports", "summary.txt.report.html"),
        # expand(os.path.join(rmats_output_path,"{{condition}}","html_reports", "filtered_{AS_event}.MATS.JC.txt.report.html"),
        #        AS_event=["A3SS", "A5SS", "SE", "MXE", "RI"]),
        expand(os.path.join(rmats_output_path,"{{condition}}","html_reports", "filtered_{AS_event}.MATS.JCEC.txt.report.html"),
               AS_event=["A3SS", "A5SS", "SE", "MXE", "RI"]),
    params:
        input_files=lambda wildcards,input: input,
        data_separators=lambda wildcards, output: ["\t" for out_file in output[1:]],
        data_titles=lambda wildcards, output:
            [
                ("rMATS: Significant results for selection: "
                    + os.path.basename(out_file).replace(".txt.report.html", "")
                 + " , Condition: " + wildcards.condition
                ) for out_file in output[1:]
            ],
        info_texts=
            # Cannot use lambda function with workflow.source_path -> ToDo Buggy
            [workflow.source_path("../../../report_source/module4_splicing_patterns_analysis/rmats/info.html")
                for i in range(6)
            ],
        html_output_dir=lambda wildcards, output: output[0],
        html_output_file_basenames=lambda wildcards, output: [os.path.basename(out_file) for out_file in output[1:]]
    threads:
        1
    conda:
        "../../../report_source/report_env.yaml"
    script:
        "../../../scripts/create_report_html_files.R"
