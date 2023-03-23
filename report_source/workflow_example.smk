# Global workflow description
report: "report_source/workflow_simple.rst"

rule create_final_report:
    output:
        report(
            directory("some_dir_being_created"),    # Mark output...
            patterns=["{name}.txt"],
            caption="report_source/caption_test.rst", # Path relative to Snakefile
            category="Snakefile Main Report")
    shell:
        "mkdir {output}; "
        "for i in 1 2 3; do echo $i > {output}/$i.txt; done"


rule create_final_report2:
    output:
        report(
            directory("some_dir_being_created2"),    # Mark output...
            patterns=["{name}.txt"],        # To dynamically include files
            caption="report_source/caption_test.rst", # Path relative to Snakefile
            category="Snakefile Main Report",
            subcategory="Subcategory"
        )
    shell:
        "mkdir {output}; "
        "for i in 1 2 3; do echo $i > {output}/$i.txt; done"


rule generate_html_hierarchy:
    output:
        report(
            directory("test3"),
            caption="report_source/caption_test.rst",
            # htmlindex="/Users/oliverkuchler/Studium/Master_BioInfo_Berlin/Masterthesis/Topic_Splicing_Variants/splice-prediction/snakemake_workflows/Snakemake_Main/report_source/test_html/test.html",
            # htmlindex="test.html",
            patterns=["{name}.txt"],# To dynamically include files
            labels={                        # Set labels manually
                "model": "MODLEY",
                "figure": "some test text"
            }
        )
    shell:
        """
        # mimic writing of an HTML hierarchy
        mkdir {output}
        """