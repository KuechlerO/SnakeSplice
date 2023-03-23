# Abstract
This workflow (module) performs a detection and quantification of transcripts on RNAseq data.

# Workflow specific requirements
None

TODO....
# Workflow Details
## Input - `input_dir_of_bam_files` 
- Aligned and sorted BAM files

## Output
- GATK-Workflow: Genomic - `output/gatk/calls/{sample}.g.vcf`:\
Results from the GATK-Workflow are saved here.
g.vcf files are VCF-like files, which contain a record for every position in the genome, even if there is no variant.\
Further information [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format)


## Settings
The respective module settings can be adjusted via the `config_files/config_module3_variant_calling.yaml`-file.
In the configuration file you can find 2 different sections: 
1. `switch_variables`: Here you can switch on/off the individual functions of this module.
2. function specific variables: Here you can adjust the parameters of the individual functions of this module.
The following table lists the different variables with example entries.

For more information explore the `config_files/config_module2_gene_fusion.yaml`-file.
Every setting is decorated with a more detailed comment.

| Variable name                             | Example entry        | Explanation                                                                |
|-------------------------------------------|----------------------|----------------------------------------------------------------------------|
| arriba_settings/arriba_optional_arguments | ""                   | Additional settings for the execution of arriba can be included here       |
| arriba_settings/arriba_blacklist_file     | "/path/to/blacklist" | Optional blacklist file for arriba (see previous sections for explanation) |


## Run workflow via Snakemake 
See main Snakemake workflow for details concerning the execution of the overall workflow.
Single modules can be activated or deactivated via the `config_files/config_main.yaml` file.