# Abstract
This workflow (module) assesses the presence of gene fusions in the RNA-seq data. 


# Workflow specific requirements
## Arriba
### Blacklist file - `arriba_blacklist_file` 
Arriba allows to use a blacklist file to filter out known false positive fusions. 
However, this is only a recommendation, and not strictly enforced. 
For more information, visit the [Arriba documentation](https://arriba.readthedocs.io/en/latest/input-files/#blacklist).


# Workflow Details
## Input - `input_dir_of_bam_files` 
- Aligned and sorted BAM files

## Output
- Arriba: Gene Fusion Detection - `output/arriba/`:\
Results from Arriba are saved here.


## Settings
The respective module settings can be adjusted via the `config_files/config_module2_gene_fusion.yaml`-file.
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