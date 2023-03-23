# Abstract
Preprocessing and assessment steps are essential before further downstream processing.
Thus, this workflow (module) covers the quality control part, but includes also preprocessing 
(trimming) and the subsequent alignment of the input reads against a reference genome.

# Workflow specific requirements
## Alignment tools
This workflow offers the usage of two different alignment algorithms: \
[STAR](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html) & [OLego](https://zhanglab.c2b2.columbia.edu/index.php/OLego) \
The automatic installation of the aligner of your choice is handled by this workflow.
However, the Olego aligner can also be installed manually and its path to the installation directory
can be specified in the config file.
The required steps for their installation are given here:

### Installation of OLego
Link: https://zhanglab.c2b2.columbia.edu/index.php/OLego_Documentation \
1. `cd {to_dir_where_you_want_to_save_all}`
2. `git clone https://github.com/chaolinzhanglab/olego.git`
3. `cd olego`
4. `make`
5. Set parent directory of the Olego installation directory in the config file via `olego_installation_dir`.
6. Insert `olego_installed.done` into the module's output directory to indicate that the installation was successful
   (i.e. output/{module_output_dir}/output/olego_installed.done).


## Required data sets
Furthermore, in addition to the data provided by the Main-Snakemake config-file,
additional specific data-files are required as inputs. 
The automatic download of these files is also facilitated by this workflow.
-> TODO remove this part, since it is not implemented yet
However, if desired, the download can also be done manually and the path to the downloaded files
can be specified in the config file.

### 1. Adaptor sequences: Download Illumina adaptor sequences (in data-folder) to remove them in the trimming steps
Download illumina's adaptor sequences\
`curl -O https://gist.githubusercontent.com/photocyte/3edd9401d0b13476e60f8b104c2575f8/raw/7e1e7e3dac674c196d8a05169558b675a57b7e23/Sequencing_adaptors.fasta`

### 2. (Mini-)Kraken2 database: Kraken2 database for taxonomic classification of reads (check for contamination)
`curl ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/old/minikraken2_v1_8GB_201904.tgz` \
`tar -xvzf minikraken2_v1_8GB_201904.tgz`

### 3. ALFA: Relies on a GTF-file in ENSEMBL format (downloaded from UCSC)
Links:
- [GTF for GRCh37](http://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz)
- [GTF for GRCh38](http://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.107.gtf.gz)


# Workflow Details
## Input
- RNAseq paired-end reads in FASTQ-files (file.fq)


## Output
- Check strandedness of reads - `output/check_strandedness/`:\
Results for the strandedness check of each read-pair/sample are saved here.
- Trimmomatic - `output/trimmomatic/`:\
Result for the trimming of each read-pair/sample are saved here.
- Kraken2-reports in `"output/trimmomatic/{sample}_report_summary.report"`:\
Results of checks for potential contamination in the studied samples.
- FASTQC - `output/fastqc_original/` & `output/fastqc_after_trimming/`:\
Trimmed reads are saved here
- STAR alignment - `output/star/{sample}/{sample}.sorted.bam`:\
Aligned reads are saved here
- OLego alignment - `output/olego/{sample}/{sample}.sorted.bam`:\
Aligned reads are saved here
- QualiMap - `output/qualimap2/star/{sample}_qualimapReport.html` or `output/qualimap2/olego/{sample}_qualimapReport.html`:\
Quality assessment of alignment results is saved here.
- ALFA - `output/alfa/{sample}_alfa_report.txt`:\
Creates an assessment of the feature distribution in the aligned reads.
- MultiQC-report in `output/multiQC/multiqc.html`:\
Summary HTML file of all Quality Assessment results.


## Settings
The respective module settings can be adjusted via the `config_files/config_module1_qc_preproc_alignment.yaml`-file.
In the configuration file you can find 2 different sections: 
1. `switch_variables`: Here you can switch on/off the individual functions of this module.
2. function specific variables: Here you can adjust the parameters of the individual functions of this module.
The following table lists the different variables with example entries.

For more information explore the `config_files/config_module1_qc_preproc_alignment.yaml`-file.
Every setting is decorated with a more detailed comment.

| Variable name                                                  | Example entry               | Explanation                                                                                                                            |
|----------------------------------------------------------------|-----------------------------|----------------------------------------------------------------------------------------------------------------------------------------|
| trimming_settings/file_sequencing_adaptors:sequencing_adaptors | "sequencing_adaptors.fasta" | Fasta-file of adaptors for trimmomatic                                                                                                 |
| contamination_check_settings/file_kraken2_db                   | "/data/kraken2_db"          | Path leading to required Kraken-database                                                                                               |
| olego_alignment_settings/olego_installation_dir                | "lib"                       | Directory where Olego code is saved.                                                                                                   |
| olego_alignment_settings/olego_regression_model                | "/olego/models/hg.cfg"      | Directory where Olego's precomputed regression model is saved. If empty: Standard model in olego-directory will be used: models/hg.cfg |
| olego_alignment_settings/olego_allowed_missmatches             | 4                           | Allowed nr of mismatches for alignemnt                                                                                                 |
| olego_alignment_settings/olego_min_exon_size                   | 6                           | Minimal length of exon (read-fragments)                                                                                                |
| alfa_settings/ensemble_sourced_gtf_file                        | "/path/to/ensembl.gtf"      | Path to GTF-file that is provided by Ensembl (ALFA relies on that specific format)                                                     |
| alfa_settings/require_gtf_conversion                           | True                        | If True, the chromosome names in the GTF-file are converted from Ensembl to non-Ensembl format                                         |
| qualimap_settings/qualimap_java_heap_size                      | "131G"                      | Java heap size for qualimap2                                                                                                           |


## Run workflow via Snakemake 
See main Snakemake workflow for details concerning the execution of the overall workflow.
Single modules can be activated or deactivated via the `config_files/config_main.yaml` file.