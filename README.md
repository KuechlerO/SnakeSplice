# SnakeSplice
A Snakemake based modular Workflow that facilitates RNA-Seq analyses with a 
special focus on splicing


## Modular Build-up
The given parent workflow is a wrapper workflow, which includes the following sub-workflows (called modules): \
1. Module1: Quality Control, Preprocessing and Alignment
2. Module2: Gene Fusion Detection
3. Module3: Transcript Quantification & Expression Analysis
4. Module4: Splice Pattern Analysis


## Software Requirements
- Conda: [Conda Webpage](https://docs.conda.io/en/latest/miniconda.html)
- Snakemake: [Snakemake Webpage](https://snakemake.readthedocs.io/en/stable/index.html)
- For PEP required:
	1. peppy is required and can be installed via Conda:  `conda install -c conda-forge peppy`
	2. eido required is required and can be installed via Conda: `conda install -c conda-forge eido`

## Input Data
The input data for this workflow is provided via a sample sheet (default location: `input_data/input_samples.csv`), 
whereby the structure of the sample sheet is defined by the PEP (file `pep/pep_schema_config.yaml`) file.

### Tabular structure for the sample sheet
The sample sheet is a tabular file, which consists of the following columns:
1. `sample_name`: Name of the sample
2. `sample_directory`: Path to the directory, where the sample data (FASTQ-files) are located. This
information is only used if the FASTQ-files are needed.
3. `read1`: Name of the FASTQ-file for read1
4. `read2`: Name of the FASTQ-file for read2
5. `control`: true or false (if true, the sample is treated as control sample)
6. `condition`: Name of the condition (e.g. treatment group)
7. `protocol`: Name of the protocol (e.g. RNAseq-PolyA). This information is not yet used...
8. `stranded`: No, if library is unstranded, yes if library is stranded, reverse if library is reverse stranded
9. `adaptors_file`: Path to the file, which contains the adaptors for the sample
10. `additional_comment`: Additional comment for the sample

Currently, the entries for the columns `protocol` and `additional_comment` are not used.\
The entries "read1", read2" and "adaptors_file" are marked as mandatory, as they are needed for the
execution of the alignment workflow. However, if the user has already aligned the samples, these columns
can be either filled with dummy data (make sure the references files exist!), or one can manipulate the
PEP-file (path: `pep/pep_schema_config.yaml`) to make these columns optional.

## Configurations
The respective workflow settings can be adjusted via the configuration files, which
are placed in the directory `config_files`.
In this folder is a `config_main.yaml`-file, which holds the general settings for the 
workflow.
Additionally, every sub-workflow/module has its own 
`config_module{X}_{module_name}.yaml`-file, which lists the settings for the 
respective sub-workflow.


### Main Configuration File - `config_files/config_main.yaml`
This configuration file holds the general settings for this master workflow.
It consists of 2 parts:
1. Module switches - `module_swiches`:\
Here, the user can switch on/off the sub-workflows/modules, which should be executed.
**Note**: Submodule 1 has to be run first alone, as the output of this submodule is 
used as input for the other submodules. Subsequently, the other modules can be run in (almost) any order.
2. Module output directory names - `module_output_dir_names`:\
Every submodule saves their output in a separate sub-directory of the main output directory `output`.\
The names of these sub-directories can be adjusted here.


### Cluster Profile - `profile_cluster\config.yaml`
A profile configuration file can be used to summarize all desired settings for the snakemake execution, that
are normally passed to the snakemake command via the command line (with arguments).
This workflow offers a predefined profile configuration file for the execution on a cluster (using SLURM).
The respective setting options are listed and explained below.\
**Note**: Go to the bottom of this file to find out, how to execute Snakemake using this profile-settings file.


| Variable name     | Default entry                                                                                                                                                                    | Explanation                                                                          |
|-------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------|
| use-conda         | True                                                                                                                                                                             | Enables the use of conda environments (and Snakemake wrappers)                       |
| keep-going        | True                                                                                                                                                                             | Go on with independent jobs, if one job fails                                        |
| latency-wait      | 60                                                                                                                                                                               | Wait given seconds if an output file of a job is not present after the job finished. |
| rerun-incomplete  | True                                                                                                                                                                             | Rerun all jobs where the output is incomplete                                        |
| printshellcmds    | True                                                                                                                                                                             | Printout shell commands that will be executed                                        |
| jobs              | 50                                                                                                                                                                               | Number of jobs / rules to run (maximal) in parallel                                  |    
| default-resources | [cpus=1, mem_mb=2048, time_min=60]                                                                                                                                               | Default resources for each job (can be overwritten in the rule definition)           |
| resources         | [cpus=100, mem_mb=500000]                                                                                                                                                        | Resource constraints for the whole workflow                                          |
| cluster           | "sbatch -t {resources.time_min} --mem={resources.mem_mb} -c {resources.cpus} -o logs_slurm/{rule}.%j.out -e logs_slurm/{rule}.%j.out --mail-type=FAIL --mail-user=user@mail.com" | Cluster command for the execution on a cluster (here: SLURM)                         |


## Where to get data from? 
### Reference Genome
We recommend an analysis set reference genome. Its advantages over other common forms of reference genomes can
be read [here](https://gatk.broadinstitute.org/hc/en-us/articles/360041155232-Reference-Genome-Components).\
Download suitable FASTA-file of reference genome (e.g. analysis set reference genome for hg19):\
`curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.fa.gz`\
`gunzip hg19.p13.plusMT.no_alt_analysis_set.fa.gz`

### Gene Annotation File
Download gene annotation file (GTF-file) for reference genome (hg19) [no unzip needed]:\
1. `curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.ncbiRefSeq.gtf.gz`
2. `gunzip hg19.ncbiRefSeq.gtf.gz`


## Execute workflow via Snakemake

### Steps to simple execution of Snakemake workflow
1. Activate Conda-Snakemake environment\
`conda activate snakemake`
2. Execute Workflow (you can adjust the passed number of cores to your desire...) \
`snakemake -s Snakefile --cores 4 --use-conda`
3. Run Workflow in background\
`rm nohup.out && nohup snakemake -s Snakefile --cores 4 --use-conda &`

### Visualization & Dry Runs 
- Visualize DAG of jobs\
`snakemake --dag | dot -Tsvg > dag.svg`
- Dry run -> Get overview of job executions, but no real output is generated\
`snakemake -n -r --cores 4`

### Cluster: Execute Snakemake workflow on a HPC cluster
1. Adjust settings in profile-settings file (e.g. here in `profiles/profile_cluster/config.yaml`). 
2. Execute workflow\
mkdir -p logs_slurm && rm nohup.out || true && nohup snakemake --profile profiles/profile_cluster &

### Monitor execution stats on a HPC cluster with SLURM
`sacct -a --format=JobID,User,Group,Start,End,State,AllocNodes,NodeList,ReqMem,MaxVMSize,AllocCPUS,ReqCPUS,CPUTime,Elapsed,MaxRSS,ExitCode -j <job-ID>`
\
Explanation:
- `-a`: Show jobs for all users
- `--format=JobID...`: Format output

### Kill cluster jobs
`killall -TERM snakemake`

### Node stats on SLURM cluster
`sinfo -o "%n %e %m %a %c %C"`