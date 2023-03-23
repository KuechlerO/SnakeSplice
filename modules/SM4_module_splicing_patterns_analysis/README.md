# Abstract
This workflow consists of two steps:
1. Identify events of differential splicing behaviour via [Leafcutter](https://davidaknowles.github.io/leafcutter/).
2. Visualize the Leafcutter output via [LeafViz](https://davidaknowles.github.io/leafcutter/articles/Visualization.html)

# Workflow specific requirements
This workflow utilizes the software tools Leafcutter in order to run its analysis steps.
Their respective installation steps are given below:


## Installation of Leafcutter
Links:\
https://davidaknowles.github.io/leafcutter/articles/Installation.html \
https://githubmemory.com/repo/davidaknowles/leafcutter/issues/188 \
1. `cd_to_dir_where_you_want_to_save_all`
2. `git clone https://github.com/davidaknowles/leafcutter`
3. `conda env create -f envs/original_leafcutter_env.yaml`
4. `conda activate original_leafcutter_env`
5. `R -e 'devtools::install_github("davidaknowles/leafcutter/leafcutter")`
6. __*Attention:*__ Since the devtools-step is needed for installing additional dependencies, this 4. step is also required when running a rule with its hidden conda-env (see leafcutter rules).


# Details

## Input
This workflow takes the previously created alignment results in form of BAM-files as input.
See `input_bams_dir` in setting section.

## Output

### 2. Leafcutter (see LeafViz for compact results)
Output-file: `output/leafcutter/<group>/diff_splicing_analysis/sample_analysis_cluster_significance.txt`\
Lists per cluster **(with matched gene)** its p-value & adjusted p-value -> significance in differential splicing.

Output-file: `output/leafcutter/all_perind_numers.counts.gz`\
This file lists all the clusters with their respective introns (rows) and how the given samples (columns) match with their number of supporting split reads.


### 3. LeafcutterMD
Output-file: `output/leafcutter/diff_splicing_analysis/outliers/filtered_all_introns_pvals.tsv`\
Output-file: `output/leafcutter/diff_splicing_analysis/outliers/filtered_all_clusters_pvals.tsv"`\
List for all given samples the UNADJUSTED! p-values. 

Output-file: `output/leafcutter/diff_splicing_analysis/outliers/filtered_bud13_intron_pvals.tsv`\
Output-file: `output/leafcutter/diff_splicing_analysis/outliers/filtered_supt7l_intron_pvals.tsv`\
List of filtered entries (p-value <= 0.05) for the respective assay.


### 4. LeafViz 
Output-file: `output/leafviz/{comparison_group}/leafcutter_ds_sig_clusters.tsv`\
Output-file: `output/leafviz/{comparison_group}/leafcutter_ds_sig_introns.tsv`\
Lists only significant clusters/introns for selected comparison.

Output-file: `output/leafviz/supt7l_analysis/supt7l_analysis.RData`\
Output-file: `output/leafviz/bud13_analysis/bud13_analysis.RData`\
Output-files for LeafViz-Visualization-Tool



## Settings
The respective workflow settings can be adjusted via the `config.yaml`-file. 
The following table lists the different variables with example entries.

Variable name  | Example entrey | Explanation
-------------- | -------------- | ------------
run_arriba | True | Switch for Arriba: Detect gene fusions from chimeric STAR output
run_leafcutter_diff_splicing_analysis | True | Switch for Leafcutter differential splicing analysis
run_leafcutter_leafviz_report | True | Switch for visualization of Leafcutter output via Leafviz. __Attention:__ Leafcutter output must be existant!!!
run_leafcutter_outlier_analysis | True | Outlier detection via Leafcutter
input_bams_dir | splice-prediction/Snakemake_QC_and_Alignment/output/star/ | Input directory of bam- and bai-files. __Attention:__ BAM-files need to be indexed -> .bai
RNA_seq_read_file_input_table | input_info/input_files_all.tsv | TSV-file with info for all samples
reference_seq | data/reference_genome/hg19_analysis_set/hg19.p13.plusMT.no_alt_analysis_set.fa | FASTA-file of reference genome
reference_gtf_file | data/reference_genome/annotation_files/hg19.ncbiRefSeq.gtf | Annotation file for reference genome
arriba_blacklist | arriba/arriba_v2.1.0/database/blacklist_hg19_hs37d5_GRCh37_v2.1.0.tsv.gz | Black list of genomic areas to be ignored
leafcutter_dir | /sc-projects/sc-proj-btg/olik_splicing_project/leafcutter | Directory where Leafcutter package is saved
leafcutter_output_files_prefix | controls_vs_supt7l | Prefix for Leafcutter's output files
leafcutter_group_files | input_info/leafcutter_groups_controls_vs_supt7l.tsv | TSV-file which annotates the group associations for each sample (needed for Leafcutter)
regtools_junctions_anchor_length | 8 | Anchor length that is utilized to detect splice junctions with regtools
regtools_junctions_mininimum_intron_length | 50 | Minimum intron length for detection of splice junctions via Regtools
regtools_junctions_maximum_intron_length | 500000 | Maximum intron length for detection of splice junctions via Regtools
leafcutter_clustering_supporting_split_reads | 1 | Minimum number of splitted reads supporting a cluster (Leafcutter)
leafcutter_clustering_maximum_intron_length | 500000 | Maximum intron length (Leafcutter)
leafcutter_min_samples_per_intron | 1 | Minimum number of samples needed to support on intron (Leafcutter)
leafcutter_min_samples_per_group | 1 | Require this many samples in each group to have at least min_coverage reads (Leafcutter)
leafcutter_min_coverage | 10 | Require min_samples_per_group samples in each group to have at least this many reads [default 20] (Leafcutter)
false_discovery_rate | 0.10 | Filter-Setting for significant results
