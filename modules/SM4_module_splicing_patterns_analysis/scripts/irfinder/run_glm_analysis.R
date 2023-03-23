library(DESeq2)
source(snakemake@input[["deseq2_constructor_script"]])  				# Load IRFinder-related function

results = read.table(snakemake@input[["irfinder_results_file_paths_collection"]])
paths = as.vector(results$V1)                                            # File names must be saved in a vector
experiment = read.table(snakemake@input[["sample_condition_mapping_file"]], header=T)
experiment$Condition=factor(experiment$Condition,levels=c(snakemake@wildcards[["condition"]], "None"))    # Set WT as the baseline in the analysis
rownames(experiment)=NULL                                                # Force removing rownames

# WARNING: make sure the rownames of `experiment` is set to NULL.
# WARNING: users MUST check if the order of files in the `path` matches the order of samples in `experiment` before continue

metaList=DESeqDataSetFromIRFinder(filePaths=paths, designMatrix=experiment, designFormula=~1)
# The above line generates a meta list containing four slots
# First slot is a DESeq2 Object that can be directly passed to DESeq2 analysis.
# Second slot is a matrix for trimmed means of intron depth
# Third slot is a matrix for correcting splicing depth flanking introns
# Fourth slot is a matrix for maximum splicing reads at either ends of introns
# We build a “null” regression model on the interception only.
# A “real” model can be assigned either here directly, or in the downstream step. See below

dds = metaList$DESeq2Object                       # Extract DESeq2 Object with normalization factors ready
print("Check design of matrix")
colData(dds)                                      # Check design of matrix


# Please note that sample size has been doubled and one additional column "IRFinder" has been added.
# This is because IRFinder considers each sample has two sets of counts: one for reads inside intronic region
# and one for reads at splice site, indicating by "IR" and "Splice" respectively.
# "IRFinder" is considered as an additional variable in the GLM model.
# Please also be aware that size factors have been set to 1 for all samples. Re-estimation of size factors is NOT recommended and is going to bias the result.
# More details at the end of the instruction.

design(dds) = ~Condition + Condition:IRFinder     # Build a formula of GLM. Read below for more details.
dds = DESeq(dds)                                  # Estimate parameters and fit to model

print("Check actual variable names assigned by DeSeq2")
resultsNames(dds)                                 # Check the actual variable name assigned by DESeq2


res.WT = results(dds, name = "ConditionWT.IRFinderIR")
# This tests if the number of IR reads are significantly different from normal spliced reads, in the WT samples.
# We might only be interested in the "log2FoldChange" column, instead of the significance.
# This is because "log2FoldChange" represents log2(number of intronic reads/number of normal spliced reads).
# So we have the value of (intronic reads/normal spliced reads) by

WT.IR_vs_Splice=2^res.WT$log2FoldChange

# As IR ratio is calculated as (intronic reads/(intronic reads+normal spliced reads))
# We can easily convert the above value to IR ratio by

IRratio.WT = WT.IR_vs_Splice/(1+WT.IR_vs_Splice)

# Similarly, we can get IR ratio in the KO samples
res.KO = results(dds, name = "ConditionKO.IRFinderIR")
KO.IR_vs_Splice=2^res.KO$log2FoldChange
IRratio.KO = KO.IR_vs_Splice/(1+KO.IR_vs_Splice)

# Finally we can test the difference of (intronic reads/normal spliced reads) ratio between WT and KO
res.diff = results(dds, contrast=list("ConditionKO.IRFinderIR","ConditionWT.IRFinderIR"))
write.csv(dxr1, file=snakemake@output[["output_results_csv_file"]], row.names=TRUE)

# We can plot the changes of IR ratio with p values
# In this example we defined significant IR changes as
# 1) IR changes no less than 10% (both direction) and
# 2) with adjusted p values less than 0.05

IR.change = IRratio.KO - IRratio.WT
# Create plot and save it as JPEG
output_plot_file = snakemake@output[["output_plot_file"]]
jpeg(file=output_plot_file)
print(plot(IR.change,col=ifelse(res.diff$padj < 0.05 & abs(IR.change)>=0.1, "red", "black")))
dev.off()


# PLEASE NOTE
# You might see dots at the same horizontal level (same IR change) either marked as significant (red) or not (black)
# This is because the Wald test applied above is testing on fold changes instead of absolute changes
# For example, both changes from 0.01 to 0.11 and from 0.8 to 0.9 have absolute changes of 0.1.
# But the fold changes for them are 11 folds and 1.1 folds respectively, and lead to different p values