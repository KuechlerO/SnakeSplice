SnakeSplice
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The presented report was automatically created by
`SnakeSplice <https://github.com/KuechlerO/SnakeSplice/>`_
and displays the analysis results for the given input samples.


Report overview
----------------
The report is divided into the following sections:
	1. **Module: Quality Control, Preprocessing, and Alignment**
		RNAseq data often contains a large amount of low quality reads.
		To avoid these reads from being included in the downstream analysis,
		quality control and preprocessing steps are performed.
		Afterwards, the reads are aligned to the reference genome.

	2. **Module: Detection of possible Gene Fusion Events**
		Gene fusions can be a common phenomenon in e.g. cancer cells.
		To detect possible gene fusions, the Arriba tool is used.

	3. **Module: Gene expression quantification and analysis**
		In order to quantify the expression of genes and transcripts,
		this modules employs the Salmon, and Kallisto tools.
		Furthermore, downstream analysis of the quantification results
		is possible, including differential expression analysis,
		and gene set enrichment analysis.
		Additionally, an outlier analysis can be performed.

	4. **Module: Analysis of Splicing Patterns**
		The analysis of Splicing patterns is performed using the
		various tools, which overall results are presented in a
		final summary table.
		Analysis tools include: Leafcutter, rMATs, FRASER, DexSeq,
		and PJD (Private Junction Detection).