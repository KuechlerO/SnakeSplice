`OUTRIDER <https://bioconductor.org/packages/release/bioc/html/OUTRIDER.html>`_:
Identification of aberrant gene expression in RNA-seq data.
Read count expectations are modeled by an autoencoder to control for confounders in the data.
Given these expectations, the RNA-seq read counts are assumed to follow a negative binomial distribution with a gene-specific dispersion.
Outliers are then identified as read counts that significantly deviate from this distribution.
Furthermore, OUTRIDER provides useful plotting functions to analyze and visualize the results.