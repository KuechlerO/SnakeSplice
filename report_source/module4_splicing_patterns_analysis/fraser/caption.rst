`FRASER <http://www.bioconductor.org/packages/release/bioc/html/FRASER.html>`_: Detection of rare aberrant splicing events in transcriptome profiles.
Read count ratio expectations are modeled by an autoencoder to control for confounding factors in the data.
Given these expectations, the ratios are assumed to follow a beta-binomial distribution with a junction specific dispersion.
Outlier events are then identified as read-count ratios that deviate significantly from this distribution.
FRASER is able to detect alternative splicing, but also intron retention.
The package aims to support diagnostics in the field of rare diseases where RNA-seq is performed to identify aberrant splicing defects.