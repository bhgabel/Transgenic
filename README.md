# Transgenic Paper - Haricharan Lab

## TCGA

Data obtained from cbioportal.org, TCGA July 2025 release

Gene expression is fpkm z-scores

- MLH1 low cutoff is bottom 12th percentile of gene expression
- MSH2 low cutoff is bottom 8th percentile of gene expression

MSI predictor from MANTIS - <https://ascopubs.org/doi/10.1200/PO.17.00073>

- MSI high cutoff is a MANTIS score of 0.4 as per author
- MSI low cutoff is top 5th percentile of scores (including MSI-H)