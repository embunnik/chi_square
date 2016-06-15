# chi_square
R code script for performing the chi-square test and BH FDR correction on a dataframe (Bunnik et al., Genome Biology, 2016)

This script performs a chi-square test on each row of a dataframe with 4 columns. The four values of each row are re-arranged into a 2x2-table and then subjected to the chi-square test. The list of raw p-values is then adjusted using Benjamini-Hochberg false discovery rate correction.
