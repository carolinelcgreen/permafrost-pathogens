# README.md #
# permafrost-pathogens

Effects of a thaw gradient on the abundance of virulence-like factors in Alaskan permafrost samples

Caroline Green

### Structure ###
* contrast-files directory contains .csv files comparing thawing/frozen and thawed/frozen at each location
* L03_Virulence.tsv/.csv are the count matrix files for SEED subsystem "Virulence" and "Virulence, Disease and Defense" for the shotgun metagenomic data
* count_information.csv is the metadata for the count matrix
* framework.R is the bioinformatics pipeline for the DESeq2 analysis
* relevantPathogens.R is a list of clinically-relevant human pathogens found in the full list of identified bacterial species
