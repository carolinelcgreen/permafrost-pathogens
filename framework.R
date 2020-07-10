# Data processing framework
# a working outline of steps that need to be completed to
# statistically analyze the permafrost data

# WORKFLOW FOR FUNCTIONAL DATA (ie virulence factors)
# decide on relevant VFs and sort data

# download DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")

BiocManager::install("DESeq2")

#read in matrix count data as a tsv
counts <- as.matrix(read.csv("~/Documents/CRREL/permafrost-pathogens/L03_Virulence.tsv", 
                                      sep = "\t", row.names=1))

# read in matrix count data
sample_information <- read.csv("~/Documents/CRREL/permafrost-pathogens/sample_information.csv", row.names=1)

# format collection data
metadata <- sample_information[,c("location","thaw.temp","replicate")]
metadata$location <- factor(metadata$location)
metadata$thaw.temp <- factor(metadata$thaw.temp)
metadata$replicate <- factor(metadata$replicate)

head(counts,2)

# verify correct labeling
all(rownames(metadata) %in% colnames(counts)) # TRUE
all(rownames(metadata) == colnames(counts)) #TRUE

# open DESeq2
library("DESeq2")

# make DESeq Data Set from count matrix
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ location)
dds

# re-leveling location data based on control group, 
  # then age of sample (depth)
dds$location <- factor(dds$location, levels = c("control",
                      "new tunnel", "35 meters", "45 meters", 
                      "60 meters", "83 meters"))

# re-leveling thaw state over temp gradient (frozen to thawed)
dds$thaw.temp <- factor(dds$thaw.temp, levels = c("frozen", "thawing", "thawed"))

# can use dds$condition <- droplevels(dds$condition) to drop 
  # unneeded samples/metadata

# should we pre-filter?

# differential expression analysis

# save results to file

# output graphs






