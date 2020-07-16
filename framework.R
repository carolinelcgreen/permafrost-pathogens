
# WORKFLOW FOR FUNCTIONAL DATA DESEQ2 ANALYSIS
# Author: C. Green
# Last edited: 7/15/2020

# Citation: M. I. Love, W. Huber, S. Anders: Moderated 
# estimation of fold change and dispersion for RNA-Seq data 
# with DESeq2. bioRxiv (2014). doi:10.1101/002832

# download DESeq2 and BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")

BiocManager::install("DESeq2")

# ~~~ data set creation ~~~
# using count matrix input data workflow

# read in count matrix
  # only including "Virulence" and "Virulence, Disease and Defense" SEED subsystems
  # from L03 sequence data
counts <- as.matrix(read.csv("~/Documents/CRREL/permafrost-pathogens/L03_Virulence.tsv", 
                                      sep = "\t", row.names=1))

# read in count matrix metadata
count_information <- read.csv("~/Documents/CRREL/permafrost-pathogens/count_information.csv", row.names=1)

# format count matrix metadata
metadata <- count_information[,c("location","thaw.temp","combined_factor")]
metadata$location <- factor(metadata$location)
metadata$thaw.temp <- factor(metadata$thaw.temp)
metadata$combined_factor <- factor(metadata$combined_factor) 


# verify correct labeling
all(rownames(metadata) %in% colnames(counts)) # TRUE
all(rownames(metadata) == colnames(counts)) #TRUE

# open DESeq2
library("DESeq2")

install.packages("DelayedArray")

# make DESeq Data Set from count matrix
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ combined_factor)
dds

# ~~ re-level data ~~~
# re-level location data based on control group, 
  # then depth of sample
dds$location <- factor(dds$location, levels = c("new_tunnel", 
                    "35_meters", "45_meters", "60_meters", "83_meters"))

# re-leveling thaw state over temp gradient (frozen to thawed)
dds$thaw.temp <- factor(dds$thaw.temp, levels = c("frozen", "thawing", "thawed"))

# ~~~ pre-filtering taken care of in BBTools ~~~

# ~~collapse technical replicates~~
  # technical replicates are grouped into a single row
  # does having 3 v 4 replicates for some groups affect abundance?
?collapseReplicates
ddsColl <- collapseReplicates(dds, dds$combined_factor, renameCols=TRUE)

colData(ddsColl)
colnames(ddsColl)

# check that the sum of the counts for "replicate1" is the same
# as the counts in the "replicate1" column in ddsColl
matchFirstLevel <- dds$combined_factor == levels(dds$combined_factor)[1]
stopifnot(all(rowSums(counts(dds[,matchFirstLevel])) == counts(ddsColl[,1])))



# ~~~differential expression analysis~~~
  # what are the default tests being run?
  # what are the comparisons we are looking for?
?DESeq
dds <- DESeq(dds)  
res <- results(dds)
res

# prints the list of comparisons
resultsNames(dds)

# ~~~ contrast results ~~~
# do we use?? independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=TRUE
# do we want to LFCShrink???
res35thawed_frozen <- results(dds, contrast=c("combined_factor", "thawed35_meters", "frozen35_meters"))
res35thawing_frozen <- results(dds, contrast=c("combined_factor", "thawing35_meters", "frozen35_meters"))

res45thawed_frozen <- results(dds, contrast=c("combined_factor", "thawed45_meters", "frozen45_meters"))
res45thawing_frozen <- results(dds, contrast=c("combined_factor", "thawing45_meters", "frozen45_meters"))

res60thawed_frozen <- results(dds, contrast=c("combined_factor", "thawed60_meters", "frozen60_meters"))
res60thawing_frozen <- results(dds, contrast=c("combined_factor", "thawing60_meters", "frozen60_meters"))

res83thawed_frozen <- results(dds, contrast=c("combined_factor", "thawed83_meters", "frozen83_meters"))
res83thawing_frozen <- results(dds, contrast=c("combined_factor", "thawing83_meters", "frozen83_meters"))

resNewThawed_frozen <- results(dds, contrast=c("combined_factor", "thawednew_tunnel", "frozennew_tunnel"))
resNewThawing_frozen <- results(dds, contrast=c("combined_factor", "thawingnew_tunnel", "frozennew_tunnel"))

?lfcShrink

# output site results to csv files
write.csv(as.data.frame(res35thawed_frozen), 
          file="res35thawed_frozen.csv")
write.csv(as.data.frame(res35thawing_frozen), 
          file="res35thawing_frozen.csv")

write.csv(as.data.frame(res45thawed_frozen), 
          file="res45thawed_frozen.csv")
write.csv(as.data.frame(res45thawing_frozen), 
          file="res45thawing_frozen.csv")

write.csv(as.data.frame(res60thawed_frozen), 
          file="res60thawed_frozen.csv")
write.csv(as.data.frame(res60thawing_frozen), 
          file="res60thawing_frozen.csv")

write.csv(as.data.frame(res83thawed_frozen), 
          file="res835thawed_frozen.csv")
write.csv(as.data.frame(res83thawing_frozen), 
          file="res83thawing_frozen.csv")

write.csv(as.data.frame(resNewThawed_frozen), 
          file="resNewThawed_frozen.csv")
write.csv(as.data.frame(resNewThawing_frozen), 
          file="resNewThawing_frozen.csv")



