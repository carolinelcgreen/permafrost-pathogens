#load library

#read in matrix count data as a tsv
sampleCounts <- as.matrix(read.csv("~/Documents/CRREL/permafrost-pathogens/sample_counts.tsv", 
                             sep = "\t", row.names=1))

# read in matrix count data
sampleMeta <- read.csv("~/Documents/CRREL/permafrost-pathogens/sample_data.csv", row.names=1)

# format collection data
sampleMetadata <- sampleMeta[,c("location","thaw.temp")]
sampleMetadata$location <- factor(sampleMetadata$location)
sampleMetadata$thaw.temp <- factor(sampleMetadata$thaw.temp)
#nsampleMetadata$replicate <- factor(sampleMetadata$replicate)

head(sampleCounts,2)

# verify correct labeling
all(rownames(sampleMetadata) %in% colnames(sampleCounts)) # TRUE
all(rownames(sampleMetadata) == colnames(sampleCounts)) #TRUE

# open DESeq2
library("DESeq2")

# make DESeq Data Set from count matrix
sample_dds <- DESeqDataSetFromMatrix(countData = sampleCounts,
                              colData = sampleMetadata,
                              design = ~ location)

sample_dds

# re-leveling location data based on control group, 
# then age of sample (depth)
sample_dds$location <- factor(sample_dds$location, levels = c("control",
                                                "new tunnel", "35 meters", "45 meters"))

# re-leveling thaw state over temp gradient (frozen to thawed)
sample_dds$thaw.temp <- factor(sample_dds$thaw.temp, levels = c("frozen", "thawing", "thawed"))

# differential expression analysis
sample_dds <- DESeq(sample_dds)
sample_res <- results(sample_dds)
sample_res

# COLLAPSE REPLICATES
dds <- makeExampleDESeqDataSet(m=12)

# make data with two technical replicates for three samples
dds$sample <- factor(sample(paste0("sample",rep(1:9, c(2,1,1,2,1,1,2,1,1)))))
dds$run <- paste0("run",1:12)

ddsColl <- collapseReplicates(dds, dds$sample, dds$run)

# examine the colData and column names of the collapsed data
colData(ddsColl)
colnames(ddsColl)

# check that the sum of the counts for "sample1" is the same
# as the counts in the "sample1" column in ddsColl
matchFirstLevel <- dds$sample == levels(dds$sample)[1]
stopifnot(all(rowSums(counts(dds[,matchFirstLevel])) == counts(ddsColl[,1])))






