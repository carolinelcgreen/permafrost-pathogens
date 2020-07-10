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






