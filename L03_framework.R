# WORKFLOW FOR FUNCTIONAL DATA DESEQ2 ANALYSIS
# Created: 7/01/2020
# Author: C. Green
# Last edited: 7/24/2020

# Citation: M. I. Love, W. Huber, S. Anders: Moderated 
# estimation of fold change and dispersion for RNA-Seq data 
# with DESeq2. bioRxiv (2014). doi:10.1101/002832

# download DESeq2 and BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")

BiocManager::install("DESeq2")
BiocManager::install("pheatmap")

install.packages("RColorBrewer")

# load libraries---------
library("DESeq2")
library("dplyr")
library("tidyr")
library("ggplot2")
library("RColorBrewer")


# DATA SET CREATION -------------------------------------
# using count matrix input data workflow

# read in count matrix
  # only including "Virulence" and "Virulence, Disease and Defense" SEED subsystems
  # from L03 sequence data. Control columns were already removed.
counts <- as.matrix(read.csv("~/Documents/CRREL/permafrost-pathogens/L03_Virulence.tsv", 
                                      sep = "\t", row.names=1))

# read in count matrix metadata
  # metadata table created based on DESeq2 requirements:
  # column 1 contains sample names in the order they appear in row one of count table
  # column 2 contains location (35, 45, 60, 83 meters or New Tunnel)
  # column 3 contains thaw.state (frozen (-3C), thawing (0C), thawed (6C))
  # column 4 contains combined_factor (both thaw.temp and location data)
count_information <- read.csv("~/Documents/CRREL/permafrost-pathogens/count_information.csv", row.names=1)

# format count matrix metadata
metadata <- count_information[,c("location","thaw.temp","combined_factor")]
metadata$location <- factor(metadata$location)
metadata$thaw.temp <- factor(metadata$thaw.temp)
metadata$combined_factor <- factor(metadata$combined_factor)


# verify correct labeling
all(rownames(metadata) %in% colnames(counts)) # TRUE
all(rownames(metadata) == colnames(counts)) # TRUE


# make DESeq Data Set from count matrix
  # design used combined_factor due to the comparison requirements
  # of the program
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ combined_factor)
dds

# ~~ re-level data ~~~
# location data leveled based on depth
dds$location <- factor(dds$location, levels = c("new_tunnel", 
                    "35_meters", "45_meters", "60_meters", "83_meters"))

# re-leveling thaw state over temp gradient (frozen to thawed)
dds$thaw.temp <- factor(dds$thaw.temp, levels = c("frozen", "thawing", "thawed"))

# ~~~ pre-filtering taken care of in BBTools ~~~



# DESeq ANALYSIS ----------------------------------------

  # automatic independent filtering = TRUE
  # alpha = 0.1
dds <- DESeq(dds)  
res <- results(dds)



# prints the list of comparisons
resultsNames(dds)

# ~~review independent filtering~~~

# prints filtering results
metadata(res)$filterThreshold
# 51.70751% of genes were filtered out
# 0.72585 mean count of genes filtered out

# print cook's distance
assays(dds)[["cooks"]]

# plots cook's distance for each sample
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

# plots dispersion estimates
plotDispEsts(dds)

# plot frequency and removal of p-val [1,0]
use <- res$baseMean > metadata(res)$filterThreshold
h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")

barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("top", fill=rev(colori), legend=rev(names(colori)))




# CONTRAST RESULTS --------------------------------------

  # creates specific comparison of frozen->thawing and frozen -> at each location
  # output is a DESEqResults table for the contrast
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




# EXPORT to .csv FILE -----------------------------------
# output is a csv file containing
  # L03 subsystem name
  # Log2 Fold Change for frozen->thawing
  # Adjusted p-value for frozen->thawing
  # Log2 Fold Change for frozen->thawed
  # Adjusted p-value for frozen->thawed

res35 <- c(as.data.frame(res35thawed_frozen@rownames), 
           as.data.frame(res35thawing_frozen$log2FoldChange),
           as.data.frame(res35thawing_frozen$padj),
           as.data.frame(res35thawed_frozen$log2FoldChange), 
           as.data.frame(res35thawed_frozen$padj))
write.csv(res35, file="res35.csv")
head(res35)

res45 <- c(as.data.frame(res45thawed_frozen@rownames),
           as.data.frame(res45thawing_frozen$log2FoldChange),
           as.data.frame(res45thawing_frozen$padj),
           as.data.frame(res45thawed_frozen$log2FoldChange),
           as.data.frame(res45thawed_frozen$padj))
write.csv(res45, file="res45.csv")

res60 <- c(as.data.frame(res60thawed_frozen@rownames), 
           as.data.frame(res60thawing_frozen$log2FoldChange),
           as.data.frame(res60thawing_frozen$padj),
           as.data.frame(res60thawed_frozen$log2FoldChange), 
           as.data.frame(res60thawed_frozen$padj))
write.csv(res60, file="res60.csv")

res83 <- c(as.data.frame(res83thawed_frozen@rownames),
           as.data.frame(res83thawing_frozen$log2FoldChange),
           as.data.frame(res83thawing_frozen$padj),
           as.data.frame(res83thawed_frozen$log2FoldChange), 
           as.data.frame(res83thawed_frozen$padj))
write.csv(res83, file="res83.csv")

resNew <- c(as.data.frame(resNewThawed_frozen@rownames),
            as.data.frame(resNewThawing_frozen$log2FoldChange),
            as.data.frame(resNewThawing_frozen$padj),
            as.data.frame(resNewThawed_frozen$log2FoldChange), 
            as.data.frame(resNewThawed_frozen$padj))
write.csv(resNew, file="resNew.csv")

# output ALL site results to csv files
  # (includes all info, not just log2 and padj)
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




# NO FILTER TEST ----------------------------------------
  # outputs one contrast file for 35meters frozen->thawed WITHOUT 
  # independent filtering and creates heatmap for those results.

# make results file
res35noFilt <- results(dds, contrast=c("combined_factor", "thawed35_meters", "frozen35_meters"), 
                    
                          independentFiltering = FALSE)

res35noFilt_df <- c(as.data.frame(res35noFilt@rownames),
             as.data.frame(res35noFilt$log2FoldChange), 
             as.data.frame(res35noFilt$padj))
write.csv(res35noFilt_df, file="res35noFilt.csv")

# read in result csv files as a tibble
res35noFiltcsv <- readr::read_csv("res35noFilt.csv", col_names = TRUE, col_types= "icdd")
res35noFiltcsv <- select(res35noFiltcsv, -1) # delete first column of numbers

# re name columns
colnames(res35noFiltcsv) <- c("L03", "thawed", "thawed.padj")

# gather and remove log2 NA rows
df <- res35noFiltcsv %>%
  gather(key="thaw.state", value= "Log2", 
        thawed, -thawed.padj, 
         na.rm=TRUE)

meters35noFilt <- ggplot(df, aes(x=L03, y = thaw.state , fill=Log2)) +
  # tile with black contour
  geom_tile(colour="black") +
  # B&W theme, no grey background
  theme_bw() +
  # get rid of y axis titles
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  # Green color theme for `fill`
  scale_fill_distiller(palette="RdYlBu", direction=-1) +
  # since there is no legend, adding a title
  labs(title = "35 Meters: NO FILTER Log Fold Change from Frozen")
meters35noFilt
# save to file in figures directory
ggsave(
  "35NoFilt.png",
  plot = meters35noFilt,
  device = png(),
  path = "~/Documents/CRREL/permafrost-pathogens/permafrost-figures",
  scale = 1,
  width = 30,
  height = 10,
  units = "cm",
  dpi = 300)








# DATA VISUALIZATION ----------

# ~~ 35 METERS ~~~

# read in result csv files as a tibble
res35csv <- readr::read_csv("res35.csv", 
                            col_names = TRUE, 
                            col_types="icdddd")
res35csv <- select(res35csv, -1) # delete first column of numbers

# re name columns
colnames(res35csv) <- c("L03", "thawing", "thawing.padj", "thawed", "thawed.padj")

# gather and remove log2 NA rows
df <- res35csv %>%
  gather(key="thaw.state", value= "Log2", 
         thawing, thawed, -thawing.padj, -thawed.padj, 
         na.rm=TRUE)

meters35 <- ggplot(df, aes(x=L03, y = thaw.state , fill=Log2)) +
  # tile with black contour
  geom_tile(colour="black") +
  # B&W theme, no grey background
  theme_bw() +
  # get rid of y axis titles
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  # Green color theme for `fill`
  scale_fill_distiller(palette="RdYlBu", direction=-1) +
  # since there is no legend, adding a title
  labs(title = "35 Meters: Log Fold Change from Frozen")
meters35

# save to file in figures directory
ggsave(
  "35meters.png",
  plot = meters35,
  device = png(),
  path = "~/Documents/CRREL/permafrost-pathogens/permafrost-figures",
  scale = 1,
  width = 30,
  height = 10,
  units = "cm",
  dpi = 300)



# ~~ 45 METERS ~~~

# read in result csv files as a tibble
res45csv <- readr::read_csv("res45.csv", col_names = TRUE, col_types="icdddd")
res45csv <- select(res45csv, -1) # delete first column of numbers

# re name columns
colnames(res45csv) <- c("L03", "thawing", "thawing.padj", "thawed", "thawed.padj")

# gather and remove log2 NA rows
df <- res45csv %>%
  gather(key="thaw.state", value= "Log2", 
         thawing, thawed, -thawing.padj, -thawed.padj, 
         na.rm=TRUE)

meters45 <- ggplot(df, aes(x=L03, y = thaw.state , fill=Log2)) +
  # tile with black contour
  geom_tile(colour="black") +
  # B&W theme, no grey background
  theme_bw() +
  # get rid of y axis titles
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  # Green color theme for `fill`
  scale_fill_distiller(palette="RdYlBu", direction=-1) +
  # since there is no legend, adding a title
  labs(title = "45 Meters: Log Fold Change from Frozen")
meters45

# save to file in figures directory
ggsave(
  "45meters.png",
  plot = meters45,
  device = png(),
  path = "~/Documents/CRREL/permafrost-pathogens/permafrost-figures",
  scale = 1,
  width = 30,
  height = 10,
  units = "cm",
  dpi = 300)



# ~~ 60 METERS ~~~

# read in result csv files as a tibble
res60csv <- readr::read_csv("res60.csv", col_names = TRUE, col_types="icdddd")
res60csv <- select(res60csv, -1) # delete first column of numbers

# re-name columns
colnames(res60csv) <- c("L03", "thawing", "thawing.padj", "thawed", "thawed.padj")

# gather and remove log2 NA rows
df <- res60csv %>%
  gather(key="thaw.state", value= "Log2", 
         thawing, thawed, -thawing.padj, -thawed.padj, 
         na.rm=TRUE)

meters60 <- ggplot(df, aes(x=L03, y = thaw.state , fill=Log2)) +
  # tile with black contour
  geom_tile(colour="black") +
  # B&W theme, no grey background
  theme_bw() +
  # get rid of y axis titles
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  # Green color theme for `fill`
  scale_fill_distiller(palette="RdYlBu", direction=-1) +
  # since there is no legend, adding a title
  labs(title = "60 Meters: Log Fold Change from Frozen")
meters60

# save to file in figures directory
ggsave(
  "60meters.png",
  plot = meters60,
  device = png(),
  path = "~/Documents/CRREL/permafrost-pathogens/permafrost-figures",
  scale = 1,
  width = 30,
  height = 10,
  units = "cm",
  dpi = 300)


# ~~ 83 METERS ~~~

# read in result csv files as a tibble
res83csv <- readr::read_csv("res83.csv", col_names = TRUE, col_types="icdddd")
res83csv <- select(res83csv, -1) # delete first column of numbers

# re name columns
colnames(res83csv) <- c("L03", "thawing", "thawing.padj", "thawed", "thawed.padj")

# gather and remove log2 NA rows
df <- res83csv %>%
  gather(key="thaw.state", value= "Log2", 
         thawing, thawed, -thawing.padj, -thawed.padj, 
         na.rm=TRUE)

meters83 <- ggplot(df, aes(x=L03, y = thaw.state , fill=Log2)) +
  # tile with black contour
  geom_tile(colour="black") +
  # B&W theme, no grey background
  theme_bw() +
  # get rid of y axis titles
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  # Green color theme for `fill`
  scale_fill_distiller(palette="RdYlBu", direction=-1) +
  # since there is no legend, adding a title
  labs(title = "83 Meters: Log Fold Change from Frozen")
meters83

# save to file in figures directory
ggsave(
  "83meters.png",
  plot = meters83,
  device = png(),
  path = "~/Documents/CRREL/permafrost-pathogens/permafrost-figures",
  scale = 1,
  width = 30,
  height = 10,
  units = "cm",
  dpi = 300)


# ~~ NEW TUNNEL ~~~

# read in result csv files as a tibble
resNewcsv <- readr::read_csv("resNew.csv", col_names = TRUE, col_types="icdddd")
resNewcsv <- select(resNewcsv, -1) # delete first column of numbers

# re name columns
colnames(resNewcsv) <- c("L03", "thawing", "thawing.padj", "thawed", "thawed.padj")

# gather and remove log2 NA rows
df <- resNewcsv %>%
  gather(key="thaw.state", value= "Log2", 
         thawing, thawed, -thawing.padj, -thawed.padj, 
         na.rm=TRUE)

metersNew <- ggplot(df, aes(x=L03, y = thaw.state , fill=Log2)) +
  # tile with black contour
  geom_tile(colour="black") +
  # B&W theme, no grey background
  theme_bw() +
  # get rid of y axis titles
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  # Green color theme for `fill`
  scale_fill_distiller(palette="RdYlBu", direction=-1) +
  # since there is no legend, adding a title
  labs(title = "New Tunnel: Log Fold Change from Frozen")
metersNew

# save to file in figures directory
ggsave(
  "NewTunnel.png",
  plot = metersNew,
  device = png(),
  path = "~/Documents/CRREL/permafrost-pathogens/permafrost-figures",
  scale = 1,
  width = 30,
  height = 10,
  units = "cm",
  dpi = 300)






