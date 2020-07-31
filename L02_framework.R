# WORKFLOW FOR FUNCTIONAL DATA DESEQ2 ANALYSIS of L02 SUBSYSTEM DATA
  # a more informative view of LFC between frozen and thaw.state
# Created: 7/27/2020
# Author: C. Green
# Last edited: 7/27/2020

# load libraries
library("DESeq2")
library("dplyr")
library("tidyr")
library("ggplot2")
library("scales")
library("RColorBrewer")


# DATA SET CREATION -------------------------------------

# using count matrix input data workflow

# read in count matrix
# only including "Virulence" and "Virulence, Disease and Defense" SEED subsystems
# from L03 sequence data. Control columns were already removed.
L02_counts <- as.matrix(read.csv("~/Documents/CRREL/permafrost-pathogens/L02_Virulence.tsv", 
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
all(rownames(metadata) %in% colnames(L02_counts)) # TRUE
all(rownames(metadata) == colnames(L02_counts)) # TRUE

# open DESeq2
library("DESeq2")

# make DESeq Data Set from count matrix
# design used combined_factor due to the comparison requirements
# of the program
dds2 <- DESeqDataSetFromMatrix(countData = L02_counts,
                              colData = metadata,
                              design = ~ combined_factor)
dds2

# ~~ re-level data ~~~
# location data leveled based on depth
dds2$location <- factor(dds2$location, levels = c("new_tunnel", 
                                                "35_meters", "45_meters", "60_meters", "83_meters"))

# re-leveling thaw state over temp gradient (frozen to thawed)
dds2$thaw.temp <- factor(dds2$thaw.temp, levels = c("frozen", "thawing", "thawed"))

# ~~~ pre-filtering taken care of in BBTools ~~~


# DESeq ANALYSIS ----------------------------------------


  # automatic independent filtering = TRUE
  # alpha = 0.1
dds2 <- DESeq(dds2)  
res2 <- results(dds2)

# set range for color gradient
colorRange <- c(-6,5)

# contrast results --------------------------------------
  # creates specific comparison of frozen->thawing and frozen -> at each location
  # output is a DESEqResults table for the contrast

res35thawed_frozen2 <- results(dds2, contrast=c("combined_factor", "thawed35_meters", "frozen35_meters"))
res35thawing_frozen2 <- results(dds2, contrast=c("combined_factor", "thawing35_meters", "frozen35_meters"))

res45thawed_frozen2 <- results(dds2, contrast=c("combined_factor", "thawed45_meters", "frozen45_meters"))
res45thawing_frozen2 <- results(dds2, contrast=c("combined_factor", "thawing45_meters", "frozen45_meters"))

res60thawed_frozen2 <- results(dds2, contrast=c("combined_factor", "thawed60_meters", "frozen60_meters"))
res60thawing_frozen2 <- results(dds2, contrast=c("combined_factor", "thawing60_meters", "frozen60_meters"))

res83thawed_frozen2 <- results(dds2, contrast=c("combined_factor", "thawed83_meters", "frozen83_meters"))
res83thawing_frozen2 <- results(dds2, contrast=c("combined_factor", "thawing83_meters", "frozen83_meters"))

resNewThawed_frozen2 <- results(dds2, contrast=c("combined_factor", "thawednew_tunnel", "frozennew_tunnel"))
resNewThawing_frozen2 <- results(dds2, contrast=c("combined_factor", "thawingnew_tunnel", "frozennew_tunnel"))


# export to .csv ----------------------------------------
  # output is a csv file containing:
    # L03 subsystem name
    # Log2 Fold Change for frozen->thawing
    # Adjusted p-value for frozen->thawing
    # Log2 Fold Change for frozen->thawed
    # Adjusted p-value for frozen->thawed

res35_2 <- c(as.data.frame(res35thawed_frozen2@rownames), 
             as.data.frame(res35thawing_frozen2$log2FoldChange),
             as.data.frame(res35thawing_frozen2$padj),
             as.data.frame(res35thawed_frozen2$log2FoldChange), 
             as.data.frame(res35thawed_frozen2$padj))
write.csv(res35_2, file="res35_2.csv")

res45_2 <- c(as.data.frame(res45thawed_frozen2@rownames), 
             as.data.frame(res45thawing_frozen2$log2FoldChange),
             as.data.frame(res45thawing_frozen2$padj),
             as.data.frame(res45thawed_frozen2$log2FoldChange), 
             as.data.frame(res45thawed_frozen2$padj))
write.csv(res45_2, file="res45_2.csv")

res60_2 <- c(as.data.frame(res60thawed_frozen2@rownames), 
             as.data.frame(res60thawing_frozen2$log2FoldChange),
             as.data.frame(res60thawing_frozen2$padj),
             as.data.frame(res60thawed_frozen2$log2FoldChange), 
             as.data.frame(res60thawed_frozen2$padj))
write.csv(res60_2, file="res60_2.csv")

res83_2 <- c(as.data.frame(res83thawed_frozen2@rownames), 
             as.data.frame(res83thawing_frozen2$log2FoldChange),
             as.data.frame(res83thawing_frozen2$padj),
             as.data.frame(res83thawed_frozen2$log2FoldChange), 
             as.data.frame(res83thawed_frozen2$padj))
write.csv(res83_2, file="res83_2.csv")

resNew_2 <- c(as.data.frame(resNewThawed_frozen2@rownames), 
              as.data.frame(resNewThawing_frozen2$log2FoldChange),
              as.data.frame(resNewThawing_frozen2$padj),
              as.data.frame(resNewThawed_frozen2$log2FoldChange), 
              as.data.frame(resNewThawed_frozen2$padj))
write.csv(resNew_2, file="resNew_2.csv")




# 35 Meters ---------------------------------------------

# read in result csv files as a tibble
res35csv2 <- readr::read_csv("res35_2.csv", 
                            col_names = TRUE, 
                            col_types= "icdddd")
res35csv2 <- select(res35csv2, -1) # delete first column of numbers

# re-name columns
colnames(res35csv2) <- c("L02", "thawing", "thawing.padj", "thawed", "thawed.padj")

# gather and remove log2 NA rows
df2 <- res35csv2 %>%
  gather(key="thaw.state", value= "Log2", 
         thawing, thawed, -thawing.padj, -thawed.padj, 
         na.rm=TRUE)

# re-level thaw states
df2$thaw.state <- factor(df2$thaw.state,levels = c("thawing", "thawed"))

meters35_2 <- ggplot(df2, aes(x=thaw.state, y =L02 , fill=Log2)) +
  # tile with black contour
  geom_tile(colour="black") +
  # B&W theme, no grey background
  theme_bw() +
  # get rid of y axis titles
  theme(axis.ticks.y=element_blank()) +
  # set palette
  scale_fill_distiller(palette="RdBu", 
                       direction =-1,
                       na.value = "grey50",
                       limits=c(-7.1, 7.1)) +
  # since there is no legend, adding a title
  labs(title = "35 Meters: LFC from frozen")

meters35_2


# save to file in figures directory
ggsave(
  "35meters_L02.png",
  plot = meters35_2,
  device = png(),
  path = "~/Documents/CRREL/permafrost-pathogens/permafrost-figures",
  scale = 1,
  width = 50,
  height = 30,
  units = "cm",
  dpi = 300)

# 45 Meters ---------------------------------------------


# read in result csv files as a tibble
res45csv2 <- readr::read_csv("res45_2.csv", 
                             col_names = TRUE, 
                             col_types= "icdddd")
res45csv2 <- select(res45csv2, -1) # delete first column of numbers

# re-name columns
colnames(res45csv2) <- c("L02", "thawing", "thawing.padj", "thawed", "thawed.padj")

# gather and remove log2 NA rows
df2 <- res45csv2 %>%
  gather(key="thaw.state", value= "Log2", 
         thawing, thawed, -thawing.padj, -thawed.padj, 
         na.rm=TRUE)

# re-level thaw states
df2$thaw.state <- factor(df2$thaw.state,levels = c("thawing", "thawed"))

meters45_2 <- ggplot(df2, aes(x=thaw.state, y = L02 , fill=Log2)) +
  # tile with black contour
  geom_tile(colour="black") +
  # B&W theme, no grey background
  theme_bw() +
  # get rid of y axis titles
  theme(axis.ticks.y=element_blank()) +
  # Green color theme for `fill`
  scale_fill_distiller(palette="RdBu", 
                       direction=-1,
                       na.value = "grey50",
                       limits= c(-7.1, 7.1)) +
  # since there is no legend, adding a title
  labs(title = "45 Meters: LFC from frozen")
meters45_2

# save to file in figures directory
ggsave(
  "45meters_L02.png",
  plot = meters45_2,
  device = png(),
  path = "~/Documents/CRREL/permafrost-pathogens/permafrost-figures",
  scale = 1,
  width = 50,
  height = 30,
  units = "cm",
  dpi = 300)



# 60 Meters ---------------------------------------------

# read in result csv files as a tibble
res60csv2 <- readr::read_csv("res60_2.csv", 
                             col_names = TRUE, 
                             col_types= "icdddd")
res60csv2 <- select(res60csv2, -1) # delete first column of numbers

# re-name columns
colnames(res60csv2) <- c("L02", "thawing", "thawing.padj", "thawed", "thawed.padj")

# gather and remove log2 NA rows
df2 <- res60csv2 %>%
  gather(key="thaw.state", value= "Log2", 
         thawing, thawed, -thawing.padj, -thawed.padj, 
         na.rm=TRUE)

# re-level thaw states
df2$thaw.state <- factor(df2$thaw.state,levels = c("thawing", "thawed"))

meters60_2 <- ggplot(df2, aes(x=thaw.state, y = L02 , fill=Log2)) +
  # tile with black contour
  geom_tile(colour="black") +
  # B&W theme, no grey background
  theme_bw() +
  # get rid of y axis titles
  theme(axis.ticks.y=element_blank()) +
  # Green color theme for `fill`
  scale_fill_distiller(palette="RdBu", 
                       direction=-1,
                       na.value = "grey50",
                       limits= c(-7.1, 7.1)) +
  # since there is no legend, adding a title
  labs(title = "60 Meters: Log Fold Change from Frozen")
meters60_2

# save to file in figures directory
ggsave(
  "60meters_L02.png",
  plot = meters60_2,
  device = png(),
  path = "~/Documents/CRREL/permafrost-pathogens/permafrost-figures",
  scale = 1,
  width = 50,
  height = 30,
  units = "cm",
  dpi = 300)


# 83 Meters ---------------------------------------------

# read in result csv files as a tibble
res83csv2 <- readr::read_csv("res83_2.csv", 
                             col_names = TRUE, 
                             col_types= "icdddd")
res83csv2 <- select(res83csv2, -1) # delete first column of numbers

# re-name columns
colnames(res83csv2) <- c("L02", "thawing", "thawing.padj", "thawed", "thawed.padj")

# gather and remove log2 NA rows
df2 <- res83csv2 %>%
  gather(key="thaw.state", value= "Log2", 
         thawing, thawed, -thawing.padj, -thawed.padj, 
         na.rm=TRUE)

# re-level thaw states
df2$thaw.state <- factor(df2$thaw.state,levels = c("thawing", "thawed"))

meters83_2 <- ggplot(df2, aes(x=thaw.state, y = L02 , fill=Log2)) +
  # tile with black contour
  geom_tile(colour="black") +
  # B&W theme, no grey background
  theme_bw() +
  # get rid of y axis titles
  theme(axis.ticks.y=element_blank()) +
  # Green color theme for `fill`
  scale_fill_distiller(palette="RdBu", 
                       direction=-1,
                       na.value = "grey50",
                       limits= c(-7.1, 7.1)) +
  # since there is no legend, adding a title
  labs(title = "83 Meters: Log Fold Change from Frozen")
meters83_2
summary(res2$Log2)
# save to file in figures directory
ggsave(
  "83meters_L02.png",
  plot = meters83_2,
  device = png(),
  path = "~/Documents/CRREL/permafrost-pathogens/permafrost-figures",
  scale = 1,
  width = 50,
  height = 30,
  units = "cm",
  dpi = 300)

# New Tunnel ---------------------------------------------

# read in result csv files as a tibble
resNewcsv2 <- readr::read_csv("resNew_2.csv", 
                             col_names = TRUE, 
                             col_types= "icdddd")
resNewcsv2 <- select(resNewcsv2, -1) # delete first column of numbers

# re-name columns
colnames(resNewcsv2) <- c("L02", "thawing", "thawing.padj", "thawed", "thawed.padj")

# gather and remove log2 NA rows
df2 <- resNewcsv2 %>%
  gather(key="thaw.state", value= "Log2", 
         thawing, thawed, -thawing.padj, -thawed.padj, 
         na.rm=TRUE)

# re-level thaw states
df2$thaw.state <- factor(df2$thaw.state,levels = c("thawing", "thawed"))

metersNew_2 <- ggplot(df2, aes(x=thaw.state, y = L02 , fill=Log2)) +
  # tile with black contour
  geom_tile(colour="black") +
  # B&W theme, no grey background
  theme_bw() +
  # get rid of y axis titles
  theme(axis.ticks.y=element_blank()) +
  # Green color theme for `fill`
  scale_fill_distiller(palette="RdBu", 
                       direction=-1,
                       na.value = "grey50",
                       limits= c(-7.1, 7.1)) +
  # since there is no legend, adding a title
  labs(title = "New Tunnel: Log Fold Change from Frozen")
metersNew_2

# save to file in figures directory
ggsave(
  "NewTunnel_L02.png",
  plot = metersNew_2,
  device = png(),
  path = "~/Documents/CRREL/permafrost-pathogens/permafrost-figures",
  scale = 1,
  width = 50,
  height = 30,
  units = "cm",
  dpi = 300)





