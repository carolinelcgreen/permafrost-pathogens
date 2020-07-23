# PERMAFROST-PATHOGENS DATA VISUALIZATION
# Created: 7/20/2020
# Author: C. Green
# Last Edited: 7/22/2020


# read in location results
# res35 <- read.csv("~/Documents/CRREL/permafrost-pathogens/res35.csv", row.names=1)
# res45 <- read.csv("~/Documents/CRREL/permafrost-pathogens/res45.csv", row.names=1)
# res60 <- read.csv("~/Documents/CRREL/permafrost-pathogens/res60.csv", row.names=1)
# res83 <- read.csv("~/Documents/CRREL/permafrost-pathogens/res83.csv", row.names=1)
# resNew <- read.csv("~/Documents/CRREL/permafrost-pathogens/resNew.csv", row.names=1)

# ~~ heat map ~~~

# load library
library("dplyr")
library("tidyr")
library("ggplot2")

# read in result csv files as a tibble
res35csv <- readr::read_csv("res35.csv", col_names = TRUE, col_types="cdddd")

# re name columns
colnames(res35csv) <- c("L03", "thawing", "thawing.padj", "thawed", "thawed.padj")

# gather and remove log2 NA rows
df <- res35csv %>%
  gather(key="thaw.state", value= "Log2", 
         thawing, thawed, -thawing.padj, -thawed.padj, 
         na.rm=TRUE)

ggplot(df, aes(x=L03, y = thaw.state , fill=Log2)) +
    # tile with black contour
    geom_tile(colour="black") 
    # B&W theme, no grey background
    # theme_bw() +
    # get rid of y axis titles
    # theme(axis.title.x=element_blank(),
          # axis.ticks.x=element_blank()) +
    # Green color theme for `fill`
    # scale_fill_distiller(palette="RdYlBu", direction=-1) + 
    # since there is no legend, adding a title
    # labs(title = "35 Meters: Log Fold Change from Frozen")






