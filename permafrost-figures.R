# PERMAFROST-PATHOGENS DATA VISUALIZATION
  # creates heat maps looking at frozen to thawed and frozen to
  # thawing LFC for all L03 virulence-related genes at each site
# Created: 7/20/2020
# Author: C. Green
# Last Edited: 7/24/2020

# ~~ heat map ~~~

# load library
library("dplyr")
library("tidyr")
library("ggplot2")

# ~~ 35 METERS ~~~

# read in result csv files as a tibble
res35csv <- readr::read_csv("res35.csv", 
                            col_names = TRUE, 
                            col_types= list("icdddd"))
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


# ~~ New METERS ~~~

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

