# intersect code between two lists of species

# install reading packages
install.packages("readxl")
library("readxl")

# install.packages("dplyr")


# upload xlxs file
sampleA <- read_excel("sampleA_table.xlsx")
sampleB <- read_excel("sampleB_table.xlsx")

# test for individual species first
sp1 <- "Helicobacter pylori"
sp2 <- "Campylobacter jejuni"
sp3 <- "Helicobacter coli"
grepl(sp1, sp2, fixed = TRUE)
grepl(sp1, sp3, fixed = TRUE)
grepl(sp3, sp1, fixed = TRUE)


                       