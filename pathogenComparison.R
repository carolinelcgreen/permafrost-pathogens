

# read relevant Pathogens from Relevant Pathogens A-Z into a dataframe
# pathogens <-cbind(readLines('relevantPathogens.txt'))

# install reading packages
install.packages("readxl")
library("readxl")

install.packages("dplyr")


# upload xlxs file
Full_table <- read_excel("Full_table.xlsx")
Pathogen_table <- read_excel("pathogens.xlsx")

# find intersection
relevantPathogens <- c(intersect(Full_table[,7], Pathogen_table[,3])

# make a list
write.csv(relevantPathogens,"relevantPathogens_in_fullTable.csv", row.names = FALSE)






