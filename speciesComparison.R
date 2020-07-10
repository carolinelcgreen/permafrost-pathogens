# intersect code between two lists of species

# install reading packages
install.packages("readxl")
library("readxl")

# upload xlxs file
sampleA <- read_excel("sampleA_table.xlsx")
sampleB <- read_excel("sampleB_table.xlsx")



for(i in sampleA){
  ifelse(sample[i,1] == "Brucella melitensis", "true", "false")
}

print(sampleA[3,1])

i=0
while(i <= length(sampleA)){
  print(sampleA[i,1])
  i+1;
}




                       