################################################################################
# MAIN PROJECT
################################################################################

# Libraries

getwd()
# setwd("Bureau/Sapienza/BioInf/BioInf-Project") # <- Modify this
getwd() # <- To check you're in the right directory

################################################################################
# 1. Get Files
################################################################################

# /!\ REQUIREMENT : Pre-process the data in `utils.R`

# 1.1 Read BIOGRID Human Interactome
biogrid <- read.delim("Files/Biogrid.txt")
View(biogrid)

# 1.2 Read Huri
huri <- read.delim("Files/Huri.txt")
View(huri)

# 1.3 Read String
string <- read.delim("Files/String.txt")
View(string)

# 1.4 Read Reactome
reactome <- read.delim("Files/Reactome.txt")
View(reactome)

# 1.4 Cardiomyopathy
cardiomyopathy <- read.delim("Files/Cardiomyopathy.txt")
View(cardiomyopathy)

