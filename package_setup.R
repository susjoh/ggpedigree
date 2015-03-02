
library(kinship2)
library(ggplot2)
library(grid)
library(reshape)
library(plyr)



pedigree <- read.table("data/SoayPedigree.txt", header = T)

# pedigree <- rbind(cbind(Family = 1, read.table("data/test_78.ped", header = T)),
#                   cbind(Family = 2, read.table("data/test_38.ped", header = T)),
#                   cbind(Family = 3, read.table("data/test_31.ped", header = T)))

pedigree <- read.table("data/test_78.ped", header = T)
pedigree <- read.table("../../SNP_pedigree_2014-10-23.txt", header = T)

ggpedigree(pedigree)