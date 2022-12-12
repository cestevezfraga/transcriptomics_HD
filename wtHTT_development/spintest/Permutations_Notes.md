# Run spin permutations
# Developed by Angelika Zarkali https://github.com/AngelikaZa

## In Rstudio 

### Set working directory
setwd("C:/Users/Angelika/Dropbox/UCL/02_EXPERIMENTS/OTHER/PreclinicalHD_Nov2022/") 

###  Load necessary functions
library(matrixStats)
source("rotate.parcellation.R")
source("perm.sphere.p.R")

### Load coordinates as arrays
###### left coordinates
coord_l = scan("coordinates_DesikanCerebellum_L.txt")
coord_l = matrix(coord_l, ncol = 3, byrow = TRUE)
###### right coordinates
coord_r = scan("coordinates_DesikanCerebellum_R.txt")
coord_r = matrix(coord_r, ncol = 3, byrow = TRUE)

### Run 1000 spin permutations
perm_matrix = rotate.parcellation(coord.l = coord_l,coord.r = coord_r,nrot = 1000)

### Export to txt
write.table(perm_matrix, file="permutations_Glasser.txt", row.names=FALSE, col.names=FALSE)
