library(SpatialPack)
library(spdep)
library(ade4)
library(deldir)
library(data.table)
library(feather)
library(plot3D)
library(arrow)
spm_map="002"
#check if adespatial is available, if not load a few helper functions
has_adespatial = require(adespatial)
if (has_adespatial){
  require(adegraphics)
} else {
  #contains the key functions for adespatial for this analysis
  source("./helper_functions.R")
}


### define some functions ###
#compute test for correlation between rank(imgv) and rank(gene expression)
#also computes the p-values for Moran's I (test for spatial autocorrelation)
adj.test <- function(g, mydf, form_txt, max_mem){
  #using rank gives us spearman rank correlation
  formula.txt <- form_txt
  if (max_mem > 0){
    formula.txt <- paste(form_txt, " + ", paste("MEM",1:max_mem, collapse=" + ", sep=""))
  }
  m2 <- lm(as.formula(formula.txt), data=mydf)
  r2 <- m2$residuals
  moran.p <- moran.test(r2, mylistw)$p.value
  cor.t <- summary(m2)$coeff["rank(imgv)",3:4]
  return(c(cor.t, moran.p))
}

#runs the above function with different numbers of spatial eigenvectors (MEM)
adj.test.screen <- function(g, mydf, form_text="", max_mem=100){
  mems <- c(1,seq(5,max_mem,5))
  for (m in mems){
    res <- adj.test(g, mydf, form_text, m)
    if (res[3] > 0.05){
      return( c(res, m))
    }
  }
  return(c(res, m))
}

plot3dBrain <- function(coords, adj, pcol="black", lcol="red", my_theta=45, my_phi=45){

  #points3D(x=coords[,1], y=coords[,2], z=coords[,3], pch=19, col=pcol, cex=0.5, theta=my_theta, phi=my_phi, lwd=2)
  xxx <- c()
  for(i in 1:nrow(coords)){
    p0 <- unlist(coords[i,])
    #print(p0)
    js <- which(adj[i,]>0)
    for(j in js){
      p1 <- unlist(coords[j,])
      xxx <- rbind(xxx, c(p0,p1))
    }
  }
  segments3D(x0=xxx[,1], y0=xxx[,2], z0=xxx[,3], x1=xxx[,4], y1=xxx[,5], z1=xxx[,6], col=lcol, theta=my_theta, phi=my_phi, lwd=1, add=F)
  #points3D(x=coords[,1], y=coords[,2], z=coords[,3], pch=19, col=pcol, cex=0.5, theta=my_theta, phi=my_phi, lwd=2)

}

### options ###

#these vars can be set remotely
if (!exists("remotely")){
  #selected samples
  spm_map <- "002"

  #pick on of these:
  samples <- c("H0351_1009", "H0351_1012", "H0351_1015", "H0351_1016", "H0351_2001", "H0351_2002")
  #sid = "H0351_2001"
  sid = ""
}

sample_file  <- "/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/data/selected_samples_lh_for_FTD_AIBS.csv"
probe_file   <- "/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/data/probes_for_analysis.csv"
img_val_file <- paste("/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/data/results/trackon/dwi/image_values_", spm_map, ".csv", sep="")



knn_set <- "10"
max_mem <- 150

tmp <- ""
sid_str <- "all"
if (sid != ""){
  tmp <- "_subj"
  sid_str <- sid
}

#adjecency matrix to be used
#sample_adj_file <- "/Users/charlie/data/sample_adjecency.csv"
#sample_adj_file <- paste("/Users/charlie/Desktop/data/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/data/sample_adjecency_knn",knn_set, "_geo", tmp, ".csv",sep="")
sample_adj_file <- paste("/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/data/sample_adjecency_knn",knn_set, "_geo", tmp, ".csv",sep="")

#location of the feather file generated by the the first python script
GE_file <- arrow::read_feather("/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/data/AIBS/GE.feather")
GE <- as.data.frame(GE_file)
rownames(GE) <- GE[,"probe_id"]

ofname <- paste("/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/data/results/trackon/dwi/results_rankcor_", spm_map, "_", sid_str, "_knn", knn_set, "_geo",tmp, ".csv", sep="")


### main ###
#simplest model: rank-correlation between gene expression and imaging value
model_str <- "rank(g) ~ rank(imgv)"

#load sample information
sample_info <- read.csv(sample_file, as.is=T, row.names=1)
rownames(sample_info) <- sample_info$well_id

#graph based on sample_adjecency
#sa <- fread(sample_adj_file, data.table=F)
sa <- fread('/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/data/sample_adjecency_knn10_geo.csv', data.table=F, header=TRUE)
rownames(sa) <- sa$well_id
sa <- sa[,2:ncol(sa)]
diag(sa) <- 0

### load image values ###
img_vals <- fread(img_val_file, data.table=F)
rownames(img_vals) <- img_vals$well_id

### work with genetics data ###
probe_info <- read.csv(probe_file, as.is=T)
rownames(probe_info) <- probe_info$probe_id

## reduce to samples from a given subject:
sample_info2 <- sample_info
if (sid != ""){
  sample_info2 <- subset(sample_info, subset=sampleID==sid)
} else {
  model_str <- paste(model_str, "+ factor(sampleID)")
}

wids <- paste(sample_info2$well_id)
##build graph
#extract the coordinates
coords <- sample_info2[,c("corrected_mni_x","corrected_mni_y","corrected_mni_z")]

mylistw <- mat2listw(as.matrix(sa[wids,wids]))
#s.label(coords, nb=mylistw, pnb.edge.col="red", labels=c())

#compute MEM
aibs.mem <- mem(mylistw)

## select genetics data
GEuse <- GE[paste(probe_info$probe_id), paste(sample_info2$well_id)]

## select imaging data
imgv <- img_vals[paste(sample_info2$well_id),"original"]

my_df <- data.frame(imgv, sampleID=sample_info2$sampleID, aibs.mem)


#message('WARNING: dummy run: only the first 10 probes processed')
#message('uncomment next block to run the script for all probes')
#result <- t(apply(head(GEuse,10), 1, function(gg){
#  adj.test.screen(gg, my_df, model_str, max_mem)
#}))

####### USE THIS FOR THE REAL ANALYSIS #######
#that's the actual run!
result <- t(apply(GEuse, 1, function(gg){
  adj.test.screen(gg, my_df, model_str, max_mem)
}))

pinfo <- probe_info[rownames(result),c("probe_name","gene_symbol")]
result <- data.frame(pinfo, result)
colnames(result) <- c("Probe","Gene","T","P","moranP","nMEM")

message("writing results to: ", ofname)
write.csv(result, ofname)
