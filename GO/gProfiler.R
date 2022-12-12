
setwd("/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/gene_ontology")

library(ggplot2)
library(cowplot)
library(limma) 
library(readxl)
library(gprofiler2)
library(tidyverse)


# Get the results table from the spatial statistics script. 
genes1 = read.table("/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/gene_ontology/databases/trackon/dwi/results_rankcor_002_all_knn10_geo.csv", header = TRUE, sep=",")

#Calculate the average
avg_expression<-aggregate(genes1[, -c(1:3)], by = list(genes1$Gene),mean, na.rm = TRUE)

# Filter values that are positively or negatively correlated
pos_avg_expression<-avg_expression %>% filter(T>0) 
neg_avg_expression<-avg_expression %>% filter(T<0) 

#Get the genes from AIBS list
genes2 = read.table("AIBS_gene_list.csv", header = FALSE, stringsAsFactors=FALSE, sep=",")
aibs_genes = genes2$V1

#Keep only the genes that are also in AIBS
myfiltered_pos_genes<-subset(pos_avg_expression, (pos_avg_expression$Group.1 %in% aibs_genes))
myfiltered_neg_genes<-subset(neg_avg_expression, (neg_avg_expression$Group.1 %in% aibs_genes))


# Obtain the genes that are positively correlated (upper 10, 20, 30%) AND SORT THEM BY T VALUE
q90 = quantile(myfiltered_pos_genes$T, 0.90)
unsorted_upper10<-filter(myfiltered_pos_genes, T>q90)
write.csv(unsorted_upper10, file = "/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/gene_ontology/databases/trackon/dwi/unsorted_filtered_list_upper10.csv")
upper10 <- unsorted_upper10[order(unsorted_upper10[, "T"], decreasing = TRUE, na.last = FALSE), , drop = FALSE]
write.csv(upper10, file = "/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/gene_ontology/databases/trackon/dwi/sorted_filtered_list_upper10.csv")


q80 = quantile(myfiltered_pos_genes$T, 0.80)
unsorted_upper20<-filter(myfiltered_pos_genes, T>q80)
write.csv(unsorted_upper20, file = "/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/gene_ontology/databases/trackon/dwi/unsorted_filtered_list_upper20.csv")
upper20 <- unsorted_upper20[order(unsorted_upper20[, "T"], decreasing = TRUE, na.last = FALSE), , drop = FALSE]
write.csv(upper20, file = "/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/gene_ontology/databases/trackon/dwi/sorted_filtered_list_upper20.csv")

q70 = quantile(myfiltered_pos_genes$T, 0.70)
unsorted_upper30<-filter(myfiltered_pos_genes, T>q70)
write.csv(unsorted_upper30, file = "/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/gene_ontology/databases/trackon/dwi/unsorted_filtered_list_upper30.csv")
upper30 <- unsorted_upper30[order(unsorted_upper30[, "T"], decreasing = TRUE, na.last = FALSE), , drop = FALSE]


# Obtain the genes that are more negatively correlated (lower 10, 20, 30%)
q10 = quantile(myfiltered_neg_genes$T, 0.10)
unsorted_lower10<-filter(myfiltered_neg_genes, T<q10)
write.csv(unsorted_lower10, file = "/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/gene_ontology/databases/trackon/dwi/unsorted_filtered_list_down10.csv")
lower10 <- unsorted_lower10[order(unsorted_lower10[, "T"], decreasing = FALSE, na.last = FALSE), , drop = FALSE]
write.csv(lower10, file = "/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/gene_ontology/databases/trackon/dwi/sorted_filtered_list_down10.csv")



q20 = quantile(myfiltered_neg_genes$T, 0.20)
unsorted_lower20<-filter(myfiltered_neg_genes, T<q20)
write.csv(unsorted_lower20, file = "/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/gene_ontology/databases/trackon/dwi/filtered_list_down20.csv")
lower20 <- unsorted_lower20[order(unsorted_lower20[, "T"], decreasing = FALSE, na.last = FALSE), , drop = FALSE]
write.csv(lower20, file = "/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/gene_ontology/databases/trackon/dwi/sorted_filtered_list_down20.csv")


q30 = quantile(myfiltered_neg_genes$T, 0.30)
unsorted_lower30<-filter(myfiltered_neg_genes, T<q30)
write.csv(unsorted_lower30, file = "/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/gene_ontology/databases/trackon/dwi/filtered_list_down30.csv")
lower30 <- unsorted_lower30[order(unsorted_lower30[, "T"], decreasing = FALSE, na.last = FALSE), , drop = FALSE]
write.csv(lower30, file = "/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/gene_ontology/databases/trackon/dwi/sorted_filtered_list_down30.csv")




#####This script assumes that the average expression between probes has been calculated (one gene=one T value)
#####Also assumes that there are *.csv files with 10%, 20%, 30% of genes up/down regulated RANKED BY T VALUE 
##### I did this already in the gp_profiler script and also saved in preprocessing analysis

################################################################## Up-regulated #######################################################
## TOP 10%

lhits = upper10$Group.1
gostres <- gost(lhits, organism = "hsapiens", correction_method = "g_SCS", sources = c("GO:BP", "GO:MF","GO:CC", "KEGG", "REAC", "CORUM", "HPA"))
head(gostres$result)
data<-data.frame(lapply(gostres$result, as.character), stringsAsFactors=FALSE)
df = as.matrix(data)
###write.table(data, file="Gprofiler_CSV/gp_10_UR.csv", row.names=FALSE, na="", col.names = TRUE, sep=",")
write.table(data, file="/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/gene_ontology/GO_results/trackon/dwi/gp_10_UR.csv", row.names=FALSE, na="", col.names = TRUE, sep=",")

##TO PLOT IT
graph_upper10 <- gostplot(gostres, capped = FALSE, interactive = FALSE)
tiff("/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/gene_ontology/GO_results/trackon/dwi/graph_upper10.tiff")
plot (graph_upper10)
dev.off()



## TOP 20%
lhits = upper20$Group.1
gostres <- gost(lhits, organism = "hsapiens", correction_method = "g_SCS", sources = c("GO:BP", "GO:MF","GO:CC", "KEGG", "REAC", "CORUM", "HPA"))
head(gostres$result)
data<-data.frame(lapply(gostres$result, as.character), stringsAsFactors=FALSE)
df = as.matrix(data)
write.table(data, file="/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/gene_ontology/GO_results/trackon/dwi/gp_20_UR.csv", row.names=FALSE, na="", col.names = TRUE, sep=",")

##TO PLOT IT
graph_upper20 <- gostplot(gostres, capped = FALSE, interactive = FALSE)
tiff("/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/gene_ontology/GO_results/trackon/dwi/graph_upper20.tiff")
plot (graph_upper20)
dev.off()



## TOP 30%
lhits = upper30$Group.1
gostres <- gost(lhits, organism = "hsapiens", correction_method = "g_SCS", sources = c("GO:BP", "GO:MF","GO:CC", "KEGG", "REAC", "CORUM", "HPA"))
head(gostres$result)
data<-data.frame(lapply(gostres$result, as.character), stringsAsFactors=FALSE)
df = as.matrix(data)
write.table(data, file="/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/gene_ontology/GO_results/trackon/dwi/gp_30_UR.csv", row.names=FALSE, na="", col.names = TRUE, sep=",")

##TO PLOT IT
graph_upper30 <- gostplot(gostres, capped = FALSE, interactive = FALSE)
tiff("/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/gene_ontology/GO_results/trackon/dwi/graph_upper30.tiff")
plot (graph_upper30)
dev.off()


################################################################## Down-regulated #######################################################



## TOP 10%
lhits = lower10$Group.1
gostres <- gost(lhits, organism = "hsapiens", correction_method = "g_SCS", sources = c("GO:BP", "GO:MF","GO:CC", "KEGG", "REAC", "CORUM", "HPA"))
head(gostres$result)
data<-data.frame(lapply(gostres$result, as.character), stringsAsFactors=FALSE)
df = as.matrix(data)
###write.table(data, file="Gprofiler_CSV/gp_10_UR.csv", row.names=FALSE, na="", col.names = TRUE, sep=",")
write.table(data, file="/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/gene_ontology/GO_results/trackon/dwi/gp_10_DR.csv", row.names=FALSE, na="", col.names = TRUE, sep=",")

##TO PLOT IT
graph_lower10 <- gostplot(gostres, capped = FALSE, interactive = FALSE)
tiff("/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/gene_ontology/GO_results/trackon/dwi/graph_lower10.tiff")
plot (graph_lower10)
dev.off()

## TOP 20%
lhits = lower20$Group.1
gostres <- gost(lhits, organism = "hsapiens", correction_method = "g_SCS", sources = c("GO:BP", "GO:MF","GO:CC", "KEGG", "REAC", "CORUM", "HPA"))
head(gostres$result)
data<-data.frame(lapply(gostres$result, as.character), stringsAsFactors=FALSE)
df = as.matrix(data)
write.table(data, file="/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/gene_ontology/GO_results/trackon/dwi/gp_20_DR.csv", row.names=FALSE, na="", col.names = TRUE, sep=",")

##TO PLOT IT
graph_lower20 <- gostplot(gostres, capped = FALSE, interactive = FALSE)
tiff("/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/gene_ontology/GO_results/trackon/dwi/graph_lower20.tiff")
plot (graph_lower20)
dev.off()


## TOP 30%
lhits = lower30$Group.1
gostres <- gost(lhits, organism = "hsapiens", correction_method = "g_SCS", sources = c("GO:BP", "GO:MF","GO:CC", "KEGG", "REAC", "CORUM", "HPA"))
head(gostres$result)
data<-data.frame(lapply(gostres$result, as.character), stringsAsFactors=FALSE)
df = as.matrix(data)
write.table(data, file="/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/gene_ontology/GO_results/trackon/dwi/gp_30_DR.csv", row.names=FALSE, na="", col.names = TRUE, sep=",")

##TO PLOT IT
graph_lower30 <- gostplot(gostres, capped = FALSE, interactive = FALSE)
tiff("/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/gene_ontology/GO_results/trackon/dwi/graph_lower30.tiff")
plot (graph_lower30)
dev.off()


