##VOLCANO PLOT##
##From https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html
##More info https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/rna-seq-viz-with-volcanoplot/tutorial.html#:~:text=A%20volcano%20plot%20is%20a,the%20most%20biologically%20significant%20genes.


#THIS ONE IS TO DO WITH THE SPREADSHEETS WHERE DIVERGENT T VALUES HAVE BEEN REMOVED
#Spreadsheet obtained with the script in /Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/volcano_plots/no_divergent_T/

#You need a .csv file with the columns "P" (p value) "T" (T value) and dwi(gene name) 

library(ggrepel)
library(ggplot2)
library(data.table)



dwi <- read.table("/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/volcano_plots/no_divergent_T/dwi/avg_dwi.csv", header = TRUE, sep=",")

dwi$diffexpressed <- 'NO'
dwi$diffexpressed[dwi$T > 2 & dwi$P < 0.05 ] <- 'UP'
dwi$diffexpressed[dwi$T < -2 & dwi$P < 0.05 ] <- 'DOWN'
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
dwi$delabel <- NA
dwi$delabel[dwi$diffexpressed != 'NO'] <-dwi$dwi[dwi$diffexpressed != 'NO']


p<-ggplot(data=dwi, (aes(x=T, y=-log10(P), col=diffexpressed, label = delabel))) + geom_point() + theme_classic()  + geom_hline(yintercept=-log10(0.05), col="grey") + geom_vline(xintercept=c(-2, 2), col="grey") +
 scale_colour_manual(values = mycolors) + geom_text_repel() + ggtitle ('Cortical mean diffusivity') + theme(plot.title= element_text(hjust = 0.54, vjust=0.4, size = 20, face="bold")) + xlab("T value") + theme(legend.title = element_blank())

ggsave(p, filename = 'dwi_volcano_plot.pdf', device = 'pdf', width = 11.69, height = 8.27, dpi = 300, path ="/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/volcano_plots/no_divergent_T/dwi/graph/" )



