#SCATTERPLOT DWI - GMV

library(ggplot2)
library (readxl)
library(tidyverse)
library(stringr)
library(RColorBrewer)
library(qdapTools)
library (dplyr)
library(forcats)
library(wesanderson)

data <- read.table("/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/corr_dwi-gmv/dwi_gm_nozero.csv", header=TRUE, sep=",",)

p<-ggplot(data=data, aes(x=data$gm_original, y=data$dwi_original)) + geom_point(size=2, shape=18) + geom_smooth(method=lm) + 
  labs (x="Grey matter volume", y = "Mean Diffusivity") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave(p, filename = 'correlation_gmv_md_nozero.pdf', device = 'pdf', width = 11.69, height = 8.27, dpi = 300, path ="/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/corr_dwi-gmv/graphs/" )


cor.test (data$gm_original, data$dwi_original )