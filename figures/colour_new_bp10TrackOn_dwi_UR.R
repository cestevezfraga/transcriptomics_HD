# install.packages
# install.packages("radiant")
# install.packages("Rtools")
install.packages("wesanderson")

#library
library(ggplot2)
library (readxl)
library(tidyverse)
library(stringr)
library(RColorBrewer)
library(qdapTools)
library (dplyr)
library(forcats)
library(wesanderson)

# library(radiant)
# library(Rtools)

# load dataset

data <- read.table("/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/sorting_out_resutls/graphs/bp/trackon_dwi_ur.csv", header=TRUE, sep=",",)


# Data wrangling

#Functional_block_summary <- data 
 # dplyr:: filter (! Functional_Block %in% c("general","metabolism","physiology", "intracellular organisation")) %>% 
  #unique() %>%
  #group_by(Functional_Block) %>%
  #summarise(number = n()) %>% 

sort(data$number)

# Stacked plot
# geom_bar(position="fill") for percent on y-axis

ggplot <- 
  ggplot(data = data, aes(y = data$number, x = reorder(data$function.,data$number)), colour = "black") + 
  geom_col(aes(fill = data$number, y = data$number), colour = "black") +
  coord_flip()+   #change orientation of the plot
  scale_x_discrete(name = " ", ) +
  scale_y_continuous(name = "Number of enriched GO terms") +
  theme_classic() +   # with no grids in the background
  theme(axis.text.x = element_text(size = 25, angle = 0, hjust = 0.5), legend.position = "none")+
  theme(axis.text.y = element_text(size = 25, angle = 0), legend.position = "none")+
  ggtitle("GO terms positively associated with cortical mean diffusivity") +
  theme(plot.title = element_text(face = "bold", size = 20, hjust=1.4))+
  scale_fill_gradient(low = "red", high = "red") # only for continous values

print(ggplot)

# save plot

ggsave(filename = "COLOUR_NEW_BP_TrackOn_dwi_UR.png", 
       plot = ggplot, width = 10, 
       path = "/Users/charlie/Desktop/phd_fromaug/cs2_ANALYSIS/VBM/AIBS_FTD/sorting_out_resutls/graphs/", 
       height = 10, 
       dpi = 300)