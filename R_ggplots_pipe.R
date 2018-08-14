library(readr)
library(dplyr)
library(ggplot2)
library(forcats)

args <- commandArgs(trailingOnly = T)

directory <- toString(args[4])
name <- toString(args[5])

"%+%" <- function(...){
  paste0(...)
}


X_clones <- read_delim(directory %+% "/" %+% name %+% "_clones.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)

table_clones <- tbl_df(X_clones)

png(file="plot1.png", width=1500, height=900, res=120)

table_clones %>%
  group_by(bestVGene) %>% 
  summarise(cloneFraction=sum(cloneFraction)) %>%
  ggplot(aes(fct_reorder(bestVGene,cloneFraction,.desc = T) ,cloneFraction))+
  geom_bar(stat="identity",fill="orange") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

