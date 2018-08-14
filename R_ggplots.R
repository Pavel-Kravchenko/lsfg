library(readr)
library(ggplot2)
library(dplyr)

table_clones_103 <- read_delim("~/Desktop/TCR/data/103_clones.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)


table_clones_103 %>%
  group_by(bestVGene) %>% 
  summarise(cloneFraction=sum(cloneFraction)) %>%
  ggplot(aes(fct_reorder(bestVGene,cloneFraction,.desc = T) ,cloneFraction))+
  geom_bar(stat="identity",fill="orange") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

table_clones_103 %>%
  group_by(bestVGene) %>% 
  summarise(cloneCount=sum(cloneCount)) %>%
  ggplot(aes(fct_reorder(bestVGene,cloneCount,.desc = T) ,cloneCount))+
  geom_bar(stat="identity",fill="darkblue") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

table_clones_113 <- read_delim("~/Desktop/TCR/data/113_clones.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)

table_clones_135 <- read_delim("~/Desktop/TCR/data/135_clones.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)

table_clones_137 <- read_delim("~/Desktop/TCR/data/137_clones.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)

table_clones_143 <- read_delim("~/Desktop/TCR/data/143_clones.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)

table_clones_165 <- read_delim("~/Desktop/TCR/data/165_clones.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)

table_clones_204 <- read_delim("~/Desktop/TCR/data/204_clones.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)

table_clones_205 <- read_delim("~/Desktop/TCR/data/205_clones.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)


# aaCDR3+Vgene
out_103 <- data.frame(aaSeqCDR3=table_clones_103$aaSeqCDR3, bestVGene=table_clones_103$bestVGene)
out_113 <- data.frame(aaSeqCDR3=table_clones_113$aaSeqCDR3, bestVGene=table_clones_113$bestVGene)
out_135 <- data.frame(aaSeqCDR3=table_clones_135$aaSeqCDR3, bestVGene=table_clones_135$bestVGene)
out_137 <- data.frame(aaSeqCDR3=table_clones_137$aaSeqCDR3, bestVGene=table_clones_137$bestVGene)
out_143 <- data.frame(aaSeqCDR3=table_clones_143$aaSeqCDR3, bestVGene=table_clones_143$bestVGene)
out_165 <- data.frame(aaSeqCDR3=table_clones_165$aaSeqCDR3, bestVGene=table_clones_165$bestVGene)
out_204 <- data.frame(aaSeqCDR3=table_clones_204$aaSeqCDR3, bestVGene=table_clones_204$bestVGene)
out_205 <- data.frame(aaSeqCDR3=table_clones_205$aaSeqCDR3, bestVGene=table_clones_205$bestVGene)

print(length(out_113[[1]]) + length(out_137[[1]]) + length(out_143[[1]]) + length(out_204[[1]]))
merged_SF <- merge(out_113, out_137, by=c("aaSeqCDR3", "bestVGene"))
merged_SF <- merge(merged_SF, out_143, by=c("aaSeqCDR3", "bestVGene"))
merged_SF <- merge(merged_SF, out_204, by=c("aaSeqCDR3", "bestVGene"))
merged_SF
print(length(merged_SF[[1]]))

print(length(out_103[[1]]) + length(out_135[[1]]) + length(out_165[[1]]) + length(out_205[[1]]))
merged_PB <- merge(out_103, out_135, by=c("aaSeqCDR3", "bestVGene"))
merged_PB <- merge(merged_PB, out_165, by=c("aaSeqCDR3", "bestVGene"))
merged_PB <- merge(merged_PB, out_205, by=c("aaSeqCDR3", "bestVGene"))
merged_PB
print(length(merged_PB[[1]]))

merged <- merge(merged_PB, merged_SF, by=c("aaSeqCDR3", "bestVGene"))
print(length(merged[[1]]))

# ntCDR3+Vgene
out_103_2 <- data.frame(nSeqCDR3=table_clones_103$nSeqCDR3, bestVGene=table_clones_103$bestVGene)
out_113_2 <- data.frame(nSeqCDR3=table_clones_113$nSeqCDR3, bestVGene=table_clones_113$bestVGene)
out_135_2 <- data.frame(nSeqCDR3=table_clones_135$nSeqCDR3, bestVGene=table_clones_135$bestVGene)
out_137_2 <- data.frame(nSeqCDR3=table_clones_137$nSeqCDR3, bestVGene=table_clones_137$bestVGene)
out_143_2 <- data.frame(nSeqCDR3=table_clones_143$nSeqCDR3, bestVGene=table_clones_143$bestVGene)
out_165_2 <- data.frame(nSeqCDR3=table_clones_165$nSeqCDR3, bestVGene=table_clones_165$bestVGene)
out_204_2 <- data.frame(nSeqCDR3=table_clones_204$nSeqCDR3, bestVGene=table_clones_204$bestVGene)
out_205_2 <- data.frame(nSeqCDR3=table_clones_205$nSeqCDR3, bestVGene=table_clones_205$bestVGene)

print(length(out_113_2[[1]]) + length(out_137_2[[1]]) + length(out_143_2[[1]]) + length(out_204_2[[1]]))
merged_SF_2 <- merge(out_113_2, out_137_2, by=c("nSeqCDR3", "bestVGene"))
merged_SF_2 <- merge(merged_SF_2, out_143_2, by=c("nSeqCDR3", "bestVGene"))
merged_SF_2 <- merge(merged_SF_2, out_204_2, by=c("nSeqCDR3", "bestVGene"))
merged_SF_2
print(length(merged_SF_2[[1]]))

print(length(out_103_2[[1]]) + length(out_135_2[[1]]) + length(out_165_2[[1]]) + length(out_205_2[[1]]))
merged_PB_2 <- merge(out_103_2, out_135_2, by=c("nSeqCDR3", "bestVGene"))
merged_PB_2 <- merge(merged_PB_2, out_165_2, by=c("nSeqCDR3", "bestVGene"))
merged_PB_2 <- merge(merged_PB_2, out_205_2, by=c("nSeqCDR3", "bestVGene"))
merged_PB_2
print(length(merged_PB_2[[1]]))

merged_2 <- merge(merged_PB_2, merged_SF_2, by=c("nSeqCDR3", "bestVGene"))
print(length(merged_2[[1]]))
