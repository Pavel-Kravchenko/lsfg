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

merged_1 <- merge(out_103, out_113, by=c("aaSeqCDR3", "bestVGene"))
n_1 <- length(out_103[[1]]) + length(out_113[[1]])
print(n_1)
x_1 <- length(merged_1[[1]])
print(x_1)
bt_1 = binom.test(x = x_1, n = n_1, p = 0.5, alternative = 'less')
bt_1$p.value

merged_2 <- merge(out_135, out_137, by=c("aaSeqCDR3", "bestVGene"))
n_2 <- length(out_135[[1]]) + length(out_137[[1]])
print(n_2)
x_2 <- length(merged_2[[1]])
print(x_2)
bt_2 = binom.test(x = x_2, n = n_2, p = 0.5, alternative = 'less')
bt_2$p.value

merged_3 <- merge(out_143, out_165, by=c("aaSeqCDR3", "bestVGene"))
n_3 <- length(out_143[[1]]) + length(out_165[[1]])
print(n_3)
x_3 <- length(merged_3[[1]])
print(x_3)
bt_3 = binom.test(x = x_3, n = n_3, p = 0.5, alternative = 'less')
bt_3$p.value

merged_4 <- merge(out_204, out_205, by=c("aaSeqCDR3", "bestVGene"))
n_4 <- length(out_204[[1]]) + length(out_205[[1]])
print(n_4)
x_4 <- length(merged_4[[1]])
print(x_4)
bt_4 = binom.test(x = x_4, n = n_4, p = 0.5, alternative = 'less')
bt_4$p.value



# ntCDR3+Vgene
out_103 <- data.frame(nSeqCDR3=table_clones_103$nSeqCDR3, bestVGene=table_clones_103$bestVGene)
out_113 <- data.frame(nSeqCDR3=table_clones_113$nSeqCDR3, bestVGene=table_clones_113$bestVGene)
out_135 <- data.frame(nSeqCDR3=table_clones_135$nSeqCDR3, bestVGene=table_clones_135$bestVGene)
out_137 <- data.frame(nSeqCDR3=table_clones_137$nSeqCDR3, bestVGene=table_clones_137$bestVGene)
out_143 <- data.frame(nSeqCDR3=table_clones_143$nSeqCDR3, bestVGene=table_clones_143$bestVGene)
out_165 <- data.frame(nSeqCDR3=table_clones_165$nSeqCDR3, bestVGene=table_clones_165$bestVGene)
out_204 <- data.frame(nSeqCDR3=table_clones_204$nSeqCDR3, bestVGene=table_clones_204$bestVGene)
out_205 <- data.frame(nSeqCDR3=table_clones_205$nSeqCDR3, bestVGene=table_clones_205$bestVGene)

merged_1 <- merge(out_103, out_113, by=c("nSeqCDR3", "bestVGene"))
n_1 <- length(out_103[[1]]) + length(out_113[[1]])
#print(n_1)
x_1 <- length(merged_1[[1]])
#print(x_1)
bt_1 = binom.test(x = x_1, n = n_1, p = 0.5, alternative = 'less')
bt_1$p.value

merged_2 <- merge(out_135, out_137, by=c("nSeqCDR3", "bestVGene"))
n_2 <- length(out_135[[1]]) + length(out_137[[1]])
#print(n_2)
x_2 <- length(merged_2[[1]])
#print(x_2)
bt_2 = binom.test(x = x_2, n = n_2, p = 0.5, alternative = 'less')
bt_2$p.value

merged_3 <- merge(out_143, out_165, by=c("nSeqCDR3", "bestVGene"))
n_3 <- length(out_143[[1]]) + length(out_165[[1]])
#print(n_3)
x_3 <- length(merged_3[[1]])
#print(x_3)
bt_3 = binom.test(x = x_3, n = n_3, p = 0.5, alternative = 'less')
bt_3$p.value

merged_4 <- merge(out_204, out_205, by=c("nSeqCDR3", "bestVGene"))
n_4 <- length(out_204[[1]]) + length(out_205[[1]])
#print(n_4)
x_4 <- length(merged_4[[1]])
#print(x_4)
bt_4 = binom.test(x = x_4, n = n_4, p = 0.5, alternative = 'less')
bt_4$p.value

