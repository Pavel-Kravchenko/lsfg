library(readr)
library(dplyr)
library(ggplot2)
library(forcats)


args <- commandArgs(trailingOnly = T)
print(args)

directory <- toString(args[1])
name1 <- toString(args[2])
name2 <- toString(args[3])
save_name <- toString(args[4])


#print(directory, name1, name2, save_name)

`%+%` <- function(a, b) paste0(a, b)


source("/raid/users/pk/kmer_dig/kmer_full_withV.R")


target = directory %+% "/mixcr_results/export/" %+% name1
db = directory %+% "/mixcr_results/export/" %+% name2

Patient<-count_kmers(target = target, db = db, with.gaps = F,  target_preCalculated=T, db_precalculated=T)
data.sf<-fread(target)
data.pb<-fread(db)

png(file= directory %+% "/mixcr_results/export/" %+% save_name  %+% ".png", width=1500, height=900, res=120)
ggplot(Patient,aes(x=log2fc,y=minus_log10p,color=count.target,size=count.target))+geom_point()
dev.off()


