library(readr)
library("dplyr")

args <- commandArgs(trailingOnly = T)
print(args)

directory <- toString(args[1])
name <- toString(args[2])


`%+%` <- function(a, b) paste0(a, b)


df <- read_delim(directory %+% "/" %+% name, "\t")

#View(df)
df <- select(df, -c("V.genes","J.genes","J.allele")) #-c("J.genes","J.allele"))

colnames(df) <- c("count",	"frequency",	"CDR3nt",	"CDR3aa",	"V", "D",	"J",	"Vend",	"Dstart",	"Dend",	"Jstart")
#View(df)

write.table(df, directory %+% "/" %+% name, quote = F, sep = "\t", row.names = FALSE)

