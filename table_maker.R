library(readr)
library(stringr)
library(data.table)
metadata <- read_delim("~/Desktop/TCR/data/metadata.csv", 
                       " ", escape_double = FALSE, trim_ws = TRUE)

metadata_SMART <- read_delim("~/Desktop/TCR/data/metadata_SMART.csv", " ", escape_double = FALSE, trim_ws = TRUE)

metadata_SMART$SMART <- gsub("MK-", "", metadata_SMART$name)
metadata_SMART

merged <- merge(metadata, metadata_SMART, "SMART")
merged
out <- data.frame(SMART=merged$SMART, barcode=merged$`barcode (demultiplex)`, Empty=NaN, Read_1=merged$Read_1, Read_2=merged$Read_2)
out
fwrite(out, file = "barcodes.txt", sep = "\t")


out2 <- data.frame(MK=merged$MK, Seq=merged$Seq)
out2
fwrite(out2, file = "out.txt", sep = "\t")
