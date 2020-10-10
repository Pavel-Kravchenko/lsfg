library(readr)
library(dplyr)
library(data.table)

args <- commandArgs(trailingOnly = T)


if (toString(args[1]) == "-d"){
  CD8_dir <- toString(args[2]) # first dir
  print("CD8_dir")
  print(CD8_dir)
  
  CD4_dir <- toString(args[3]) # second dir
  print("CD4_dir")
  print(CD4_dir)
  
  CD8_dir_files <- dir(CD8_dir)
  print("CD8_dir_files")
  print(CD8_dir_files)

  CD4_dir_files <- dir(CD4_dir)
  print("CD4_dir_files")
  print(CD4_dir_files)

  
  output_dir <- toString(args[4])
  
  
} else {

  CD8_dir_files <- toString(args[1])
  print("CD8_dir_files")
  print(CD8_dir_files)

  CD4_dir_files <- toString(args[2])
  print("CD4_dir_files")
  print(CD4_dir_files)

  input_dir <- toString(args[3])

  output_dir <- toString(args[4])

  CD8_dir <- input_dir
  print("CD8_dir")
  print(CD8_dir)
  
  CD4_dir <- input_dir
  print("CD4_dir")
  print(CD4_dir)
}


"%+%" <- function(...){
  paste0(...)
}


output_dir_CD8_cleaned <- output_dir %+% "CD8_cleaned/"
output_dir_CD4_cleaned <- output_dir %+% "CD4_cleaned/"
output_dir_CD8_dropped <- output_dir %+% "CD8_dropped/"
output_dir_CD4_dropped <- output_dir %+% "CD4_dropped/"


if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Dir already exists!")
}

if (!dir.exists(output_dir_CD8_cleaned)){
  dir.create(output_dir_CD8_cleaned)
} else {
  print("Dir 'CD8_cleaned' already exists!")
}

if (!dir.exists(output_dir_CD4_cleaned)){
  dir.create(output_dir_CD4_cleaned)
} else {
  print("Dir 'CD4_cleaned' already exists!")
}

if (!dir.exists(output_dir_CD8_dropped)){
  dir.create(output_dir_CD8_dropped)
} else {
  print("Dir 'CD8_dropped' already exists!")
}

if (!dir.exists(output_dir_CD4_dropped)){
  dir.create(output_dir_CD4_dropped)
} else {
  print("Dir 'CD4_dropped' already exists!")
}



c = 1
for (file_name1 in CD8_dir_files){
  print(file_name1)
  file_name1_paste = unlist(strsplit(file_name1, split='_', fixed=TRUE))
  #print(file_name1_paste)
  if (length(file_name1_paste) == 4){
    print("== 4")
    file_name1_paste = unlist(strsplit(paste0(file_name1_paste[c(1,2,3)], collapse = "_"), split='.', fixed=TRUE))[1]
  }
  if (length(file_name1_paste) == 5){
    print("== 5")
    file_name1_paste = unlist(strsplit(paste0(file_name1_paste[c(1,2,3,5)], collapse = "_"), split='.', fixed=TRUE))[1]
  }
  
  print(file_name1_paste)
  for (file_name2 in CD4_dir_files){
    
    file_name2_paste = unlist(strsplit(file_name2, split='_', fixed=TRUE))
    #print(file_name1_paste)
    if (length(file_name2_paste) == 4){
      print("== 4")
      file_name2_paste = unlist(strsplit(paste0(file_name2_paste[c(1,2,3)], collapse = "_"), split='.', fixed=TRUE))[1]
    }
    if (length(file_name2_paste) == 5){
      print("== 5")
      file_name2_paste = unlist(strsplit(paste0(file_name2_paste[c(1,2,3,5)], collapse = "_"), split='.', fixed=TRUE))[1]
    }
    
    print(file_name2_paste)
    
    if (file_name1_paste == file_name2_paste){
      print("Found")
      print(file_name2)
      print(c)
      c = c + 1
      CD8_table <- read_delim(CD8_dir %+% "/" %+% file_name1, "\t", escape_double = FALSE, trim_ws = TRUE)[,1:7]
      CD4_table <- read_delim(CD4_dir %+% "/" %+% file_name2, "\t", escape_double = FALSE, trim_ws = TRUE)[,1:7]

      CD8_cleaned <- merge(CD8_table, CD4_table, by=c("CDR3nt", "V", "J"), all = TRUE) %>% 
        filter((is.na(frequency.y) & is.na(CDR3aa.y) & is.na(count.y)) | frequency.x>10*frequency.y) %>% 
        select(-c(count.y, frequency.y, CDR3aa.y, D.y))
        colnames(CD8_cleaned) <- c("CDR3nt", "V", "J", "count", "frequency", "CDR3aa", "D")
        
      CD4_cleaned <- merge(CD8_table, CD4_table, by=c("CDR3nt", "V", "J"), all = TRUE) %>% 
        filter((is.na(frequency.x) & is.na(CDR3aa.x) & is.na(count.x)) | frequency.y>10*frequency.x) %>% 
        select(-c(count.x, frequency.x, CDR3aa.x, D.x))
      colnames(CD4_cleaned) <- c("CDR3nt", "V", "J", "count", "frequency", "CDR3aa", "D")
      
      CD8_dropped <- merge(CD8_table, CD4_table, by=c("CDR3nt", "V", "J"), all = TRUE) %>% 
        filter(!((is.na(frequency.y) & is.na(CDR3aa.y) & is.na(count.y)) | frequency.x>10*frequency.y))
      colnames(CD8_dropped) <- c("CDR3nt", "V", "J", "count_CD8", "frequency_CD8", "CDR3aa_CD8", "D_CD8", "count_CD4", "frequency_CD4", "CDR3aa_CD4", "D_CD4")
      
      CD4_dropped <- merge(CD8_table, CD4_table, by=c("CDR3nt", "V", "J"), all = TRUE) %>% 
        filter(!((is.na(frequency.x) & is.na(CDR3aa.x) & is.na(count.x)) | frequency.y>10*frequency.x))
      colnames(CD4_dropped) <- c("CDR3nt", "V", "J", "count_CD8", "frequency_CD8", "CDR3aa_CD8", "D_CD8", "count_CD4", "frequency_CD4", "CDR3aa_CD4", "D_CD4")
      
      fwrite(CD8_cleaned, file = output_dir_CD8_cleaned %+% file_name1_paste %+% "_CD8.csv", sep = "\t", row.names=FALSE, quote = FALSE)
      fwrite(CD4_cleaned, file = output_dir_CD4_cleaned %+% file_name1_paste %+% "_CD4.csv", sep = "\t", row.names=FALSE, quote = FALSE)
      fwrite(CD8_dropped, file = output_dir_CD8_dropped %+% file_name1_paste %+% "_CD8.csv", sep = "\t", row.names=FALSE, quote = FALSE)
      fwrite(CD4_dropped, file = output_dir_CD4_dropped %+% file_name1_paste %+% "_CD4.csv", sep = "\t", row.names=FALSE, quote = FALSE)
    }
  }
}





