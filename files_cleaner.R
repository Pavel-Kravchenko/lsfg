library(readr)
library(dplyr)
library(data.table)
library(rlist)

args <- commandArgs(trailingOnly = T)

first_dir <- toString(args[1]) # first dir
#first_dir <- "/raid/users/pk/working_dir/CD8_seq/"
print(first_dir)

second_dir <- toString(args[2]) # second dir
#second_dir <- "/raid/users/pk/working_dir/CD4_seq/"
print(second_dir)

output_dir <- toString(args[3])
#output_dir <- "/raid/users/pk/working_dir/result/"

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


first_dir_files <- dir(first_dir)
print("first_dir_files")
print(first_dir_files)
first_dir_names <- strsplit(first_dir_files, split='_', fixed=TRUE)
#print(first_dir_names)

second_dir_files <- dir(second_dir)
print("second_dir_files")
print(second_dir_files)
second_dir_names <- strsplit(second_dir_files, split='_', fixed=TRUE)
#print(second_dir_names)


#name_list = list()
c = 1
for (file_name1 in first_dir_files){
  print(file_name1)
  file_name1_paste = paste0(unlist(strsplit(file_name1, split='_', fixed=TRUE))[1:3], collapse = "_")
  #print(file_name1_paste)
  for (file_name2 in second_dir_files){
    file_name2_paste = paste0(unlist(strsplit(file_name2, split='_', fixed=TRUE))[1:3], collapse = "_")
    #print(file_name2_paste)
    if (file_name1_paste == file_name2_paste){
      print("Found")
      print(file_name2)
      print(c)
      c = c + 1
      file1_table <- read_delim(first_dir %+% "/" %+% file_name1, "\t", escape_double = FALSE, trim_ws = TRUE)[,1:7]
      #file1_table_D <- file1_table %>% select(-"D")
      file2_table <- read_delim(second_dir %+% "/" %+% file_name2, "\t", escape_double = FALSE, trim_ws = TRUE)[,1:7]
      #file2_table_D <- file2_table %>% select(-"D")
      #common <- merge(file1_table, file2_table, by=c("CDR3nt", "V", "J"), all = TRUE)
      
      
      CD8_cleaned <- merge(file1_table, file2_table, by=c("CDR3nt", "V", "J"), all = TRUE) %>% 
        filter(is.na("frequency.y") == TRUE & is.na("CDR3nt.y") == TRUE & is.na("count.y") == TRUE | frequency.x>10*frequency.y) %>% 
        select(-c(count.y, frequency.y, CDR3aa.y, D.y))
        colnames(CD8_cleaned) <- c("CDR3nt", "V", "J", "count", "frequency", "CDR3aa", "D")
      CD4_cleaned <- merge(file1_table, file2_table, by=c("CDR3nt", "V", "J"), all = TRUE) %>% 
        filter(is.na("frequency.x") == TRUE & is.na("CDR3nt.x") == TRUE & is.na("count.x") == TRUE | frequency.y>10*frequency.x) %>% 
        select(-c(count.y, frequency.y, CDR3aa.y, D.y))
      colnames(CD4_cleaned) <- c("CDR3nt", "V", "J", "count", "frequency", "CDR3aa", "D")
      CD8_dropped <- merge(file1_table, file2_table, by=c("CDR3nt", "V", "J"), all = TRUE) %>% 
        filter(!(is.na("frequency.y") == TRUE & is.na("CDR3nt.y") == TRUE & is.na("count.y") == TRUE | frequency.x>10*frequency.y))
      colnames(CD8_dropped) <- c("CDR3nt", "V", "J", "count_CD8", "frequency_CD8", "CDR3aa_CD8", "D_CD8", "count_CD4", "frequency_CD4", "CDR3aa_CD4", "D_CD4")
      CD4_dropped <- merge(file1_table, file2_table, by=c("CDR3nt", "V", "J"), all = TRUE) %>% 
        filter(!(is.na("frequency.x") == TRUE & is.na("CDR3nt.x") == TRUE & is.na("count.x") == TRUE | frequency.y>10*frequency.x))
      colnames(CD4_dropped) <- c("CDR3nt", "V", "J", "count_CD8", "frequency_CD8", "CDR3aa_CD8", "D_CD8", "count_CD4", "frequency_CD4", "CDR3aa_CD4", "D_CD4")
      write.csv(CD8_cleaned, file = output_dir_CD8_cleaned %+% file_name1_paste %+% "_CD8_cleaned.csv")
      write.csv(CD4_cleaned, file = output_dir_CD4_cleaned %+% file_name1_paste %+% "_CD4_cleaned.csv")
      write.csv(CD8_dropped, file = output_dir_CD8_dropped %+% file_name1_paste %+% "_CD8_dropped.csv")
      write.csv(CD4_dropped, file = output_dir_CD4_dropped %+% file_name1_paste %+% "_CD4_dropped.csv")
    }
  }
}






