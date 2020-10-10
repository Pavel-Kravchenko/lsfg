library("ggplot2")

source("/raid/users/pk/kmer_dig/kmer_full_withV.R")

directory="/raid/users/pk/working_dir/mixcr_results/export/"
k_mer_files <- dir(directory)
k_mer_files

`%+%` <- function(a, b) paste0(a, b)

name1 <- grep("as_.+_PB_._k_mers.txt", k_mer_files, value=TRUE)
name1
name2 <- grep("as_.+_SF_._k_mers.txt", k_mer_files, value=TRUE)
name2

for (i in name1){
  for (j in name2){
    one=unlist(strsplit(i, split='_', fixed=TRUE))[-3]
    two=str_c(unlist(strsplit(j, split='_', fixed=TRUE))[-3], sep = "_")
    one=paste(one, collapse="_")
    two=paste(two, collapse="_")
    if(one == two){
      print(i)
      print(i)
      target = directory %+% j
      db = directory %+% i
      
      Patient<-count_kmers(target = target, db = db, with.gaps = F,  target_preCalculated=T, db_precalculated=T)
      Patient
      data.sf<-fread(target)
      data.pb<-fread(db)
      names<-strrep(unlist(strsplit(i, split='_', fixed=TRUE))[2], nrow(Patient))
      Patient["name"] = names
    }
  }
}

View(Patient)
# head(Patient)

write.table(Patient, directory %+% "file.csv")
