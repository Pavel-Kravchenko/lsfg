library(readr)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggrepel)


directory="/raid/users/pk/working_dir/mixcr_results_from_KK_2/"
files <- dir(directory)
files

`%+%` <- function(a, b) paste0(a, b)

file_names1 <- grep(".+_8.c.txt", files, value=TRUE)
file_names2 <- grep(".+_8_p..c.txt", files, value=TRUE)
file_names <- append(file_names1, file_names2)
file_names
print(length(file_names)/2)

Patients_common_table <- data.frame(count=NaN, frequency=NaN, CDR3nt=NaN, CDR3aa=NaN, V=NaN, D=NaN, J=NaN, Vend=NaN, Dstart=NaN, Dend=NaN, Jstart=NaN, name=NaN) 
Patients_common_table
p = 0
for (i in file_names){
  print(i)
  p = p + 1
  curr <- read_delim(directory %+% "/" %+% i, "\t", escape_double = FALSE, trim_ws = TRUE)[,1:11]
  name <- rep( paste0(unlist(strsplit(i, split='_', fixed=TRUE))[1:3], collapse = "_"), nrow(curr))
  curr["name"] = name
  Patients_common_table <- rbind(Patients_common_table, curr)
}


print(p)
head(Patients_common_table, 10)
#dim(curr)
#View(table)


K_mers_df <- data.frame(count=NaN, frequency=NaN, CDR3aa=NaN, CDR3nt=NaN, kmer=NaN, V=NaN, name_count=NaN, name=NaN)
K_mers_df

k_mers_table_B27 <- read_csv("~/working_dir/mixcr_results_from_KK_2/total_B27_from_KK_10_k_mers.csv") %>% select(-N_as, -N_psa)
k_mers_table_B27 <- filter(k_mers_table_B27, N > 2)
k_mers_table_B38 <- read_csv("~/working_dir/mixcr_results_from_KK_2/total_B38_from_KK_10_k_mers.csv") %>% select(-N_as, -N_psa)
k_mers_table_B38 <- filter(k_mers_table_B38, N > 2)
k_mers_table_C12 <- read_csv("~/working_dir/mixcr_results_from_KK_2/total_C12_from_KK_10_k_mers.csv") %>% select(-N_as, -N_psa)
k_mers_table_C12 <- filter(k_mers_table_C12, N > 2)


# for (row_i in 1:nrow(k_mers_table)) {
#   kmer_value <- k_mers_table[row_i, "kmer"][[1]]
#   v_value  <- k_mers_table[row_i, "v"][[1]]
#   name_count  <- k_mers_table[row_i, "N"][[1]]
#   # N_B27_value  <- k_mers_table[row_i, "N_B27"][[1]]
#   # N_B27_neg_value  <- k_mers_table[row_i, "N_B27_neg"][[1]]
#   # N_as_value  <- k_mers_table[row_i, "N_as"][[1]]
#   # N_psa_value  <- k_mers_table[row_i, "N_psa"][[1]]
#   print(kmer_value) 
#   #print(v_value)
#   #print(name)
#   tmp <- filter(Patients_common_table, V == v_value) # Можно ускорить, если подавать список сразу на сравнение
#   for (row_j in 1:nrow(tmp)) {
#     CDR3aa_value <- tmp[row_j, "CDR3aa"][[1]]
#     CDR3nt_value <- tmp[row_j, "CDR3nt"][[1]]
#     count_value  <- tmp[row_j, "count"][[1]]
#     frequency_value  <- tmp[row_j, "frequency"][[1]]
#     name_value  <- tmp[row_j, "name"][[1]]
#     #print(CDR3aa_value) 
#     #print(count_value)
#     #print(frequency_value)
#     found_k_mer <- grepl(kmer_value, CDR3aa_value)
#     if (found_k_mer == TRUE){
#       K_mers_df[nrow(K_mers_df) + 1,] = list(count_value, frequency_value, CDR3aa_value, CDR3nt_value, kmer_value, v_value, name_count, name_value)
#     }
#   }
# }

#colnames(tmp)


paired_pb_sf_cd8 <- read_delim("~/working_dir/paired_pb_sf_cd8.txt", "\t", escape_double = FALSE, trim_ws = TRUE)  %>%
  select(status, donor, source, B27, B38, C12) %>%
  mutate(name = apply(paired_pb_sf_cd8[1:3], 1, paste, collapse="_"))

# total <- merge(K_mers_df, paired_pb_sf_cd8, by.x = "name_value", by.y = "name") %>% 
# select("name_value", "count_value", "N", "N_B27", "N_B27_neg", "N_as", "N_psa")
# colnames(total)
# 
# View(total)

#write.csv(total,directory %+% "total_df_with_clones.csv", row.names = FALSE)




total.apply <- apply(X = k_mers_table_B27, MARGIN = 1,function(x){
  
  Patients_common_table %>%
    filter(V==x[2]) %>%
    filter(   grepl(x = Patients_common_table %>% 
                      filter(V==x[2]) %>% 
                      pull(CDR3aa), pattern=x[1])) %>% 
    mutate(kmer=x[1])
}) %>% rbindlist() %>% separate(name, c("status", "donor", "source"), sep = "_") %>% mutate(len=nchar(CDR3aa))

k_mer_positions <- mapply(FUN=function(cdr3,kmer){
  regexpr(kmer,cdr3)[1]
}, total.apply$CDR3aa,total.apply$kmer) 


total.apply$k_mer_position = k_mer_positions

#write.csv(total.apply,directory %+% "total.apply.csv", row.names = FALSE)

PB_SF_freq_with_Na <- dcast(total.apply, CDR3nt+donor+V~source, value.var="frequency",fun.aggregate = mean)
PB_SF_freq <- replace(PB_SF_freq_with_Na, is.na(PB_SF_freq_with_Na), 0)
total.apply_with_freq <- merge(PB_SF_freq, total.apply %>% select(-frequency, -donor, -V, -source), by="CDR3nt") %>%  group_by(donor, kmer, V) %>% summarise(sum_SF=sum(SF),sum_PB=sum(PB))
#total.apply_with_freq %<>% mutate( dif=sum_SF - sum_PB) %>% mutate( sum_f=sum_SF + sum_PB)
#write.csv(total.apply_with_freq,directory %+% "total.apply_with_freq.csv", row.names = FALSE)

total.apply_with_freq %<>% merge(paired_pb_sf_cd8, by="donor")


pseudocount<-Patients_common_table %>% filter(!is.nan(count)) %>% group_by(name) %>% summarise(UMI_sum=sum(count)) %>% mutate(pseudocount=0.1/UMI_sum) %>% separate("name",into=c("condition","donor","source"))

total.apply_with_pseudocount <- merge(total.apply_with_freq,pseudocount %>% dcast(donor~source,value.var="pseudocount"),by="donor") %>% mutate(sum_SF=ifelse(sum_SF==0,SF,sum_SF),sum_PB=ifelse(sum_PB==0,PB,sum_PB) )

total.apply_with_pseudocount %<>% select(-SF,-PB) 
total.apply_with_pseudocount %<>% mutate(fold_change=sum_SF/sum_PB)
write.csv(total.apply_with_pseudocount,directory %+% "total.apply_with_pseudocount.csv", row.names = FALSE)


total.apply_with_fold_change <- total.apply_with_pseudocount %>% group_by(kmer,V) %>% summarise(fold_change_m=mean(fold_change))

fold_change_pos <- total.apply_with_pseudocount %>%
  filter(B27==1)%>% 
  group_by(kmer,V) %>% summarise(fold_change_m=mean(fold_change))

fold_change_neg <- total.apply_with_pseudocount %>%
  filter(B27==0)%>% 
  group_by(kmer,V) %>% summarise(fold_change_m=mean(fold_change))


#View(diff_freq)
fold_change_res <- total.apply_with_pseudocount %>% merge(fold_change_neg, by=c("kmer","V"))  %>% merge(fold_change_pos, by=c("kmer","V"))
# n <- colnames(fold_change_res)
# n[3] = "fold_change_B27_neg"
# n[4] = "fold_change_B27_pos"
# n
# colnames(fold_change_res) <- n

final_sure <- total.apply_with_fold_change %>% select(kmer,V) %>% merge(fold_change_res %>% select(-sum_PB,-sum_SF,-fold_change_m.x,-fold_change_m.y) ,  by = c("kmer","V"))

final_sure %<>% merge(k_mers_table_C12 , by=c("kmer","V")) %>% merge(k_mers_table_B38 , by=c("kmer","V")) %>% merge(k_mers_table_B27 , by=c("kmer","V"))

names <- final_sure %>% group_by(kmer,V,donor) %>% summarise(N=n())



write.csv(final_sure,directory %+% "final_sure.csv", row.names = FALSE)
