library(ggplot2, lib.loc ="/home/rx32940/Rlibs")
library(dplyr, lib.loc ="/home/rx32940/Rlibs")
library(tidyr)

# this is sample file
# file1 <- "/scratch/rx32940/minION/polyA_cDNA/pychopper/output/Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_NoPolyATail/trim_pos_sorted.bed"
# file2 <- "/scratch/rx32940/minION/polyA_cDNA/pychopper/output/Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_NoPolyATail/aln_hits_kept.bed" 

main_path <- "/scratch/rx32940/minION/polyA_cDNA/pychopper/output"

files1 <- list.files(file.path(main_path),recursive=TRUE, pattern= "*trim_pos_sorted.bed", full.name=TRUE)
files2 <- list.files(file.path(main_path), recursive=TRUE,pattern= "*aln_hits_kept.bed", full.name=TRUE)

combined_df <- data.frame(primer=character(), order=character(), avg_trim_len=numeric(), avg_Q=numeric(), num_primers = numeric(), sample=character())
for(i in 1:length(files1)){

file1 <- files1[i]
file2 <- files2[i]

print(file1)

str_ls <- unlist(strsplit(file1, split="/", fixed=TRUE))
sample <- str_ls[length(str_ls)-1]

trim_pos <- read.csv(file1, sep="\t", header=FALSE,col.names=c("readid", "start", "end", "readLen"))
trim_primer <- read.csv(file2, sep="\t", header=FALSE, col.names=c("read_id", "primer_start", "primer_end", "primer", "Q", "strand")) %>% separate(col="read_id", into=c("readid", "others"), extra="merge", sep=" ") %>% select(-others)

trim_primer_pos <- left_join(trim_primer, trim_pos, by="readid") 

trim_len <- trim_primer_pos %>% mutate(trim_length = mapply(function(x,y,z,k,c){
    if( TRUE %in% (c(z,k) %in% c(x,y)) ){ # if cur mapped primer star and end pos prsent in the final trim posi
        if(z %in% c(x, y)){ # if trim at the start of a primer's mapped position, what after the primer was be trimmed
            c - z + 1
        }else{ # if trim at the end of a primer;s map position, what before the primer was trimmed 
            k
        }
    }else{ # this map position was not selected
        NA
    }
}, start, end, primer_start, primer_end, readLen)) %>% subset(!is.na(trim_length))

trim_len.1 <- trim_len %>% mutate(order = mapply(function(x,y,z,k,c){
if(x == k & c == "SSP"){
    "1st"
}else if(y == z & c == "VNP"){
    "1st"
}else{
    "2nd"
}
}, start, end, primer_start, primer_end, primer))

df <- trim_len.1 %>% group_by(primer, order) %>% summarise(avg_trim_len=mean(trim_length), avg_Q=mean(Q), num_primers = n_distinct(readid))
df$sample <- sample

combined_df <- rbind(combined_df, df)
}

write.csv(combined_df, file.path(main_path, "..","aln_pychopper_trimmed_len.csv"), quote=FALSE, row.names=FALSE)