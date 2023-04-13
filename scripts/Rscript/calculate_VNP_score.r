library(ggplot2, lib.loc ="/home/rx32940/Rlibs")
library(dplyr, lib.loc ="/home/rx32940/Rlibs")


main_path="/scratch/rx32940/minION/polyA_cDNA/pychopper/output"

vnp_files <- list.files(file.path(main_path), pattern = "aln_hits_kept.bed", recursive = TRUE, full.names=TRUE)
file <- vnp_files[1]

read_aln_file <- function(file){
    pd <- read.csv(file, sep="\t", header=FALSE) 
    file_ls <- unlist(strsplit(file, split="/", fixed=TRUE))
    pd$file <- file_ls[length(file_ls)-1]
    pd
}
combined_vnp <- do.call(rbind,lapply(vnp_files, read_aln_file))

aln_summary <- combined_vnp %>% group_by(file, V4) %>% summarise(avg_score = mean(V5), number_distinct=n_distinct())

write.csv(aln_summary, file.path(main_path, "../aln_primers_summary.csv"))