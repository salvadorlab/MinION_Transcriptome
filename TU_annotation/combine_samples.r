library(dplyr, lib.loc ="/home/rx32940/Rlibs")
library(tidyr, lib.loc ="/home/rx32940/Rlibs")

args <- commandArgs(trailingOnly = TRUE)

tu_path <- args[1] # input
out_path <- args[2]
feature <- args[3] # tss or tts

# tu_path <- "/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/direct_output/tss"
# out_path <- "/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/direct_output/tss"
# feature <- "TSS"

tabs <-list.files(tu_path, pattern=paste0(toupper(feature), ".merged.tab"), full.names = TRUE)

raw_tss_dfs <- lapply(tabs, function(x){
  short_file <- basename(x)
  short_file_name <- unlist(strsplit(short_file, split=".", fixed = TRUE))[1]
  file <- read.csv(x, sep="\t", header = FALSE) 
  colnames(file) <- c("chrom", "start", "end", "gene", "cov", "strand")
  file$sample <- short_file_name 
  file <- file %>% mutate(sample=
  case_when(
    sample =="LIC_NOPOLYA_rna_filtered"~ "LIC-dRNA-nonpolyA",
    sample =="LIC_POLYA_DRNA_CDNA_rna_filtered" ~ "LIC-dRNA-cDNA-polyA",
    sample =="LIC_POLYA_rna_filtered" ~ "LIC-dRNA-polyA",
    sample =="Q29_Copenhageni_Basecalled_May_22_2020_Direct-RNA_rna_filtered" ~ "Q29-dRNA-nonpolyA-R1",
    sample =="Q29_Copenhageni_Basecalled-June_11_2020_Repeat_Direct-RNA_rna_filtered" ~ "Q29-dRNA-nonpolyA-R2",
    sample =="Q36_Copenhageni_Basecalled_June_9_2020-Repeat_Direct-RNA_rna_filtered" ~ "Q36-dRNA-nonpolyA-R2",
   sample =="Q36_Copenhageni_Basecalled_May_31_2020_Direct-RNA_rna_filtered" ~ "Q36-dRNA-nonpolyA-R1")
)
  file
  })


dRNA_df_combined <- do.call(rbind, raw_tss_dfs)
dRNA_df_combined %>% group_by(chrom, strand) %>% arrange()

dRNA_df_combined$present <- paste(dRNA_df_combined$end, dRNA_df_combined$cov, sep=";")

all_samples_combined_raw <- dRNA_df_combined %>% select(!c(end, cov)) %>% pivot_wider(names_from="sample", values_from = "present") %>% arrange(start)

chromIlead <- all_samples_combined_raw %>% subset(chrom != "NC_005824.1" & strand == "+")  
chromIlag <- all_samples_combined_raw %>% subset(chrom != "NC_005824.1" & strand == "-")
chromIIlead <- all_samples_combined_raw %>% subset(chrom != "NC_005823.1" & strand == "+")  
chromIIlag <- all_samples_combined_raw %>% subset(chrom != "NC_005823.1" & strand == "-")  

chrom_strand_df <- c("chromIlead", "chromIlag", "chromIIlead", "chromIIlag")

for(df in chrom_strand_df){
  # print(get(df))
  write.table(get(df), file.path(out_path,paste0(df, ".combined.", toupper(feature), ".tab")), sep="\t", quote = FALSE, row.names=FALSE)
 
}