library(dplyr, lib.loc ="/home/rx32940/Rlibs")
library(tidyr, lib.loc ="/home/rx32940/Rlibs")

args <- commandArgs(trailingOnly = TRUE)

tu_path <- args[1] # input
out_path <- args[2]
feature <- args[3] # tss or tts

tabs <-list.files(tu_path, pattern=paste0("*.combined_",feature, ".tab"), full.names = TRUE)

raw_tss_dfs <- lapply(tabs, function(x){read.table(x, sep="\t", header = TRUE) %>% select(!end) %>% mutate(sample=basename(x))})


raw_tss_combined <- Reduce(function(...){full_join(..., all=TRUE,by=c("chrom", "strand", "gene", "start"))},raw_tss_dfs)


all_samples_combined_raw <- raw_tss_combined %>% pivot_longer(cols=starts_with("sample"), names_to="colname", values_to = "files") %>% 
subset(!is.na(files))%>%
 mutate(samples = case_when(
    files == paste0("Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_NoPolyATail_Qiagen_rna_filtered.combined_" , feature,".tab") ~ "LIC-cDNA-nonpolyA_Q",
    files ==paste0("Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_NoPolyATail_rna_filtered.combined_" , feature,".tab") ~ "LIC-cDNA-nonpolyA",
    files ==paste0("Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_PolyATail_rna_filtered.combined_" , feature,".tab")~"LIC-cDNA-polyA",
    files ==paste0("LIC_NOPOLYA_rna_filtered.combined_" , feature,".tab")~ "LIC-dRNA-nonpolyA",
    files ==paste0("LIC_POLYA_DRNA_CDNA_rna_filtered.combined_" , feature,".tab") ~ "LIC-dRNA-cDNA-polyA",
    files ==paste0("LIC_POLYA_rna_filtered.combined_" , feature,".tab") ~ "LIC-dRNA-polyA",
    files ==paste0("Q29_Copenhageni_Basecalled_May_22_2020_Direct-RNA_rna_filtered.combined_" , feature,".tab") ~ "Q29-dRNA-nonpolyA-R1",
    files ==paste0("Q29_Copenhageni_Basecalled-June_11_2020_Repeat_Direct-RNA_rna_filtered.combined_" , feature,".tab") ~ "Q29-dRNA-nonpolyA-R2",
    files ==paste0("Q36_Copenhageni_Basecalled_June_9_2020-Repeat_Direct-RNA_rna_filtered.combined_" , feature,".tab") ~ "Q36-dRNA-nonpolyA-R2",
   files ==paste0("Q36_Copenhageni_Basecalled_May_31_2020_Direct-RNA_rna_filtered.combined_" , feature,".tab") ~ "Q36-dRNA-nonpolyA-R1")) %>%
  select(-c("colname", "files")) %>% mutate(present = "x")  %>% subset(!is.na(samples))%>% pivot_wider(names_from = samples, values_from = present) %>% arrange(start)

chromIleadrawcDNA <- all_samples_combined_raw %>% subset(chrom != "NC_005824.1" & strand == "+") %>% select(!starts_with("Q36")) %>% select(c("chrom","start", "strand", "gene", starts_with("LIC-cDNA"))) %>% subset(!(is.na(`LIC-cDNA-nonpolyA_Q`) & is.na(`LIC-cDNA-nonpolyA`) & is.na(`LIC-cDNA-polyA`)))
chromIlagrawcDNA <- all_samples_combined_raw %>% subset(chrom != "NC_005824.1" & strand == "-") %>% select(!starts_with("Q36")) %>% select(c("chrom","start", "strand", "gene", starts_with("LIC-cDNA")))%>% subset(!(is.na(`LIC-cDNA-nonpolyA_Q`) & is.na(`LIC-cDNA-nonpolyA`) & is.na(`LIC-cDNA-polyA`)))
chromIIleadrawcDNA <- all_samples_combined_raw %>% subset(chrom != "NC_005823.1" & strand == "+") %>% select(!starts_with("Q36")) %>% select(c("chrom","start", "strand", "gene", starts_with("LIC-cDNA")))%>% subset(!(is.na(`LIC-cDNA-nonpolyA_Q`) & is.na(`LIC-cDNA-nonpolyA`) & is.na(`LIC-cDNA-polyA`)))
chromIIlagrawcDNA <- all_samples_combined_raw %>% subset(chrom != "NC_005823.1" & strand == "-") %>% select(!starts_with("Q36")) %>% select(c("chrom","start", "strand", "gene", starts_with("LIC-cDNA")))%>% subset(!(is.na(`LIC-cDNA-nonpolyA_Q`) & is.na(`LIC-cDNA-nonpolyA`) & is.na(`LIC-cDNA-polyA`)))

chromIleadrawdRNA <- all_samples_combined_raw %>% subset(chrom != "NC_005824.1" & strand == "+") %>% select(!starts_with("Q36")) %>% select(c("chrom","start", "strand", "gene", !starts_with("LIC-cDNA"))) %>% subset(!(is.na(`Q29-dRNA-nonpolyA-R1`) & is.na(`LIC-dRNA-nonpolyA`) & is.na(`Q29-dRNA-nonpolyA-R2`) & is.na(`LIC-dRNA-polyA`)))
chromIlagrawdRNA <- all_samples_combined_raw %>% subset(chrom != "NC_005824.1" & strand == "-") %>% select(!starts_with("Q36")) %>% select(c("chrom","start", "strand", "gene", !starts_with("LIC-cDNA")))%>% subset(!(is.na(`Q29-dRNA-nonpolyA-R1`) & is.na(`LIC-dRNA-nonpolyA`) & is.na(`Q29-dRNA-nonpolyA-R2`) & is.na(`LIC-dRNA-polyA`)))
chromIIleadrawdRNA <- all_samples_combined_raw %>% subset(chrom != "NC_005823.1" & strand == "+") %>% select(!starts_with("Q36")) %>% select(c("chrom","start", "strand", "gene", !starts_with("LIC-cDNA")))%>% subset(!(is.na(`Q29-dRNA-nonpolyA-R1`) & is.na(`LIC-dRNA-nonpolyA`) & is.na(`Q29-dRNA-nonpolyA-R2`) & is.na(`LIC-dRNA-polyA`)))
chromIIlagrawdRNA <- all_samples_combined_raw %>% subset(chrom != "NC_005823.1" & strand == "-") %>% select(!starts_with("Q36")) %>% select(c("chrom","start", "strand", "gene", !starts_with("LIC-cDNA")))%>% subset(!(is.na(`Q29-dRNA-nonpolyA-R1`) & is.na(`LIC-dRNA-nonpolyA`) & is.na(`Q29-dRNA-nonpolyA-R2`) & is.na(`LIC-dRNA-polyA`)))

chrom_strand_df <- c("chromIleadrawcDNA", "chromIlagrawcDNA", "chromIIleadrawcDNA", "chromIIlagrawcDNA", "chromIleadrawdRNA", "chromIlagrawdRNA", "chromIIleadrawdRNA", "chromIIlagrawdRNA")

for(df in chrom_strand_df){
  # print(get(df))
  write.csv(get(df), file.path(out_path,paste0(df, "_raw_", toupper(feature), ".csv")), sep=",", quote = FALSE)
  
}