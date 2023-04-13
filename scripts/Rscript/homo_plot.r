
library(ggplot2, lib.loc ="/home/rx32940/Rlibs")
library(dplyr, lib.loc ="/home/rx32940/Rlibs")
# library(plotly,lib.loc ="/home/rx32940/Rlibs")
 library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
homo.path <- args[1] #"/scratch/rx32940/minION/homopolymer_cov/output"
homo.sub <- args[2] #"homoA"
print(homo.sub)

# homo.path <-"/scratch/rx32940/minION/homopolymer_cov_genome/output/nuc_5"
# homo.sub <- "homoA"

# files_drna <- list.files(file.path(homo.path), pattern = "*PolyATail.[A-Z][0-9].CDS.relative.cov",full.names = TRUE, recursive = TRUE)
files_cdna <- list.files(file.path(homo.path), pattern = "*PolyATail.*.CDS.relative.cov",full.names = TRUE, recursive = TRUE)
files_cdna <- files_cdna[!grepl("Mankarso", files_cdna)]
files_drna <- list.files(file.path(homo.path), pattern = "*.relative.CDS.cov",full.names = TRUE, recursive = TRUE)

all_files <- c(files_cdna,files_drna)

# coverage from both strand
files <- all_files[((grepl(homo.sub, all_files)) & (!grepl("*.ss.*",all_files)) & (!grepl("Qiagen",all_files))) | (grepl("A5", all_files) & grepl("Q29", all_files)) ]
# # coverage from single strand
# files <- all_files[((grepl(homo.sub, all_files)) & (grepl("A5.ss",all_files)|grepl("*.FL.*",all_files) ) & (!grepl("Qiagen",all_files))) | (grepl("A5", all_files) & grepl("Q29", all_files)) ]

# files <- c(files_drna)[(grepl(homo.sub, files_drna) | grepl("homoT", files_drna)) & grepl("Q29", files_drna)]


# file <- "/scratch/rx32940/minION/homopolymer_cov/output/homoA/Icterohaemorrhagiae_Basecalled_Aug_16_2019_Direct-cDNA_NoPolyATail.A3.relative.cov"
all_sample_homo_cov<- data.frame(fromHomePos=numeric(), avgCov= numeric(), SD=numeric(),sample=character(), platform = character(), polyAtail = character(), cur_homo=character(), Data.Type=character())
for(file in files){

  cur_file <- read.csv(file, sep="\t")
  path_ls <- unlist(strsplit(file, split="/", fixed=TRUE))
  cur_homo <- path_ls[length(path_ls)-1]
  
  str_ls <- unlist(strsplit(basename(file), split = ".", fixed = TRUE))
  sample <- str_ls[1] # get sample name
  Data.Type <- ifelse(str_ls[3] == "FL", "FL.cDNA", ifelse(str_ls[3] == "relative", "DRS", "Raw.cDNA"))
  print(sample)
  sample_name <- switch(
    sample,
    `fast5_pass-LIC-NO-POLYA` = "LIC-dRNA-nonpolyA",
    `fast5_passLIC-POLYA-OPT-CDNA` = "LIC-dRNA-cDNA-polyA",
    `fast5_pass-LIC-PolyA` = "LIC-dRNA-polyA",
    `Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_NoPolyATail_Qiagen` = "LIC-cDNA-nonpolyA_Q",
    `Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_NoPolyATail` = "LIC-cDNA-nonpolyA",
    `Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_PolyATail` = "LIC-cDNA-polyA",
    `Icterohaemorrhagiae_Basecalled_Aug_16_2019_Direct-cDNA_NoPolyATail` = "LII-cDNA-nonpolyA",
    `Icterohaemorrhagiae_Basecalled_Aug_16_2019_Direct-cDNA_PolyATail` = "LII-cDNA-polyA",
    `Mankarso_Basecalled_Aug_16_2019_Direct-cDNA_NoPolyATail` = "LIM-cDNA-nonpolyA",
    `Mankarso_Basecalled_Aug_16_2019_Direct-cDNA_PolyATail` = "LIM-cDNA-polyA",
    `Patoc_Basecalled_Aug_16_2019_Direct-cDNA-NoPolyATail` = "LBP-cDNA-nonpolyA",
    `Patoc_Basecalled_Aug_16_2019_Direct-cDNA_PolyATail` = "LBP-cDNA-polyA",
    `Q29_Copenhageni_Direct-RNA_FAST5_Pass`="LIC-DRS-nonpolyA-R1",
    `Q29_Copenhageni_Direct-RNA_Repeat_FAST5_Pass`="LIC-DRS-nonpolyA-R2",
    `Q36_Copenhageni_Direct-RNA_Repeat_FAST5_Pass`="Q36-dRNA-nonpolyA-R2",
    `Q36_Copenhageni_Direct-RNA_FAST5_Pass`="Q36-dRNA-nonpolyA-R1"
  )
  
  # sequencing protocol, dRNA or cDNA
  platform <- unlist(strsplit(sample_name, split="-", fixed=TRUE))[2]
  
  # polyA tail added?
  polyAtail <- ifelse(unlist(strsplit(sample_name, split="-", fixed=TRUE))[3] == "nonpolyA", "nonpolyA", "polyA")
  
  # get the average cov across all transcripts mapped at the position flanking homopolymer region
  avg_pos_cov <- cur_file %>% group_by(fromHomePos) %>% summarise(avgCov = mean(relative_cov, na.rm=TRUE), SD=sd(relative_cov, na.rm=TRUE)) %>% mutate(sample=sample_name) %>% mutate(platform = platform, polyAtail = polyAtail, cur_homo = cur_homo, Data.Type = Data.Type)
  
  all_sample_homo_cov <- rbind(all_sample_homo_cov,avg_pos_cov)
}
# unique(all_sample_homo_cov$sample)

write.csv(all_sample_homo_cov, file.path(homo.path, paste0("combined_all_samples.",homo.sub,".csv")), quote = FALSE, row.names = FALSE)
# write.csv(all_sample_homo_cov, file.path(homo.path, paste0("combined_all_samples.", homo.sub,".rawSS.csv")), quote = FALSE, row.names = FALSE)

all_sample_homo_cov.1 <- all_sample_homo_cov %>% separate(sample, into=c("serovars", "protocol", "ispolyA"), remove=FALSE ) %>% group_by(Data.Type, cur_homo,fromHomePos,ispolyA ) %>% summarise(Avg.cov=mean(avgCov), sd=sd(avgCov)) %>% mutate(sd_show = ifelse(fromHomePos%%30==0 , sd, NA))

all_sample_homo_cov.1$Data.Type <- factor(all_sample_homo_cov.1$Data.Type, levels=c("Raw.cDNA", "FL.cDNA", "DRS"))
p <- ggplot(all_sample_homo_cov.1, aes(x=fromHomePos, y=Avg.cov,color = Data.Type))+
  facet_grid(~ispolyA, scales="free_y")+
      geom_errorbar(aes(ymin=`Avg.cov`-sd_show, ymax=`Avg.cov`+sd_show), width=5, alpha=0.3) +
  geom_line(size=1,alpha=0.6)+
  scale_x_continuous(breaks = seq(-100,100,25), limits = c(-100,100))+
  labs(x="5'-3'position relative to homopolymers of A's (bps)", y="Relative coverage of the nucleotide position")+
  theme(text = element_text(size = 10))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_color_manual(values=RColorBrewer::brewer.pal(12, "Dark2"))

# # ggplotly(p)
ggsave(paste0("/scratch/rx32940/minION/homopolymer_cov_genome/output/nuc_5/nuc_5_", homo.sub,"_CDS.pdf"),p,width = 8,height =4, dpi=500)
# ggsave("/scratch/rx32940/minION/homopolymer_cov_genome/output/nuc_5/nuc_5_CDS.ss.pdf",p,width = 8,height =4, dpi=500)

# test <- ggplot(all_sample_homo_cov%>% separate(sample, into=c("serovars", "protocol", "ispolyA"), remove=FALSE ) %>% subset(Data.Type %in% c("FL.cDNA", "DRS")), aes(x=fromHomePos, y=avgCov,color = sample))+
#   facet_grid(~ispolyA, scales="free_y")+
#       # geom_errorbar(aes(ymin=`Avg.cov`-sd_show, ymax=`Avg.cov`+sd_show), width=5, alpha=0.3) +
#   geom_line(size=1,alpha=0.6)+
#   scale_x_continuous(breaks = seq(-100,100,25), limits = c(-100,100))+
#   labs(x="5'-3'position relative to homopolymers of A's (bps)", y="Relative coverage of the nucleotide position")+
#   theme(text = element_text(size = 10))+
#   theme_bw()+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#     scale_color_manual(values=RColorBrewer::brewer.pal(12, "Paired"))