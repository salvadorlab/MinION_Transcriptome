import os


tu_path = "/mnt/d/Dropbox/5.Rachel-projects/Transcriptomics_PolyATail/tu_annotation"

TSS_path = os.path.join(tu_path, "TSS")
TTS_path = os.path.join(tu_path, "TTS")


with open(os.path.join(TSS_path, "output", "raw_TSS_clustered", "chromIleadrawcDNA_raw_TSS_clustered_10.csv")) as tss_all, \
    open(os.path.join(TTS_path, "combined_output", "chromIleadrawcDNA_raw_TTS_clustered_10.csv")) as tts_all:
        
        tss_df = tss_all.readlines()  # list
        tts_df = tts_all.readlines()
        cur_chrom = tss_df[1].split(",")[1] 
        cur_strand = tss_df[1].split(",")[3]

        tss_index = 1 # start reading from 2nd line, 1st line is header

with open(os.path.join(tu_path, "ORF", "chromIleadrawcDNA_orf_10.csv"), "w") as orf:
    for tts in tts_df[1:]:
        cur_tts = int(tts.split(",")[2]) # position of tts
        for tss in tss_df[tss_index:]: # check tss
            cur_tss = int(tss.split(",")[2]) # position of tss
            cur_gene = tss.split(",")[4] if tss.split(",")[4] != "" else "unaware"
            if cur_tss < cur_tts: # is curr tss < cur tts, new orf 
                orf_len = cur_tts - cur_tss
                new_orf = cur_chrom + "\t" + str(cur_tss) + "\t" + str(cur_tts) + "\t" + cur_gene + "\t" + str(orf_len + 1) + "\t" + cur_strand + "\n" 
                tss_index += 1
                a = orf.write(new_orf)

            

