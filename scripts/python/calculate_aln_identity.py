import sys

aln_bed=sys.argv[1]
aln_out=sys.argv[2]

# aln_bed="/scratch/rx32940/minION/polyA_cDNA/map/genome/bed/Mankarso_Basecalled_Aug_16_2019_Direct-cDNA_NoPolyATail.bed"
# aln_out="/scratch/rx32940/minION/polyA_cDNA/map/genome/stats/Mankarso_Basecalled_Aug_16_2019_Direct-cDNA_NoPolyATail.stats"

with open(aln_bed) as aln , open(aln_out, "w+") as out:
    for line in aln:
        cur_dict={"M":[], "I":[]}
        ln_ls=line.split("\t")
        readid=ln_ls[3]
        md = ln_ls[6]
        cur_num=""
        for i in range(len(md)-1):
            if md[i].isdigit():
                cur_num=cur_num + md[i]
            else:
                if md[i] not in cur_dict.keys():
                    cur_dict[md[i]] = []
                cur_dict[md[i]].append(int(cur_num))
                cur_num=""
        align_len = sum(cur_dict["M"]) + sum(cur_dict["I"])
        map_ident = (1-(int(ln_ls[4])/align_len)) * 100
        if "N" in cur_dict:
            stat=str(align_len) + "|" + str(int(map_ident)) + "|" + str(len(cur_dict["N"])) + "|" + str(sum(map(int, cur_dict["N"]))) 
        else:
            stat=str(align_len) + "|" + str(int(map_ident)) + "|" + str(0) + "|" + str(0)
        a= out.write("%s\t%s\n"%(line.strip("\n"),stat))
        
