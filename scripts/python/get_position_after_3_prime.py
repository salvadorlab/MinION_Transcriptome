import sys

####################################################################
#
# This script is for getting positions of 30 bps after where the reads were mapped at 3' end
#
####################################################################
bed=sys.argv[1]
after_3prime=sys.argv[2]
before_3prime=sys.argv[3]

# bed="/scratch/rx32940/minION/polyA_cDNA/vnp/out/Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_NoPolyATail/vnp_read_mapped.bed"
# after_3prime="/scratch/rx32940/minION/polyA_cDNA/vnp/out/Copenhageni_Basecalled_Aug_16_2019_Direct-cDNA_NoPolyATail/postion_after_3prime.bed"

with open(bed) as input, open(after_3prime, "w+") as out:
    for line in input:
        ln_ls=line.strip("\n").split("\t")
        strand=ln_ls[5]
        if strand == "+":
            start=int(ln_ls[2])
            end=int(ln_ls[2]) + 30
            a=out.write("%s\t%d\t%d\t%s\n"%(ln_ls[0], start, end, "\t".join(ln_ls[3:])))
        elif strand == "-": # if read mapped to the lag strand
            start=int(ln_ls[1])
            end=int(ln_ls[1]) - 30
            a=out.write("%s\t%d\t%d\t%s\n"%(ln_ls[0], end, start, "\t".join(ln_ls[3:])))

with open(bed) as input, open(before_3prime, "w+") as out:
    for line in input:
        ln_ls=line.strip("\n").split("\t")
        strand=ln_ls[5]
        if strand == "+":
            start=int(ln_ls[2]) -30
            end=int(ln_ls[2]) 
            a=out.write("%s\t%d\t%d\t%s\n"%(ln_ls[0], start, end, "\t".join(ln_ls[3:])))
        elif strand == "-": # if read mapped to the lag strand
            start=int(ln_ls[1]) + 30
            end=int(ln_ls[1]) 
            a=out.write("%s\t%d\t%d\t%s\n"%(ln_ls[0], end, start, "\t".join(ln_ls[3:])))
            