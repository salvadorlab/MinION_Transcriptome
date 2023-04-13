import sys
import re
from Bio import SeqIO

nucl=sys.argv[1]

with open(sys.argv[4], "w") as w:
    RE_STR="{" + sys.argv[2] + ",}"
    pattern = re.compile(nucl + RE_STR)
    for record in SeqIO.parse(sys.argv[3], "fasta"):
        for match in pattern.finditer(str(record.seq).upper()):
            inter_start=0
            inter_end=0
            if match.start() - 0 >= 100:
                inter_start = match.start()-100
            else:
                inter_start = 0
            if len(record) - match.end() >= 100:
                inter_end = match.end()+100
            else:
                inter_end = len(record)-1
                
                
            w.write("\t".join([record.id, str(inter_start), str(inter_end), str(match.start())+":"+str(match.end()-1), str(match.end() - match.start())]) + "\n")
