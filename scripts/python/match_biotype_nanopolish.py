import pandas as pd
import sys
from tqdm import tqdm

path = "/scratch/rx32940/minION/polyA_cDNA/"

# biotype position file processed from gtf in minION polyA analysis
gtf_pos = pd.read_csv(path + "biotype_position.csv")

sample_name=sys.argv[1]

# nanopolish polya output (input file for this script)
polya_output = pd.read_csv(path + 
"nanopolish_output/output_pass/" + sample_name + ".nanopolish.pass.tsv",
 sep="\t",header=None)

polya_pos = polya_output[2]

biotype_list=[]
for i in tqdm(polya_pos):
    current_read_biotype="Unknown"
    for j in range(len(gtf_pos)):
        if i >= gtf_pos.loc[j,"start"] and i <= gtf_pos.loc[j,"end"]:
            current_read_biotype=gtf_pos.loc[j,"biotype"]
            break
    biotype_list.append(current_read_biotype)


polya_output['biotype']=biotype_list

polya_output.to_csv(path+ "nanopolish_output/biotype_annot/" + sample_name + ".nanopolish.pass.biotype.tsv", sep="\t")