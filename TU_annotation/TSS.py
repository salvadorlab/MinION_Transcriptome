# gene file with the annotated genes, a table where chromosome is in the first column, start codon position in the fourth column, stop codon in the fifth, sign in the seventh, name in the ninth
gene_file = '/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/input/GCF_000007685.1_gene_table_test1.txt'
#input file is a bed file generated from sorted bam alignment file using bedtools bamtobed
input_file = '/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/input/Q29_Copenhageni_Basecalled_May_22_2020_Direct-RNA.sorted.test.bed'
output_file_aware = '/scratch/rx32940/minION/polyA_directRNA/TU_Annotation/input/Q29_Copenhageni_Basecalled_May_22_2020_Direct-RNA.sorted.TSS.test.tab'
#threshold is in general dependent on the sequencing depth
threshold = 10
with open(gene_file, 'r') as gen, open(input_file, 'r') as bed, open(output_file_aware, 'w') as new:
    all_reads_neg = {'NC_005823.1': [], 'NC_005824.1': []}  # both genes and reads in their chromosome
    all_reads_pos = {'NC_005823.1': [], 'NC_005824.1': []}
    genes = {'NC_005823.1': [], 'NC_005824.1': []} #only genes
    for line in gen: #store genes coordinates first
        read = line.strip().split()
        if read[6] == '-':
            all_reads_neg[read[0]].append((int(read[3]), int(read[4]), read[8], 0)) #zero in the end indicates that it is a gene not a read
        else:
            all_reads_pos[read[0]].append((int(read[3]), int(read[4]), read[8], 0))
    for line in bed: #store reads coordinates
        read = line.strip().split()
        if read[5] == '-':
            all_reads_neg[read[0]].append((int(read[1]), int(read[2]), 'read', 1, 1)) #1 to indicate it is a read not a gene
        else:
            all_reads_pos[read[0]].append((int(read[1]), int(read[2]), 'read', 1, 1))
    for ch in all_reads_pos: #process positive reads and genes seperately from negative ones
        chromosome = all_reads_pos[ch]
        chromosome.sort() #sort genes and reads according to their start site
        read_index = 0
        prev_gene = 0
        for read in chromosome: #read here can be either a gene or a read
            if read[-1] == 0: # stop if it is a gene
                before = list(chromosome[0:read_index]) #take the reads before the start codon
                before.sort(key=lambda before: before[1]) #sort them accoding to their end site
                end = before[-1][1] #the furthest end
                if end > read[0]: # if it is after the start codon i.e. read overlaps the start codon
                    after_index = len(before) - 1
                    # take reads until their end bcome less than the start codon (don't overlap the start codon anymore = )
                    while end > read[0] and after_index > 0:
                        after_index -= 1 #move to the next read
                        end = before[after_index][1] #set the end of the next read, it will be checked in the next round of the while function
                    after = list(before[after_index:])# Only the reads that overlap the start codon
                    #now resolve if their multiple peaks
                    current = list()
                    ends = list()
                    cov = len(after) # the coverage at the start codon is the number of overlapping reads
                    after.sort() #sort again according to the start site
                    prev = after[0][0]
                    for end in after: #end is actually the read not just an end
                        if end[0] - prev < threshold: #if the distance between the read start site and the next one is less than the threshold, take in the same cluster
                            current.append(end[0])
                            ends.append(end[1])
                            prev = end[0]
                        else:
                            if len(current) > (0.3 * cov):  # if the read is far from than the next one, finish the current cluster and check it if its count is important
                                ends.sort() # optional step to calculate the end of these start sites
                                tss_end = ends[-1] # the end is the most extrem end
                                new.write('%s\t%i\t%i\t%s\t+\t%i\n' % (ch, current[0] + 1, len(current), read[2], tss_end))
                            current = [end[0]] #start a new cluster with the current read
                            ends = [end[1]]
                            prev = end[0]
                    if len(current) > (0.3 * cov) and len(current) > 2: #to check the last cluster
                        # print(read[2])
                        ends.sort()
                        tss_end = ends[-1]
                        new.write('%s\t%i\t%i\t%s\t+\t%i\n' % (ch, current[0] + 1, len(current), read[2], tss_end))
            read_index += 1
    for ch in all_reads_neg: #similar to above but designed for negative reads
        chromosome = all_reads_neg[ch]
        chromosome.sort(key=lambda chromosome: chromosome[1])
        read_index = 0
        for read in chromosome:
            if read[-1] == 0:
                # after_idx = read_index
                current = list()
                for i in chromosome[read_index + 1:]:
                    if i[0] > read[1]:
                        break
                    elif i[1] > read[1]:
                        current.append((i[1], i[0]))  # i[1]
                if len(current) > 0:
                    current.sort()
                    prev = current[0][0]
                    cov = len(current)
                    curr = list()
                    ends = list()
                    for i in current:
                        if (i[0] - prev) < 10:
                            curr.append(i[0])
                            prev = i[0]
                            ends.append(i[1])
                        else:
                            if len(curr) > (0.3 * cov):
                                ends.sort()
                                tss_end = ends[0]
                                new.write('%s\t%i\t%i\t%s\t-\t%i\n' % (ch, curr[-1] + 1, len(curr), read[2], tss_end))
                            curr = [i[0]]
                            prev = i[0]
                    if len(curr) > (0.3 * cov) and len(curr) > 2:
                        ends.sort()
                        tss_end = ends[0]
                        new.write('%s\t%i\t%i\t%s\t-\t%i\n' % (ch, curr[-1] + 1, len(curr), read[2], tss_end))
            read_index += 1
#Identification of TSS with no annotation
intermediate_file = 'raw_clustering_5_prime_n101.tab'
#input file is a bed file generated from sorted bam alignment file using bedtools bamtobed
input_file = '../../n101.bed'
with open(input_file, 'r') as old, open(intermediate_file, 'w') as new:
    pos_reads = {'NC_004603.1': [], 'NC_004605.1': []} #seperate positive from negative reads
    neg_reads = {'NC_004603.1': [], 'NC_004605.1': []}
    for line in old:
        read = line.strip().split()
        if read[-1] == '+':
            pos_reads[read[0]].append(int(read[1]))
        else:
            neg_reads[read[0]].append(int(read[2]))
    for chromosome in pos_reads:
        chromo = pos_reads[chromosome]
        chromo.sort()
        prev = chromo[0]
        current = [chromo[0]]
        for site in chromo[1:]:
            if site - prev < threshold:
                current.append(site)
                prev = site
            else:
                if len(current) > 3: #if the cluster is more than three reads considere it for further analysis
                    new.write('%s\t%i\t%i\t%i\t+\n' % (chromosome, current[0]+1, len(current), current[-1] - current[0]))
                prev = site
                current = [site]

    for chromosome in neg_reads:
        chromo = neg_reads[chromosome]
        chromo.sort()
        prev = chromo[0]
        current = [chromo[0]]
        for site in chromo[1:]:
            if site - prev < 10:
                current.append(site)
                prev = site
            else:
                if len(current) > 3:
                    new.write('%s\t%i\t%i\t%i\t-\n' % (chromosome, current[-1]+1, len(current), current[-1] - current[0]))
                prev = site
                current = [site]
# seperate pos and negative reads to calculate coverage seperately
pos_bed = input_file[0:-2] + '_pos.bed'
neg_bed = input_file[0:-2] + '_neg.bed'
with open(input_file, 'r') as old, open(pos_bed, 'w') as pos, open(neg_bed, 'w') as neg:
    for line in old:
        read = line.strip().split()
        if read[-1] == '+':
            pos.write(line)
        else:
            neg.write(line)
# filter according to coverage, it was calculated with bedtools genomecoverage
output_file_unaware = 'clustering_5_prime_n101.tab'
pos_cov = 'n101_pos_cov.bed'
neg_cov = 'n101_neg_cov.bed'
with open(intermediate_file, 'r') as old, open(output_file_unaware, 'w') as new, open(pos_cov, 'r') as pos, open(neg_cov, 'r') as neg:
    coverage = {'-': {'NC_004603.1': [], 'NC_004605.1': []}, '+': {'NC_004603.1': [], 'NC_004605.1': []}} #store coverage
    for line in pos:
        read = line.strip().split()
        coverage['+'][read[0]].append(int(read[2]))
    for line in neg:
        read = line.strip().split()
        coverage['-'][read[0]].append(int(read[2]))
    for line in old:
        read = line.strip().split()
        if int(read[2]) > coverage[read[-1]][read[0]][int(read[1])]:  # if cluster meet the coverage criteria leave it, else don't print it
            new.write(line)

#combine the TSSs from the two approachs
input_file = 'n101_all_before_combining.tab' #manually combined table of the two previous files
output_file = 'n101_all_tss.tab'
with open(output_file_aware,'r') as aw, open(output_file_unaware,'r') as unaw, open(input_file,'w') as new:
    all_sites = {'NC_004603.1': [], 'NC_004605.1': []}
    for line in aw:
        read = line.strip().split()
        all_sites[read[0]].append((int(read[1]),line))
    for chr in all_sites:
        chromosome = list(all_sites[chr])
        chromosome.sort()
        for site in chromosome:
            new.write(line)
with open(input_file,'r') as old, open(output_file,'w') as new:
    common = 0
    counts = dict()
    annotated = dict()
    ends = dict()
    for line in old:
        read = line.strip().split('\t')
        if read[3][0:4] == 'gene': #our gene name all start with gene, if it is not the case you can use if not read[3].isdigit()
            name = read[0] + ';' + read[1] + ';' + read[4]
            annotated[name] = read[3]
        else:
            read[1] = str(int(read[1]) +1)
            name = read[0] + ';' + read[1] + ';' + read[4]
        if name in counts:
            common +=1
            if counts[name] < int(read[2]):
                counts[name] = int(read[2])
                ends[name] = int(read[5])
        else:
            counts[name] = int(read[2])
            ends[name] = int(read[5])
    for site in counts:
        end = ends[site]
        if site in annotated:
            name = annotated[site]
        else:
            name = 'novel'
        read = site.split(';')
        new.write('%s\t%i\t%i\t%s\t%s\t%i\n'%(read[0],int(read[1]),counts[site],read[2],name,end))
    print('number of final TSSs is %i'%(common))
# Optional: same as unaware clustering but with ends of the TSSs
with open('../Mapping/n101.bed', 'r') as old, open('clustering_5_prime_n101_test-wend.tab', 'w') as new:
    pos_reads = {'NC_004603.1': [], 'NC_004605.1': []}
    neg_reads = {'NC_004603.1': [], 'NC_004605.1': []}
    for line in old:
        read = line.strip().split()
        if read[-1] == '+':
            pos_reads[read[0]].append((int(read[1]), int(read[2])))
        else:
            neg_reads[read[0]].append((int(read[2]), int(read[1])))
    for chromosome in pos_reads:
        print(chromosome)
        chromo = pos_reads[chromosome]
        chromo.sort()
        prev = chromo[0][0]
        current = [chromo[0][0]]
        ends = [chromo[0][1]]
        for site in chromo[1:]:
            if site[0] - prev < 10:
                current.append(site[0])
                ends.append(site[1])
                prev = site[0]
            else:
                if len(current) > 3:
                    ends.sort()
                    new.write('%s\t%i\t%i\t%i\t%i\t+\n' % (
                    chromosome, current[0], len(current), current[-1] - current[0], ends[-2]))
                prev = site[0]
                current = [site[0]]
                ends = [site[1]]

    for chromosome in neg_reads:
        print(chromosome)
        chromo = neg_reads[chromosome]
        chromo.sort()
        prev = chromo[0][0]
        current = [chromo[0][0]]
        ends = [chromo[0][1]]
        for site in chromo[1:]:
            if site[0] - prev < 10:
                current.append(site[0])
                ends.append(site[1])
                prev = site[0]
            else:
                if len(current) > 3:
                    ends.sort()
                    new.write('%s\t%i\t%i\t%i\t%i\t-\n' % (
                    chromosome, current[-1], len(current), current[-1] - current[0], ends[1]))
                prev = site[0]
                current = [site[0]]
                ends = [site[1]]


