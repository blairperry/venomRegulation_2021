
TADs = []

with open('../../bed/allChrom_10kb_HierTADS_hier1.bed') as a:
    for line in a.readlines():
        line = line.rstrip().split('\t')
        TADs.append(line)

overlap = []

with open('../Cvv_SE_NearestTSS.bed') as a:
    for line in a.readlines():
        line = line.rstrip().split('\t')
        scaff = line[0]
        start_SE = int(line[1])
        end_SE = int(line[2])
        start_gene = int(min(line[5],line[6]))
        end_gene = int(max((line[5],line[6])))
        for TAD in TADs:
            if TAD[0] == scaff:
                if (int(TAD[1]) <= start_SE and end_SE <= int(TAD[2])) and (int(TAD[1]) <= start_gene and end_gene <= int(TAD[2])):
                    overlap.append(line)
                    print line
                    print scaff, start_SE, end_SE, start_gene, end_gene
                    print TAD[1], TAD[2]
                    break


zeros = 0
nonzeros = 0
total = len(overlap)

with open('Cvv_SEassocGenes_withinTADs.txt','w') as o:
    for line in overlap:
        if line[-1] == '0':
            zeros = zeros + 1
        elif line[-1] != '0':
            nonzeros = nonzeros + 1
        print >> o, '\t'.join(line)

#print 'Number of non-overlapping SE-Gene associations contained within a TAD: ' + str(nonzeros)
