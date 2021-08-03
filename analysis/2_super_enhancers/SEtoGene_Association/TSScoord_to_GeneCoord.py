
with open('Cvv_SE_NearestGeneCoords.bed','w') as o:
    with open('Cvv_SE_NearestTSS.bed') as a:
        for line in a.readlines():
            line = line.rstrip().split('\t')
            gene_start = ((line[7].split(';')[-1]).split('_'))[0]
            gene_end = ((line[7].split(';')[-1]).split('_'))[1]
            new_desc = ';'.join((line[7].split(';')[:-1]))
            print >> o, '\t'.join([line[0],line[1],line[2],line[3],line[4],gene_start,gene_end,new_desc,line[8]])
