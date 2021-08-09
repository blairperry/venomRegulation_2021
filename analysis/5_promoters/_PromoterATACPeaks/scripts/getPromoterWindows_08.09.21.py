
with open('../PromWindows_1kbUpDown_08.09.21.txt','w') as out:
    with open('/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation_NEW_Aug2021/data/annotation/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019.GeneEntriesOnly.gff') as a:
        for line in a.readlines():
            line = line.rstrip().split('\t')
            if 'trnascan' not in line[8]:

                if line[6] == '+':
                    gene_start = int(line[3])
                    prom_start = gene_start - 1000
                    prom_end = gene_start + 1000

                    if prom_start < 0:
                        prom_start = 0

                    desc = '@'.join(['#'+line[6],line[8]])

                    outlist = [line[0], str(prom_start),str(prom_end),desc]

                    print >> out, '\t'.join(outlist)

                else:
                    gene_start = int(line[4])
                    prom_start = gene_start - 1000
                    prom_end = gene_start + 1000

                    if prom_start < 0:
                        prom_start = 0

                    desc = '@'.join(['#'+line[6],line[8]])

                    outlist = [line[0], str(prom_start),str(prom_end),desc]

                    print >> out, '\t'.join(outlist)
