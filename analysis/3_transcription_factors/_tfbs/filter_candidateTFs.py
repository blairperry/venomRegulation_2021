

cand_tfs = []

with open('./candidateTFlist_dnaBinding_08.06.21.txt') as a:
    for line in a.readlines():
        cand_tfs.append(line.rstrip().upper())

cand_motifs = []
cand_motif_ids = []

with open('JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt') as a:
    for i,line in enumerate(a.readlines()):
        line = line.rstrip().split('\t')
        # print line[1].upper()
        if line[1].upper() in cand_tfs:
            cand_motifs.append(line)
            cand_motif_ids.append(line[1].upper())

with open('./candidate_filtering/firstRoundFiltered_JASPAR2020_CORE_NR_08.06.21.txt','w') as out:
    for motif in cand_motifs:
        print >> out, '\t'.join(motif[:2])
        print >> out, '\n'.join(motif[2:])

with open('./candidate_filtering/notFound_NeedToCheck_08.06.21.txt','w') as out:
    for id in cand_tfs:
        if id not in cand_motif_ids:
            print >> out, id
