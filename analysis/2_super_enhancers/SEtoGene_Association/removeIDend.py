with open('Cvv_SE_NearestTSS_simpleInfo.bed','w') as o:
    with open('Cvv_SE_NearestTSS.bed') as a:
        for line in a.readlines():
            line = line.rstrip().split('\t')
            info = line[7]
            new_info = ';'.join(info.split(';')[:-1])
            print >> o, '\t'.join([line[0],line[1],line[2],line[3],line[4],line[5],line[6],new_info,line[8]])
