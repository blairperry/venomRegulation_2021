
seen = []
lookup_dict = {}

within_TAD = []

total = 946 # number of SE-gene pairs identified

with open('Cvv_SEassocPairs_TADoverlap.sorted.bed') as a:
    for line in a.readlines():
        line = line.rstrip().split('\t')
        if line[8] not in seen:
            seen.append(line[8])
            lookup_dict[line[8]] = line
        else:
            new = [line[3],line[4]]
            other = [lookup_dict[line[8]][3],lookup_dict[line[8]][4]]
            if new[0] == other[0] and new[1] != other[1]:
                out_line = line[0],line[1],line[2],line[3],line[8]
                within_TAD.append(out_line)


print 'Total number of SE-gene pairs: ' + str(total) + '\n'

print 'Number of SE-gene pairs overlapping same TAD: ' + str(len(within_TAD))
print 'Percent of SE-gene pairs overlapping same TAD: ' + str(len(within_TAD)/float(total) * 100) + '\n'

with open('Cvv_SEassocPairs_ContainedWithinTAD_12.10.18.bed','w') as o:
    for line in within_TAD:
        print >> o, '\t'.join(line)
