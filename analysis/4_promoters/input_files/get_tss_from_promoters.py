import sys

infile = sys.argv[1]

outfile = infile[:-4] + '.tss.bed'

with open(outfile,'w') as out:
    with open(infile) as a:
        for line in a.readlines():
            line = line.rstrip().split()
            if 'trnascan' not in line[3]:
                if line[5] == '+':
                    print >> out, '\t'.join([line[0], line[1], str(int(line[1])+1),line[3]])
                elif line[5] == '-':
                    print >> out, '\t'.join([line[0], line[2], str(int(line[2])+1),line[3]])
