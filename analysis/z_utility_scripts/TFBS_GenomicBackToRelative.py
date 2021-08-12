
import sys

infile = sys.argv[1]

outfile = infile[:-4] + '.relPos.bed'

with open(infile) as a, open(outfile,'w') as out:
	for line in a.readlines():
		line = line.rstrip().split('\t')
		relInfo = line[3].split(';')
		relInfo.append(line[4])
		if len(line) > 5:
			relInfo.append(line[5])
			relInfo.append(line[6])
		print >> out, '\t'.join(relInfo)