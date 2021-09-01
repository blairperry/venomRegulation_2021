
import sys

infile = sys.argv[1]
outdir = sys.argv[2]

if outdir[-1] != '/':
	outdir = outdir + '/'

outfile = outdir + infile[:-4].split('/')[-1] + '.bed'

with open(infile) as a, open(outfile,'w') as out:
	for line in a.readlines():
		line = line.rstrip().split('\t')
		start = min([int(line[8]),int(line[9])])
		end = max([int(line[8]),int(line[9])])
		outline = [line[1],str(start),str(end),line[0]]
		print >> out, '\t'.join(outline)