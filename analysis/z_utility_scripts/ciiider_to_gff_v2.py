
## This version takes as input a table in which the MA### IDs of target motifs are in the first column
# Usage: python ciiider_to_gff.py [path/to/target tsv] [path/to/scan result csv] [path/to/outfile.gff]

import sys


target_file = sys.argv[1]
scan_file = sys.argv[2]
out_file = sys.argv[3]

target_TFBS = []

with open(target_file) as a:
    for line in a.readlines()[1:]:
        line = line.split('\t')
        target_TFBS.append(line[0])

# print target_TFBS

gff_entries = ['GFF']

with open(scan_file) as a:
    for line in a.readlines()[1:]:
        line = line.split(',')
        if line[3] in target_TFBS:
            strand = '+'
            if line[6] == '-1':
                strand = '-'
            gff_entries.append('\t'.join([line[0],'ciiider',line[2],line[4],line[5],'.',strand,'.',line[2]]))


with open(out_file,'w') as out:
    for line in gff_entries:
        print >> out, line
