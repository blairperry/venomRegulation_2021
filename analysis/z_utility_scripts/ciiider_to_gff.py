
## Example files to build the script with:
##    - SVMP vs. NonVenom Genes Enrichment: '4_VenomFamilies.vs.NonVenGenes/SVMP/Enrichment: Text1360239721433099872_MostSigDeficit.csv'
##    - TFBS scan of SVMPs + paralogs: '/Users/perryb/Desktop/ciiider_temp/Promoter Panel.csv'

# Usage: python ciiider_to_gff.py [path/to/enrichment csv] [path/to/scan result csv] [path/to/outfile.gff]

import sys

enrich_file = sys.argv[1]
scan_file = sys.argv[2]
out_file = sys.argv[3]

# enrich_file = '4_VenomFamilies.vs.NonVenGenes/SVMP/Enrichment: Text1360239721433099872_MostSigDeficit.csv'
# scan_file = '/Users/perryb/Desktop/ciiider_temp/Promoter Panel.csv'
# out_file = 'test.gff'

enrich_TFBS = []

with open(enrich_file) as a:
    for line in a.readlines()[1:]:
        line = line.split(',')
        if float(line[15]) < 0.05 and line[14] == 'Up':
            enrich_TFBS.append(line[1])

gff_entries = ['GFF']

with open(scan_file) as a:
    for line in a.readlines()[1:]:
        line = line.split(',')
        if line[2] in enrich_TFBS:
            strand = '+'
            if line[6] == '-1':
                strand = '-'
            gff_entries.append('\t'.join([line[0],'ciiider',line[2],line[4],line[5],'.',strand,'.',line[2]]))


with open(out_file,'w') as out:
    for line in gff_entries:
        print >> out, line
