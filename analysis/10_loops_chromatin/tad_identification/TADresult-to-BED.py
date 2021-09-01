
import os

# USAGE:
# Run this in a directory of rGMAP output files and it will automatically
# convert all TADs to BED format, with separate files for different hierachical
# classifications of TADs

# Assumes that files are named like this "ma1_10kb_HierTADS.txt"

curr_dir = os.getcwd()

for filename in os.listdir(curr_dir):
    if filename[-4:] == '.txt':
        openname = curr_dir + '/' + filename
        with open(openname) as a:
            outname1 = curr_dir + "/BED/" + filename[:-4] + '_hier1.bed'
            outname2 = curr_dir + "/BED/" + filename[:-4] + '_hier2.bed'
            h1 = []
            h2 = []
            scaffold_id = 'scaffold-' + filename.split('_')[0]
            for line in a.readlines()[1:]:
                line = line.rstrip().split('\t')
                order = line[2]
                newline = [scaffold_id,line[0],line[1]]
                if line[2] == '1':
                    h1.append(newline)
                elif line[2] == '2':
                    h2.append(newline)
            with open(outname1,'w') as o:
                for i, line in enumerate(h1):
                    tadID = 'TAD_' + filename[:3] + '_h1_' + str(i+1)
                    outline = [line[0],line[1],line[2],tadID]
                    print >> o, '\t'.join(outline)
            if len(h2) > 0:
                with open(outname2,'w') as o:
                    for i, line in enumerate(h2):
                        tadID = 'TAD_' + filename[:3] + '_h2_' + str(i+1)
                        outline = [line[0],line[1],line[2],tadID]
                        print >> o, '\t'.join(outline)
