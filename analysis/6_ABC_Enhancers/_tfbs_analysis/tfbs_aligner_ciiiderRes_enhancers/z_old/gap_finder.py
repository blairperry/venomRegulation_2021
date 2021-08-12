
# Outputs positions of all gap regions for each sequence in an aligned fasta

al_fasta = './example_files/aligned_input2.fa'
outfile = 'test_gaps.txt'

all_gaps = []

with open(al_fasta) as a:
    a = a.readlines()
    for i,line in enumerate(a):
        if line[0] == '>':
            seq_id = line[1:].rstrip()
            last_gap_pos = 0;
            gap = []
            for j,char in enumerate(a[i+1].rstrip()):
                if char == '-':
                    if ((j+1)-last_gap_pos) == 1:
                        gap.append(j)
                        last_gap_pos = j+1
                    elif ((j+1)-last_gap_pos) > 1:
                        gap.append(last_gap_pos)
                        all_gaps.append([seq_id,str(min(gap)),str(max(gap))])
                        gap = []
                        last_gap_pos = j+1
            if len(gap) > 0:
                gap.append(last_gap_pos)
                all_gaps.append([seq_id,str(min(gap)),str(max(gap))])

with open(outfile,'w') as out:
    for line in all_gaps:
        print >> out, '\t'.join(line)
