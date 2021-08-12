
# Given a multi-seq alignment, will calculate a simple consensus score at each position (i.e. percent of sequences with same base per position)

al_fasta = './example_files/aligned_input2.fa'
outfile = 'test_consensus.txt'

all_seqs = []

with open(al_fasta) as a:
    for line in a.readlines():
        if line[0] != '>':
            seq = []
            for pos in line.rstrip():
                seq.append(pos)
            all_seqs.append(seq)

scores = []
num_seqs = len(all_seqs)

for i in range(len(all_seqs[0])):
    a = 0
    t = 0
    c = 0
    g = 0
    gap = 0
    for seq in all_seqs:
        if seq[i] == 'A':
            a += 1
        elif seq[i] == 'T':
            t += 1
        elif seq[i] == 'C':
            c += 1
        elif seq[i] == 'G':
            g += 1
        elif seq[i] == '-':
            gap += 1
    scores.append(float(max([a,t,c,g]))/float(num_seqs))

with open(outfile,'w') as out:
    for i,score in enumerate(scores):
        print >> out, '\t'.join([str(i+1),str(score)])
