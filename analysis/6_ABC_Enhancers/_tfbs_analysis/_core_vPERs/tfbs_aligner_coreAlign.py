
# Script to convert TFBS positions in a raw fasta to their new positions an aligned/cropped fasta from Jalview

# If 'all' mode is used, also outputs all gaps per sequence and a consensus score for each position in the alignment
# Takes as input:
#   - type [tfbs,all] (defaults to tfbs)
#   - aligned fasta from Jalview
#   - Jalview GFF of tfbs positions*
#   - outfile directory

# Example usage:
# python tfbs_aligner_coreAlign.py all SVSP_CoreVperAlign_08.23.21.fa SVSP_CoreVperAlign_08.23.21.gff ./

import sys

type = sys.argv[1]
aligned_seqs = sys.argv[2]
raw_tfbs = sys.argv[3]

output_dir = sys.argv[4]
if output_dir[-1] != '/':
    output_dir = output_dir + '/'

out_prefix = output_dir + aligned_seqs.split('/')[-1][:-3]

do_tfbs = True
do_cons = False
do_gaps = False

if type == 'all':
    do_cons = True
    do_gaps = True

if (do_tfbs):

    tfbs_out = out_prefix + '_alignTFBS.txt'

    aligned_tfbs = []

    with open(aligned_seqs) as a:
        a = a.readlines()
        for j,line in enumerate(a):
            if line[0] != '>':
                seq_id = a[j-1].rstrip()[1:]
                crop_start = int(seq_id.split('/')[1].split('-')[0])
                # print crop_start
                crop_end = int(seq_id.split('/')[1].split('-')[1])
                seq_id_simple = seq_id.split('/')[0]
                coord_convert = {}
                old_pos = crop_start-1
                for i,char in enumerate(line.rstrip()):
                    new_pos = i+1
                    if char != '-':
                        old_pos += 1
                        coord_convert[old_pos] = new_pos
                # print coord_convert
                with open(raw_tfbs) as tfs:
                    for tfbs in tfs.readlines()[1:]:
                        tfbs = tfbs.rstrip().split()
                        # print tfbs
                        if tfbs[0] == seq_id_simple:
                            tf_id = tfbs[2]
                            # print tfbs[3]
                            old_start = int(tfbs[3])

                            if old_start == 0:
                                old_start = 1                       

                            old_stop = int(tfbs[4])

                            if old_start >= crop_start and old_stop <= crop_end:  
                                old_len = old_stop - old_start
                                # print tfbs
                                # print coord_convert
                                strand = tfbs[6]

                                # print old_start

                                new_start = coord_convert[old_start]
                                # print old_start, new_start
                                new_stop = coord_convert[old_stop]
                                new_len = new_stop - new_start
                                # print new_start, new_stop
                                if old_len != new_len:
                                    prev_coord = new_start-1;
                                    intervals = []
                                    seen = []
                                    for n in range(old_start,old_stop):
                                        if (coord_convert[n] - prev_coord) > 1:
                                            intervals.append(seen)
                                            seen = []
                                        prev_coord = coord_convert[n]
                                        seen.append(coord_convert[n])
                                        if n == old_stop-1:
                                            intervals.append(seen)
                                    for k,set in enumerate(intervals):
                                        set_id = tf_id + '_split' + str(k+1)
                                        set_start = min(set)
                                        set_stop = max(set)
                                        aligned_tfbs.append([seq_id, str(set_start), str(set_stop), set_id,strand,tfbs[2]])
                                elif old_len == new_len:
                                    aligned_tfbs.append([seq_id, str(new_start), str(new_stop), tf_id,strand,tfbs[2]])

    with open(tfbs_out,'w') as out:
        for entry in aligned_tfbs:
            print >> out, '\t'.join(entry)


if (do_cons):

    cons_out = out_prefix + '_consScore.txt'

    all_seqs = []

    with open(aligned_seqs) as a:
        for line in a.readlines():
            if line[0] != '>' and line[0] != '\n':
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
            if seq[i] == 'A' or seq[i] == 'a':
                a += 1
            elif seq[i] == 'T' or seq[i] == 't':
                t += 1
            elif seq[i] == 'C' or seq[i] == 'c':
                c += 1
            elif seq[i] == 'G' or seq[i] == 'g':
                g += 1
            elif seq[i] == '-':
                gap += 1
        scores.append(float(max([a,t,c,g]))/float(num_seqs))

    with open(cons_out,'w') as out:
        for i,score in enumerate(scores):
            print >> out, '\t'.join([str(i+1),str(score)])



if (do_gaps):

    gap_out = out_prefix + '_gaps.txt'

    all_gaps = []

    with open(aligned_seqs) as a:
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

    with open(gap_out,'w') as out:
        for line in all_gaps:
            print >> out, '\t'.join(line)
