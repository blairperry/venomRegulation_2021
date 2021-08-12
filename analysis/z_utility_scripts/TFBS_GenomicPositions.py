
### Takes as input the result CSVs from CIIIDER - enrichment stats and TFBS sites - and 
### converts to a bedfile format w/ chromosome-level positions rather than positions within
### input regions. 

### Usage: python TFBS_GenomicPositions.py [Encirhment CSV] [TFBS Position CSV] [Output/dir/prefix]

import sys

enr_res = sys.argv[1]
pos_res = sys.argv[2]
outfile = sys.argv[3]

# enr_res = '/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation/analyses/_new_promoterATACpeaks/_venom_promoterPeaks/_tfbs_enrichment_binding/SVMP_vs_nonVenPAP/Enrichment: Text4699459774371289374_MostSigDeficit.csv'
# pos_res = '/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation/analyses/_new_promoterATACpeaks/_venom_promoterPeaks/_tfbs_enrichment_binding/SVMP_vs_nonVenPAP/Promoter Panel.csv'

# outfile = '/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation/analyses/_new_promoterATACpeaks/_venom_promoterPeaks/_tfbs_enrichment_binding/SVMP_vs_nonVenPAP/SVMP_TFBSGenCoord.bed'



enriched = []


if enr_res == 'all':
    with open(pos_res) as a, open(outfile,'w') as out:
        for line in a.readlines()[1:]:
            line = line.rstrip().split(',')
            # print line
            if ':' not in line[0]:
                print 'Ignoring ' + line[0]
            else:
                scaff = line[0].split(':')[2]
                region_start = int(line[0].split(':')[-1][:-3].split('-')[0])
                region_end = int(line[0].split(':')[-1][:-3].split('-')[1])
                

                tfbs_start = int(line[4])
                tfbs_stop = int(line[5])

                if line[0][-2] == '-':

                    new_tfbs_start = region_end - tfbs_stop
                    new_tfbs_stop = region_end - tfbs_start

                elif line[0][-2] == '+':

                    new_tfbs_start = region_start + tfbs_start
                    new_tfbs_stop = region_start + tfbs_stop

                else:
                    new_tfbs_start = region_start + tfbs_start
                    new_tfbs_stop = region_start + tfbs_stop


                print >> out, '\t'.join([scaff, str(new_tfbs_start), str(new_tfbs_stop), ';'.join(line)])
else:
    with open(enr_res) as a:
        for line in a.readlines()[1:]:
            line = line.rstrip().split(',')
            if float(line[-4]) < 0.05 and line[7] == 'Up':
                enriched.append(line[0])

    with open(pos_res) as a, open(outfile,'w') as out:
        for line in a.readlines():
            line = line.rstrip().split(',')
            if line[3] in enriched:
                
                scaff = line[0].split(':')[2]
                region_start = int(line[0].split(':')[-1][:-3].split('-')[0])
                region_end = int(line[0].split(':')[-1][:-3].split('-')[1])
                

                tfbs_start = int(line[4])
                tfbs_stop = int(line[5])


                if line[0][-2] == '-':

                	new_tfbs_start = region_end - tfbs_stop
                	new_tfbs_stop = region_end - tfbs_start

                else :

                	new_tfbs_start = region_start + tfbs_start
                	new_tfbs_stop = region_start + tfbs_stop

                print >> out, '\t'.join([scaff, str(new_tfbs_start), str(new_tfbs_stop), ';'.join(line)])