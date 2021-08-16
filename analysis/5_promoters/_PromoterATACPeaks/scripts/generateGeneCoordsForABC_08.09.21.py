
##
## This script generates the gene coordinate file needed for input in ABC enhancer inference analyses.
## Normally, this file would have gene start-stop positions from a normal GFF, but since we have revised
## the likely TSS regions for some genes based on ATAC-seq peaks falling within 1kb of annotated TSS, we
## need to generate a new gene coord file for genes w/ ATAC peak has midpoint of ATAC peak through gene end,
## and otherwise just keeps the normal gene start-stop coordinates if an ATAC peak isn't present
##

# Desired output format:
# scaffold-ma3	127727929	127748472	crovir-transcript-5997	0	+


## Parse ATAC-peak regions first and save in dictionary

infile_atacPeaks = '/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation_NEW_Aug2021/analysis/5_promoters/_PromoterATACPeaks/venom_promoterPeaks/allThree_ATACpeaks_PromWindows_08.16.21.bed'

atac_peaks = {}

with open(infile_atacPeaks) as a1:
    for line in a1.readlines():
        line = line.rstrip().split('\t')
        if 'SVMP' in line[3] or 'SVSP' in line[3] or 'PLA2' in line[3]:
            midpoint = (int(line[1]) + int(line[2])) / 2
            atac_peaks[line[3]] = midpoint



## Parse GFF gene entries

infile_gff = '/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation_NEW_Aug2021/data/annotation/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019.GeneEntriesOnly.gff'

outfile = '/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation_NEW_Aug2021/analysis/5_promoters/CvvGeneCoords_forABC_08.16.21.bed'

with open(infile_gff) as a2, open(outfile,'w') as out:
    for line in a2.readlines():
        line = line.rstrip().split('\t')
        if 'un' not in line[0] and 'trnascan' not in line[8]:

            id = line[8].split(';')[2][21:]

            if id in atac_peaks:
                if line[6] == '+':
                    outlist =  [ line[0], str(atac_peaks[id]), line[4], id, '0', line[6] ]
                elif line[6] == '-':
                    outlist =  [ line[0], line[3], str(atac_peaks[id]), id, '0', line[6] ]
            else:
                outlist = [ line[0], line[3], line[4], id, '0', line[6] ]
            print >> out, '\t'.join(outlist)
