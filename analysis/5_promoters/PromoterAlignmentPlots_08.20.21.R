
library(ggmsa)
library(Biostrings)

svmp.align <- Biostrings::readDNAMultipleAlignment('analysis/5_promoters/_tfbs_analysis/aligned/SVMP_PromPeaks_seqs_08.16.21.aln.mafft.singleLine.fa')
svsp.align <- Biostrings::readDNAMultipleAlignment('analysis/5_promoters/_tfbs_analysis/aligned/SVSP_PromPeaks_seqs_08.16.21.aln.mafft.singleLine.fa')
pla2.align <- Biostrings::readDNAMultipleAlignment('analysis/5_promoters/_tfbs_analysis/aligned/PLA2_PromPeaks_seqs_08.16.21.aln.mafft.singleLine.fa')


ggmsa(svmp.align,170,300,color = 'Chemistry_NT',font = NULL, border = NA,by_conservation = T) + geom_msaBar()

ggmsa(pla2.align,40,190,color = 'Chemistry_NT',font = NULL, border = NA,by_conservation = T)

