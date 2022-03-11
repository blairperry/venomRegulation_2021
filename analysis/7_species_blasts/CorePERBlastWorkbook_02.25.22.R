
library(tidyverse)
library(seqinr)
library(patchwork)

# Read and parse snake genom blast results --------------------------------

all.snakeBlast <- read_tsv('analysis/7_species_blasts/offline_snakeGenomeBlast/blast_results/allSpecies_coreVPER_BlastResults_08.24.21.txt',
                           col_names = c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'sseq','species')) %>% 
  mutate(species = str_split_fixed(species,'[_]',2)[,1])
all.snakeBlast

## Filter to include top 5 hits per species based on e-value first, and secondarily by bitscore for cases where more than 5 have identical lowest e-value
all.snakeBlast.top5 <- all.snakeBlast %>% 
  group_by(qseqid,species) %>% 
  top_n(n=5,wt = -evalue) %>% 
  top_n(n=5,wt = bitscore) %>% 
  filter(bitscore > 100)

all.snakeBlast.top5%>% group_by(qseqid) %>% tally()

descrip.ncbiBlast <- readxl::read_xlsx('analysis/7_species_blasts/ncbi_allSnakeBlast/CoreVPER_DescripTable_08.23.21.xlsx') %>% janitor::clean_names() %>% 
  mutate(scientific_name_temp1 = str_split_fixed(scientific_name,' ',3)[,1],
         scientific_name_temp2 = str_split_fixed(scientific_name,' ',3)[,2],
         scientific_name = paste(scientific_name_temp1,scientific_name_temp2,sep = '_')) %>% 
  select(accession,scientific_name,description)

all.snakeBlast.top5.svmp <- all.snakeBlast.top5 %>% filter(str_detect(qseqid,'SVMP')) %>% mutate(fasta_id = paste(sseqid,' ',sstart,'-',send,' ',species,sep = ''))
all.snakeBlast.top5.svsp <- all.snakeBlast.top5 %>% filter(str_detect(qseqid,'SVSP')) %>% mutate(fasta_id = paste(sseqid,' ',sstart,'-',send,' ',species,sep = ''))
all.snakeBlast.top5.pla2 <- all.snakeBlast.top5 %>% filter(str_detect(qseqid,'PLA2')) %>% mutate(fasta_id = paste(sseqid,' ',sstart,'-',send,' ',species,sep = ''))

# write.fasta(as.list(all.snakeBlast.top5.svmp$sseq),names=all.snakeBlast.top5.svmp$fasta_id,nbchar = 300,file.out = 'analysis/7_species_blasts/top5_combined/SVMP_CoreVPER_SnakeHits_08.24.21.fa')
# write.fasta(as.list(all.snakeBlast.top5.svsp$sseq),names=all.snakeBlast.top5.svsp$fasta_id,nbchar = 300,file.out = 'analysis/7_species_blasts/top5_combined/SVSP_CoreVPER_SnakeHits_08.24.21.fa')
# write.fasta(as.list(all.snakeBlast.top5.pla2$sseq),names=all.snakeBlast.top5.pla2$fasta_id,nbchar = 300,file.out = 'analysis/7_species_blasts/top5_combined/PLA2_CoreVPER_SnakeHits_08.24.21.fa')



# Read and parse NCBI blast results against Serpentes ---------------------

all.ncbiBlast <- read_csv('analysis/7_species_blasts/ncbi_allSnakeBlast/CoreVPER_FullResults_08.23.21.csv',col_names = c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')) %>% 
  left_join(descrip.ncbiBlast,by=c('sseqid'='accession'))


all.ncbiBlast.top5 <- all.ncbiBlast %>% 
  group_by(qseqid,scientific_name) %>% 
  top_n(n=5,wt = -evalue) %>% 
  top_n(n=5,wt = bitscore)%>% 
  filter(bitscore > 100) %>% 
  group_by(qseqid,scientific_name,sstart,send) %>% 
  mutate(fasta_id = paste(sseqid,':',min(c(sstart,send)),'-',max(c(sstart,send)),' ',description,sep = ''))

all.ncbiBlast.top5.svmp <- all.ncbiBlast.top5 %>% filter(str_detect(qseqid,'SVMP'))
all.ncbiBlast.top5.svsp <- all.ncbiBlast.top5 %>% filter(str_detect(qseqid,'SVSP'))
all.ncbiBlast.top5.pla2 <- all.ncbiBlast.top5 %>% filter(str_detect(qseqid,'PLA2'))


all.ncbiBlast.top5 %>% group_by(qseqid) %>% tally()

ncbi.hitSeqs <- read.fasta('analysis/7_species_blasts/ncbi_allSnakeBlast/all_CoreVPER_NCBI_HitSeqs_08.23.21.fa',seqtype = 'DNA',whole.header = T,as.string = T)

ncbi.hitSeqs.top5.svmp <- ncbi.hitSeqs[c(which(names(ncbi.hitSeqs) %in% all.ncbiBlast.top5.svmp$fasta_id))]
ncbi.hitSeqs.top5.svsp <- ncbi.hitSeqs[c(which(names(ncbi.hitSeqs) %in% all.ncbiBlast.top5.svsp$fasta_id))]
ncbi.hitSeqs.top5.pla2 <- ncbi.hitSeqs[c(which(names(ncbi.hitSeqs) %in% all.ncbiBlast.top5.pla2$fasta_id))]

# write.fasta(ncbi.hitSeqs.top5.svmp,names(ncbi.hitSeqs.top5.svmp),nbchar = 300,as.string = T,file.out = 'analysis/7_species_blasts/top5_combined/SVMP_CoreVPER_NCBI_08.24.21.fa')
# write.fasta(ncbi.hitSeqs.top5.svsp,names(ncbi.hitSeqs.top5.svsp),nbchar = 300,as.string = T,file.out = 'analysis/7_species_blasts/top5_combined/SVSP_CoreVPER_NCBI_08.24.21.fa')
# write.fasta(ncbi.hitSeqs.top5.pla2,names(ncbi.hitSeqs.top5.pla2),nbchar = 300,as.string = T,file.out = 'analysis/7_species_blasts/top5_combined/PLA2_CoreVPER_NCBI_08.24.21.fa')



# Filtering and plotting aligned TFBS ---------------------------------------------------

svmp.boundTFBS <- read_tsv('analysis/6_ABC_Enhancers/_tfbs_analysis/_core_vPERs/SVMP/SVMP.boundTFs.tsv')

svmp.aligned <- read_tsv('analysis/7_species_blasts/top5_combined/curated/ciiider_scans/SVMP_scans/SVMP_CoreVPER_AllCurated_08.24.21_alignTFBS.txt',
                         col_names = c('seq','start','stop','tfbs_splitID','ignore1','ignore2','ignore3','MA_ID')) %>% 
  mutate(tfbs_id = ifelse(str_detect(tfbs_splitID,'split'),str_split_fixed(tfbs_splitID,'[_]',2)[,1],tfbs_splitID)) %>% 
  filter(MA_ID %in% svmp.boundTFBS$MA_ID) %>% 
  mutate(genus = ifelse(str_detect(seq,'Crotalus'),'Crotalus',    
                          ifelse(str_detect(seq,'Deinag'),'Deinagkistrodon',
                                 ifelse(str_detect(seq,'Thamn'),'Thamnophis',
                                        ifelse(str_detect(seq,'Naja'),'Naja',
                                               ifelse(str_detect(seq,'Habu'),'Protobothrops',
                                                      ifelse(str_detect(seq,'vPER'),'Core vPER','Unknown'))))))) %>% 
  arrange(desc(genus)) %>% 
  mutate(seq = factor(seq,levels=unique(.$seq)))


svmp.cons_score <- read_tsv('analysis/7_species_blasts/top5_combined/curated/ciiider_scans/SVMP_scans/SVMP_CoreVPER_AllCurated_08.24.21_consScore.txt',col_names = c('position','score'))

ggplot(svmp.aligned,aes(x=start,xend=stop,y=seq,yend=seq,color=tfbs_id)) +
  geom_segment(alpha=0.65,lwd=3) +
  ylab('')+
  xlab('Relative Position') +
  coord_cartesian(xlim =c(0,max(svmp.aligned$stop + 15)),clip = 'on',expand = T)+
  theme_classic() +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        legend.position = 'left',
        strip.text.y.left = element_text(angle = 0,color='black',face = 'bold'),
        axis.title.x = element_blank(),
        axis.text.y=element_text(size=8),
        strip.background = element_rect(fill='grey90',size = 0.5),
        panel.spacing = unit(0.25, "lines"))

p1 <- ggplot(svmp.aligned,aes(x=start,xend=stop,y=seq,yend=seq,color=tfbs_id)) +
  # geom_segment(inherit.aes = F,data=svmp.splits,aes(x=split_start,xend=split_end,y=0.5*X5*X6,yend=0.5*X5*X6,color=cluster_col),lwd=1,alpha=0.75) +
  geom_segment(alpha=0.65,lwd=2) +
  # geom_point(inherit.aes = F,aes(x=-10,y=seq,fill=genus),pch=21,size=2) +
  scale_fill_manual(values=c('viridis'='red','atrox'='goldenrod','horridus'='darkslategray1','scutulatus'='green4','tigris'='black','durissus'='grey')) +
  # geom_segment(inherit.aes = F,aes(y=0,yend=0,x=0,xend=max(svmp.ms.cons_score$X1+20)),lwd=0.75,color='grey60') +
  # geom_segment(data=svmp.ms.tf.ali,inherit.aes = F,aes(x=X2,xend=X3,y=0,yend=0),color='white',lwd=2,alpha=0.85) +
  ylab('')+
  xlab('Relative Position') +
  coord_cartesian(xlim =c(0,max(svmp.aligned$stop + 15)),clip = 'on',expand = T)+
  # scale_y_continuous(position = 'right',limits = c(-abs(max(svmp.ms.tf.ali$X6)*1.2),abs(max(svmp.ms.tf.ali$X6)*1.2))) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        legend.position = 'left',
        strip.text.y.left = element_text(angle = 0,color='black',face = 'bold'),
        axis.title.x = element_blank(),
        axis.text.y=element_text(size=8),
        strip.background = element_rect(fill='grey90',size = 0.5),
        panel.spacing = unit(0.25, "lines"))
p1

p2 <- ggplot(svmp.cons_score,aes(x=position,y=score)) +
  geom_area(fill='dodgerblue3',alpha=0.5) +
  scale_y_continuous(limits=c(0,1.1))+
  coord_cartesian(xlim =c(0,max(svmp.aligned$stop + 15)),clip = 'on',expand = T)+
  # scale_fill_continuous_sequential('Mint') +
  ylab('Consensus Score') +
  xlab('Relative Position')+
  theme_linedraw() + theme(panel.grid.major.y = element_line(),panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),legend.position = 'left')
 
p1 + p2 + plot_layout(heights = c(10,1))


# Revision work: highlighting viper-specific, elapid-specific, shared TFBS --------

unique(svmp.aligned$genus)

num.elapid

svmp.aligned.groups <- svmp.aligned %>% 
  mutate(group = case_when(
    genus %in% c('Protobothrops','Crotalus','Deinagkistrodon') ~ 'vipers',
    genus == 'Naja' ~ 'elapid',
    genus %in% c('Thamnophis','Core vPER') ~ 'other')) 

numSeqs <- svmp.aligned.groups %>% 
  select(seq,group) %>% 
  unique() %>% 
  group_by(group) %>% 
  tally()

numSeqs

num.elapid <- 5
num.viper <- 28

svmp.aligned.groups.spec <- svmp.aligned.groups %>% 
  filter(group != 'other') %>% 
  group_by(tfbs_splitID,start,stop,group) %>% 
  select(-contains('ignore')) %>% 
  unique() %>% 
  add_count() %>% 
  mutate(is_majority = case_when(
    group == 'vipers' ~ ifelse(n > (num.viper/2),T,F),
    group == 'elapid' ~ ifelse(n > (num.elapid/2),T,F)
  )) %>% 
  ungroup() %>% 
  select(tfbs_id,start,stop,group,is_majority) %>% 
  unique() %>% 
  pivot_wider(names_from=group,values_from = is_majority,values_fill = F) %>% 
  filter(vipers == T | elapid == T) %>% 
  mutate(specific_group = case_when(
    vipers == T & elapid == T ~ 'shared',
    vipers == T & elapid == F ~ 'viper-specific',
    vipers == F & elapid == T ~ 'elapid-specific'
  )) %>% 
  mutate(specific_group = factor(specific_group,levels=c('shared','viper-specific','elapid-specific')))
  
svmp.aligned.groups.spec %>% 
  group_by(specific_group) %>% 
  tally() %>% 
  ggplot(aes(x=specific_group,y=n)) +
  geom_bar(stat='identity') +
  theme_classic()
  
p3 <- ggplot(svmp.aligned.groups.spec,aes(x=start,xend=stop,y=reorder(specific_group,desc(specific_group)),yend=reorder(specific_group,desc(specific_group)),color=specific_group)) +
  geom_segment(lwd=5,show.legend = F) +
  coord_cartesian(xlim =c(0,max(svmp.aligned$stop + 15)),clip = 'on',expand = T)+
  scale_color_manual(values = c('elapid-specific'='seagreen','viper-specific'='dodgerblue3','shared'='goldenrod')) +
  theme_linedraw()

p3 + p1 + p2 + plot_layout(heights = c(2,10,1))


svmp.aligned.groups.spec.totals <- svmp.aligned.groups.spec %>% 
  group_by(tfbs_id,specific_group) %>% 
  tally() %>% 
  group_by(tfbs_id) %>% 
  mutate(total = sum(n)) %>% 
  arrange(total) %>% 
  # mutate(tfbs_id = factor(tfbs_id,levels=unique(.$tfbs_id))) %>% 
  mutate(specific_group = factor(specific_group,levels=c('shared','viper-specific','elapid-specific')))


spec.p1 <- ggplot(svmp.aligned.groups.spec.totals,aes(y=reorder(tfbs_id,desc(tfbs_id)),x=n,fill=specific_group)) +
  geom_bar(stat='identity',orientation = 'y') +
  scale_fill_manual(values = c('elapid-specific'='seagreen','viper-specific'='dodgerblue3','shared'='black')) +
  scale_x_continuous(breaks=seq(1,13,1)) +
  labs(x= 'Number of TFBS Positions',y='TFBS',fill='Group') +
  theme_linedraw() + theme(legend.position = c(0.75,0.15),legend.background = element_rect(color='black'),panel.grid.minor = element_blank())

spec.p1.v2 <- ggplot(svmp.aligned.groups.spec.totals,aes(y=reorder(tfbs_id,desc(tfbs_id)),x=specific_group,fill=specific_group)) +
  geom_point(aes(size=n),pch=23,show.legend = T) +
  scale_size_continuous(range=c(3,7)) +
  scale_fill_manual(values = c('elapid-specific'='seagreen','viper-specific'='dodgerblue3','shared'='goldenrod')) +
  guides(fill='none') +
  labs(x= '',y='TFBS',size='# TFBS Sites') +
  theme_linedraw() + theme(legend.position = 'bottom',legend.background = element_rect(color='black'),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))

tf.characterization <- read_tsv('analysis/3_transcription_factors/_functional_characterization/CandTFs_InterestingSubset_FULL_08.27.21.tsv') %>% 
  mutate(value = ifelse(id == 'EHF' & feature == 'Interacts with ERK',T,value)) %>% 
  filter(value == T) 

spec.p2 <- svmp.aligned.groups.spec %>% 
  left_join(tf.characterization,by=c('tfbs_id'='id')) %>% 
  # mutate(tfbs_id = factor(tfbs_id,levels=unique(svmp.aligned.groups.spec.totals$tfbs_id))) %>% 
  ggplot(aes(x=reorder(feature,desc(feature)),y=reorder(tfbs_id,desc(tfbs_id)),fill=feature)) +
  geom_point(size=4,show.legend = F,pch=21) +
  # scale_size_continuous(range = c(2,8)) +
  scale_size_area(max_size = 6) +
  # scale_fill_manual(values = c('PLA2'='#A23B5F','SVSP'='#7DCBAC','SVMP'='#A3C0CF','Other'='grey60'))+
  scico::scale_fill_scico_d(direction = -1) +
  xlab('') +
  ylab('') +
  # ggtitle('Venom Gene Promoters - Putatively-Bound TFBS') +
  theme_linedraw() + theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = 'bottom',
                           axis.text.y=element_blank(),
                           plot.title.position = 'plot',plot.title = element_text(face='bold'),strip.text.x = element_blank(),
                           strip.text.y = element_text(angle = 0,hjust = 0,colour = 'black'),strip.background = element_rect(fill='grey95',color='NA') )

# spec.p1 + spec.p2 + plot_layout(widths = c(3,2))
spec.p1.v2 + spec.p2 + plot_layout(widths = c(3,2))

#


## Pla2

pla2.boundTFBS <- read_tsv('analysis/6_ABC_Enhancers/_tfbs_analysis/_core_vPERs/pla2/pla2.boundTFs.tsv')

pla2.aligned <- read_tsv('analysis/7_species_blasts/top5_combined/curated/ciiider_scans/pla2_scans/pla2_CoreVPER_AllCurated_08.24.21_alignTFBS.txt',
                         col_names = c('seq','start','stop','tfbs_splitID','ignore1','ignore2','ignore3','MA_ID')) %>% 
  mutate(tfbs_id = ifelse(str_detect(tfbs_splitID,'split'),str_split_fixed(tfbs_splitID,'[_]',2)[,1],tfbs_splitID)) %>% 
  filter(MA_ID %in% pla2.boundTFBS$MA_ID) %>% 
  mutate(genus = ifelse(str_detect(seq,'Crotalus'),'Crotalus',    
                        ifelse(str_detect(seq,'Deinag'),'Deinagkistrodon',
                               ifelse(str_detect(seq,'Thamn'),'Thamnophis',
                                      ifelse(str_detect(seq,'Ovophis'),'Ovophis',
                                             ifelse(str_detect(seq,'Trimeresurus'),'Trimeresurus',
                                                    ifelse(str_detect(seq,'vPER'),'Core vPER','Unknown'))))))) %>% 
  arrange(desc(genus)) %>% 
  mutate(seq = ifelse(str_detect(seq,'[;]'),str_split_fixed(seq,'[;]',2)[,1],seq)) %>% 
  mutate(seq = factor(seq,levels=unique(.$seq)))


pla2.cons_score <- read_tsv('analysis/7_species_blasts/top5_combined/curated/ciiider_scans/pla2_scans/pla2_CoreVPER_AllCurated_08.24.21_consScore.txt',col_names = c('position','score'))

ggplot(pla2.aligned,aes(x=start,xend=stop,y=seq,yend=seq,color=tfbs_id)) +
  geom_segment(alpha=0.65,lwd=3) +
  ylab('')+
  xlab('Relative Position') +
  coord_cartesian(xlim =c(0,max(pla2.aligned$stop + 15)),clip = 'on',expand = T)+
  theme_classic() +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        legend.position = 'left',
        strip.text.y.left = element_text(angle = 0,color='black',face = 'bold'),
        axis.title.x = element_blank(),
        axis.text.y=element_text(size=8),
        strip.background = element_rect(fill='grey90',size = 0.5),
        panel.spacing = unit(0.25, "lines"))

pla2.p1 <- ggplot(pla2.aligned,aes(x=start,xend=stop,y=seq,yend=seq,color=tfbs_id)) +
  # geom_segment(inherit.aes = F,data=pla2.splits,aes(x=split_start,xend=split_end,y=0.5*X5*X6,yend=0.5*X5*X6,color=cluster_col),lwd=1,alpha=0.75) +
  geom_segment(alpha=0.65,lwd=2) +
  # geom_point(inherit.aes = F,aes(x=-10,y=seq,fill=genus),pch=21,size=2) +
  scale_fill_manual(values=c('viridis'='red','atrox'='goldenrod','horridus'='darkslategray1','scutulatus'='green4','tigris'='black','durissus'='grey')) +
  # geom_segment(inherit.aes = F,aes(y=0,yend=0,x=0,xend=max(pla2.ms.cons_score$X1+20)),lwd=0.75,color='grey60') +
  # geom_segment(data=pla2.ms.tf.ali,inherit.aes = F,aes(x=X2,xend=X3,y=0,yend=0),color='white',lwd=2,alpha=0.85) +
  ylab('')+
  xlab('Relative Position') +
  coord_cartesian(xlim =c(0,max(pla2.aligned$stop+15)),clip = 'on',expand = T)+
  # scale_y_continuous(position = 'right',limits = c(-abs(max(pla2.ms.tf.ali$X6)*1.2),abs(max(pla2.ms.tf.ali$X6)*1.2))) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        legend.position = 'left',
        strip.text.y.left = element_text(angle = 0,color='black',face = 'bold'),
        axis.title.x = element_blank(),
        axis.text.y=element_text(size=8),
        strip.background = element_rect(fill='grey90',size = 0.5),
        panel.spacing = unit(0.25, "lines"))
pla2.p1

pla2.p2 <- ggplot(pla2.cons_score,aes(x=position,y=score)) +
  geom_area(fill='dodgerblue3',alpha=0.5) +
  scale_y_continuous(limits=c(0,1.1))+
  coord_cartesian(xlim =c(0,max(pla2.aligned$stop + 15)),clip = 'on',expand = T)+
  # scale_fill_continuous_sequential('Mint') +
  ylab('Consensus Score') +
  xlab('Relative Position')+
  theme_linedraw() + theme(panel.grid.major.y = element_line(),panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),legend.position = 'left')

pla2.p1 + pla2.p2 + plot_layout(heights = c(6,1))


p1 + p2 + pla2.p1 + pla2.p2 + plot_layout(heights = c(6,1,4,1))
