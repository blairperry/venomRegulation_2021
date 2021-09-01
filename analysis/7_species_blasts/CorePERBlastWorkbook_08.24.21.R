
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
