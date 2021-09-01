
library(tidyverse)
library(ggbeeswarm)

all.bound <- read_tsv('analysis/6_ABC_Enhancers/_tfbs_analysis/_allBoundTFBS_vPERs_08.23.21.tsv')

all.tfbs <- read_tsv('analysis/3_transcription_factors/Motif_TFfamilies_08.17.21.txt',col_names = c('MA_ID','symbol','family')) %>% 
  mutate(symbol = str_split_fixed(symbol,'[(]',2)[,1]) %>% 
  select(symbol,MA_ID,family)

  
svmp.bound <- all.bound %>% 
  filter(vgfamily=='SVMP') %>% 
  select(vgfamily,tfbs_id) %>% 
  unique() %>% 
  left_join(all.tfbs,by=c('tfbs_id'='symbol')) %>% 
  select(MA_ID,tfbs_id,family)

svsp.bound <- all.bound %>% 
  filter(vgfamily=='SVSP') %>% 
  select(vgfamily,tfbs_id) %>% 
  unique() %>% 
  left_join(all.tfbs,by=c('tfbs_id'='symbol')) %>% 
  select(MA_ID,tfbs_id,family)

pla2.bound <- all.bound %>% 
  filter(vgfamily=='PLA2') %>% 
  select(vgfamily,tfbs_id) %>% 
  unique() %>% 
  left_join(all.tfbs,by=c('tfbs_id'='symbol')) %>% 
  select(MA_ID,tfbs_id,family)

# Write bound TFBS as target files for Jalview gff prep
# write_tsv(svmp.bound,'analysis/6_ABC_Enhancers/_tfbs_analysis/_core_vPERs/SVMP.boundTFs.tsv')
# write_tsv(svsp.bound,'analysis/6_ABC_Enhancers/_tfbs_analysis/_core_vPERs/SVSP.boundTFs.tsv')
# write_tsv(pla2.bound,'analysis/6_ABC_Enhancers/_tfbs_analysis/_core_vPERs/PLA2.boundTFs.tsv')


# Filter and plot CIIIDER scan of core vPER consensus seqs ----------------


# SVMP
svmp.consensusTFBS <- read_tsv('analysis/6_ABC_Enhancers/_tfbs_analysis/_core_vPERs/SVMP/SVMP_CoreVperAlign_08.23.21.singleLine_alignTFBS.txt',
                               col_names = c('PER_ID','start','stop','tfbs_splitID','strand','tfbs')) 
  
svmp.total_seqs = length(unique(svmp.consensusTFBS$PER_ID))


svmp.filteredTFBS <- svmp.consensusTFBS %>% 
  group_by(tfbs,start,stop,strand) %>% 
  tally() %>% 
  filter(n > svmp.total_seqs/2) %>% 
  mutate(roundedPosition = round((start+stop)/2,-1)) %>% 
  select(-strand) %>% 
  unique() %>% 
  left_join(all.tfbs,by=c('tfbs'='symbol')) %>% 
  mutate(vgfamily='SVMP',max=413)


ggplot(svmp.filteredTFBS,aes(x=roundedPosition,y='SVMP',color=family)) +
  # geom_quasirandom(groupOnX=FALSE,lwd=6,varwidth = TRUE) +
  geom_beeswarm(side=1L,lwd=4,cex=1) +
  theme_linedraw() 

ggplot(svmp.filteredTFBS,aes(x=roundedPosition,fill=family)) +
  geom_segment(aes(x=0,xend=max(svmp.filteredTFBS$roundedPosition+20),y=0,yend=0),color='grey70') +
  geom_dotplot(binwidth = 10,stackdir = 'center',stackgroups = T,binpositions = 'all') +
  xlab('Position') +
  theme_linedraw() + theme(panel.grid = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank())


# SVSP
svsp.consensusTFBS <- read_tsv('analysis/6_ABC_Enhancers/_tfbs_analysis/_core_vPERs/svsp/svsp_CoreVperAlign_08.23.21.singleLine_alignTFBS.txt',
                               col_names = c('PER_ID','start','stop','tfbs_splitID','strand','tfbs')) 

svsp.total_seqs = length(unique(svsp.consensusTFBS$PER_ID))


svsp.filteredTFBS <- svsp.consensusTFBS %>% 
  group_by(tfbs,start,stop,strand) %>% 
  tally() %>% 
  filter(n > svsp.total_seqs/2) %>% 
  mutate(roundedPosition = round((start+stop)/2,-1)) %>% 
  select(-strand) %>% 
  unique()%>% 
  left_join(all.tfbs,by=c('tfbs'='symbol')) %>% 
  mutate(vgfamily='SVSP',max=655)

ggplot(svsp.filteredTFBS,aes(x=roundedPosition,y='svsp',color=tfbs)) +
  # geom_quasirandom(groupOnX=FALSE,lwd=6,varwidth = TRUE) +
  geom_beeswarm(side=1L,lwd=4,cex=2) +
  theme_linedraw()


# PLA2
pla2.consensusTFBS <- read_tsv('analysis/6_ABC_Enhancers/_tfbs_analysis/_core_vPERs/pla2/pla2_CoreVperAlign_08.23.21.singleLine_alignTFBS.txt',
                               col_names = c('PER_ID','start','stop','tfbs_splitID','strand','tfbs')) 

pla2.total_seqs = length(unique(pla2.consensusTFBS$PER_ID))


pla2.filteredTFBS <- pla2.consensusTFBS %>% 
  group_by(tfbs,start,stop,strand) %>% 
  tally() %>% 
  filter(n > pla2.total_seqs/2) %>% 
  mutate(roundedPosition = round((start+stop)/2,-1)) %>% 
  select(-strand) %>% 
  unique()%>% 
  left_join(all.tfbs,by=c('tfbs'='symbol')) %>% 
  mutate(vgfamily='PLA2',max=148)

ggplot(pla2.filteredTFBS,aes(x=roundedPosition,y='pla2',color=tfbs)) +
  # geom_quasirandom(groupOnX=FALSE,lwd=6,varwidth = TRUE) +
  geom_beeswarm(side=1L,lwd=4,cex=2) +
  theme_minimal()


all.filteredTFBS <- svmp.filteredTFBS %>% 
  bind_rows(svsp.filteredTFBS,pla2.filteredTFBS) %>% 
  mutate(vgfamily=factor(vgfamily,levels=c('SVMP','SVSP','PLA2'))) %>% 
  arrange(family,tfbs) %>% 
  mutate(family = str_replace_all(family,'_',' ')) %>% 
  mutate(family = factor(family,levels=unique(.$family))) 

# ggplot(all.filteredTFBS,aes(x=roundedPosition,fill=tfbs)) +
#   geom_dotplot(binwidth = 10,stackdir = 'center',stackgroups = T,binpositions = 'all',) +
#   facet_wrap(~vgfamily,ncol = 1,strip.position = 'left') +
#   xlab('Position') +
#   theme_linedraw() + theme(panel.grid = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank())

ggplot(all.filteredTFBS,aes(x=start,xend=stop,y=tfbs,yend=tfbs,color=family)) +
  geom_segment(lwd=3) +
  facet_grid(rows=vars(vgfamily),space='free',scales = 'free') +
  theme_linedraw()


## Sequence alignments

library(Biostrings)

library(msa)

svmp.core <- readAAMultipleAlignment('./analysis/6_ABC_Enhancers/_tfbs_analysis/_core_vPERs/SVMP/SVMP_CoreVperAlign_08.23.21.singleLine.withConsensus.fa')
svsp.core <- readAAMultipleAlignment('./analysis/6_ABC_Enhancers/_tfbs_analysis/_core_vPERs/SVSP/SVSP_CoreVperAlign_08.23.21.singleLine.withConsensus.fa')
pla2.core <- readAAMultipleAlignment('./analysis/6_ABC_Enhancers/_tfbs_analysis/_core_vPERs/PLA2/PLA2_CoreVperAlign_08.23.21.singleLine.withConsensus.fa')



msaPrettyPrint(svmp.core,output = 'pdf',file = 'analysis/6_ABC_Enhancers/_tfbs_analysis/_core_vPERs/SVMP/SVMP_coreVPER_seqAlign_08.30.21.pdf',showLogo = 'none',showConsensus = 'none',askForOverwrite = F,paperWidth = 20,paperHeight = 20,showLegend = F)

msaPrettyPrint(svsp.core,output = 'pdf',file = 'analysis/6_ABC_Enhancers/_tfbs_analysis/_core_vPERs/SVSP/SVSP_coreVPER_seqAlign_08.30.21.pdf',showLogo = 'none',showConsensus = 'none',askForOverwrite = F,paperWidth = 20,paperHeight = 20,showLegend = F)

msaPrettyPrint(pla2.core,output = 'pdf',file = 'analysis/6_ABC_Enhancers/_tfbs_analysis/_core_vPERs/PLA2/PLA2_coreVPER_seqAlign_08.30.21.pdf',showLogo = 'none',showConsensus = 'none',askForOverwrite = F,paperWidth = 20,paperHeight = 20,showLegend = F)
