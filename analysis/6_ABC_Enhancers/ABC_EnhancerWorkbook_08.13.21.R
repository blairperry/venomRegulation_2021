
library(tidyverse)
library(colorspace)
library(patchwork)

venomTxIDs <- read_tsv('./data/venom_annotation/PriorityVenomGenes_08.02.21.txt',col_names = c(
  'chr','start','stop','long_name','messy_name','tx_id','short_name'
)) %>% 
  filter(str_detect(short_name, 'gII|ADAM',negate = T))

all.predics <- read_delim('./analysis/6_ABC_Enhancers/ABC_output/Predictions/EnhancerPredictionsFull.txt',delim='\t') %>% 
  mutate(geneType = ifelse(TargetGene %in% venomTxIDs$tx_id,'Venom','NonVenom'))

# highestExp <- read_tsv('../gene_expression/HighExpressed_VenGenes_10.14.20.tsv')

# Total number of unique enhancer sequences
all.predics %>% 
  select(name) %>% 
  unique() %>% 
  tally()

# Total number of unique genes
all.predics %>% 
  select(TargetGene) %>% 
  unique() %>% 
  tally()

# % of genic and intergenic PERs
all.predics %>% 
  select(name,class) %>% 
  unique() %>% 
  group_by(class) %>% 
  tally() %>% 
  mutate(n = n / sum(n))


# Median distance from gene
all.predics %>% 
  select(name,TargetGene,distance) %>% 
  unique() %>% 
  summarise(med.dist = median(distance))


geneStats <- read_delim('./analysis/6_ABC_Enhancers/ABC_output/Predictions/GenePredictionStats.txt',delim='\t') %>% 
  mutate(geneType = ifelse(TargetGene %in% venomTxIDs$tx_id,'Venom','NonVenom'))

geneStats %>% 
  group_by(geneIsExpressed) %>% 
  tally()

geneStats %>% 
  filter(geneIsExpressed == TRUE) %>% 
  mutate(PERs = ifelse(nDistalEnhancersPredicted > 0,1,0)) %>% 
  group_by(PERs) %>% 
  tally()

t <- geneStats %>% 
  filter(geneIsExpressed == TRUE) %>% 
  filter(nDistalEnhancersPredicted > 0)

mean(t$nDistalEnhancersPredicted)


venom.predics <- all.predics %>% 
  filter(geneType == 'Venom') %>% 
  left_join(venomTxIDs,by=c('TargetGene'='tx_id'))

### Number of Enhancers per Gene plots

geneStats %>% 
  group_by(chr,nDistalEnhancersPredicted) %>% 
  tally() %>% 
  ggplot(aes(x=nDistalEnhancersPredicted,y=n)) +
  geom_bar(stat='identity',fill='dodgerblue4') +
  scale_y_continuous(expand = c(0,0)) +
  ylab('Frequency') +
  xlab('# of Predicted Enhancers') +
  ggtitle('Number of Enhancers Per Gene') +
  theme_classic()


geneStats %>% 
  group_by(geneType,nDistalEnhancersPredicted) %>% 
  tally() %>% 
  mutate(proportion = n / sum(n)) %>% 
  ggplot(aes(x=nDistalEnhancersPredicted,y=proportion)) +
  geom_bar(stat='identity') +
  facet_wrap(~geneType)+
  ylab('Proportion') +
  xlab('# of Predicted Enhancers') +
  theme_classic()

# Pie chart
p.pie <- geneStats %>% 
  mutate(simple_nEnhancers = ifelse(nDistalEnhancersPredicted >= 5, '>= 5',as.character(nDistalEnhancersPredicted))) %>% 
  mutate(simple_nEnhancers = factor(simple_nEnhancers,levels=c('0','1','2','3','4','>= 5'))) %>% 
  group_by(geneType,simple_nEnhancers) %>% 
  tally() %>% 
  mutate(proportion = n / sum(n)) %>% 
  ggplot(aes(x='',y=proportion,fill=simple_nEnhancers)) +
  geom_bar(stat='identity',width = 1,color='white') +
  scale_fill_discrete_sequential('lajolla') +
  coord_polar('y',start = 0,direction = -1)+
  facet_wrap(~geneType) +
  labs(fill='Number of PERs per gene') +
  theme_void(base_size = 16) + theme(legend.position = 'bottom')
  
p.pie

### Number of genes per enhancer plots

all.predics %>% 
  group_by(name,class) %>% 
  tally() %>% 
  ggplot(aes(x=n)) +
  geom_histogram(fill='dodgerblue4',bins=40) +
  ylab('Frequency') +
  xlab('# of Genes per Enhancer') +
  scale_x_continuous(limits=c(0,20)) +
  theme_classic()


p.dens <- ggplot(all.predics,aes(x=distance)) +
  geom_density(fill='grey',alpha=0.38) +
  geom_vline(xintercept = median(all.predics$distance),lty=2) +
  geom_text(aes(label=median(all.predics$distance),x=median(all.predics$distance),y=0.5),nudge_x = .2,check_overlap = T) +
  geom_vline(xintercept = median(venom.predics$distance),lty=2,color='red') +
  geom_text(aes(label=median(venom.predics$distance),x=median(venom.predics$distance),y=0.55),nudge_x = -.2,check_overlap = T,color='red') +
  scale_x_log10(expand=c(0,0),labels = scales::comma) +
  xlab('Distance between PER and gene (bp)') +
  # scale_x_continuous(limits=c(0,50000)) +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.7)) +
  theme_classic(base_size = 12)

ggplot(all.predics,aes(x=distance,y=geneType,fill=geneType)) +
  geom_boxplot(alpha=0.38) +
  geom_jitter() +
  # geom_vline(xintercept = median(all.predics$distance),lty=2) +
  # geom_text(aes(label=median(all.predics$distance),x=median(all.predics$distance),y=0.5),nudge_x = .2,check_overlap = T) +
  # geom_vline(xintercept = median(venom.predics$distance),lty=2,color='red') +
  # geom_text(aes(label=median(venom.predics$distance),x=median(venom.predics$distance),y=0.5),nudge_x = -.2,check_overlap = T,color='red') +
  scale_x_log10(expand=c(0,0)) +
  # scale_x_continuous(limits=c(0,50000)) +
  # scale_y_continuous(expand=c(0,0),limits=c(0,0.6),labels = scales::comma) +
  theme_classic()

###

p.pie / p.dens

##



ggplot(all.predics,aes(x=distance)) +
  geom_density(fill='dodgerblue',alpha=0.38) +
  geom_vline(xintercept = median(all.predics$distance),lty=2) +
  geom_vline(xintercept = median(venom.predics$distance),lty=2) +
  geom_text(aes(label=median(all.predics$distance),x=median(all.predics$distance),y=0.5),check_overlap = T) +
  scale_x_log10(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.7)) +
  theme_classic()



# % of genic and intergenic PERs
venom.predics %>% 
  select(name,class) %>% 
  unique() %>% 
  group_by(class) %>% 
  tally() %>% 
  mutate(n = n / sum(n))

# Median distance from gene
venom.predics %>% 
  select(name,TargetGene,distance) %>% 
  unique() %>% 
  summarise(med.dist = median(distance))



###


venom.predics.simple <- venom.predics %>% 
  select(chr.x,start.x,end,TargetGene,short_name) %>% 
  mutate(newID1 = ifelse(str_detect(TargetGene,'SVMP|SVSP|PLA2'),str_remove(TargetGene,'crovir-transcript-'),paste(str_replace(TargetGene,'crovir-transcript-','cvt'),short_name,sep='-'))) %>% 
  group_by(chr.x,start.x,end) %>% 
  summarise(newID2 = paste(newID1,collapse ='.')) %>% 
  arrange(chr.x,start.x) %>% 
  rowid_to_column() %>% 
  mutate(newID3 = paste('PER',rowid,'_',newID2,sep = '')) %>% 
  select(chr.x,start.x,end,newID3)

write_tsv(venom.predics.simple,'./analysis/6_ABC_Enhancers/ABC_output/_reformat/EnhancerPredictionsFull_VenomGenes_simple_newID_08.14.21.bed',col_names = F)

write_tsv(venom.predics.simple %>% filter(str_detect(newID3,'SVMP')),'./analysis/6_ABC_Enhancers/ABC_output/_reformat/SVMP_EnhancerPredictionsFull_VenomGenes_simple_newID_08.14.21.bed',col_names = F)
write_tsv(venom.predics.simple %>% filter(str_detect(newID3,'SVSP')),'./analysis/6_ABC_Enhancers/ABC_output/_reformat/SVSP_EnhancerPredictionsFull_VenomGenes_simple_newID_08.14.21.bed',col_names = F)
write_tsv(venom.predics.simple %>% filter(str_detect(newID3,'PLA2')),'./analysis/6_ABC_Enhancers/ABC_output/_reformat/PLA2_EnhancerPredictionsFull_VenomGenes_simple_newID_08.14.21.bed',col_names = F)
write_tsv(venom.predics.simple %>% filter(str_detect(newID3,'cvt')),'./analysis/6_ABC_Enhancers/ABC_output/_reformat/Other_EnhancerPredictionsFull_VenomGenes_simple_newID_08.14.21.bed',col_names = F)

### Haven't updated this highest expressed chunk
# highexp.venom.predics.simple <- venom.predics %>% 
#   select(chr.x,start.x,end,TargetGene,short_name) %>% 
#   filter(short_name %in% highestExp$gene) %>% 
#   left_join(venom.predics.simple) %>% 
#   select(chr.x,start.x,end,newID3)
# 
# # write_tsv(highexp.venom.predics.simple,'output/ABC_Output_Dec2020/Predictions/EnhancerPredictionsFull_HighestExpVenomGenes_simple_newID.bed',col_names = F)
# # 
# # write_tsv(highexp.venom.predics.simple %>% filter(str_detect(newID3,'SVMP')),'output/ABC_Output_Dec2020/Predictions/highSVMP_EnhancerPredictionsFull_VenomGenes_simple_newID.bed',col_names = F)
# # write_tsv(highexp.venom.predics.simple %>% filter(str_detect(newID3,'SVSP')),'output/ABC_Output_Dec2020/Predictions/highSVSP_EnhancerPredictionsFull_VenomGenes_simple_newID.bed',col_names = F)
# # write_tsv(highexp.venom.predics.simple %>% filter(str_detect(newID3,'PLA2')),'output/ABC_Output_Dec2020/Predictions/highPLA2_EnhancerPredictionsFull_VenomGenes_simple_newID.bed',col_names = F)
# # write_tsv(highexp.venom.predics.simple %>% filter(str_detect(newID3,'cvt')),'output/ABC_Output_Dec2020/Predictions/highOther_EnhancerPredictionsFull_VenomGenes_simple_newID.bed',col_names = F)
# 

nonVen.predicts <- all.predics %>% 
  filter(geneType != 'Venom')

nonVen.predicts.simple <- nonVen.predicts %>% 
  select(1,2,3,4) %>% 
  unique()

write_tsv(nonVen.predicts.simple,'./analysis/6_ABC_Enhancers/ABC_output/_reformat/EnhancerPredictionsFull_NONvenomGenes_simple_noDups_08.14.21.bed',col_names = F)


## Full detail table for supplement


fulldetail <- all.predics %>% 
  left_join(venom.predics.simple,by=c('chr'='chr.x','start'='start.x','end'='end'))

write_tsv(fulldetail,'./analysis/6_ABC_Enhancers/ABC_output/_reformat/EnhancerPredictionsFull_wVPERIDs_08.14.21.tsv')
