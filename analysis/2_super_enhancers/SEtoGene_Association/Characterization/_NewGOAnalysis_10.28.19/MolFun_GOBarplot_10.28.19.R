
library(tidyverse)

setwd('/Users/perryb/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/ChIPseq/SuperEnhancers/SEtoGene_Association/Characterization/_NewGOAnalysis_10.28.19/MolFuncResults/')

mfRes <- read_delim('enrichment_results_wg_result1572292198.txt',delim='\t') %>% 
  mutate(description=str_replace_all(description,', ',',\n'))

ggplot(mfRes,aes(x=reorder(description, enrichmentRatio),y=enrichmentRatio))+
  geom_bar(stat='identity',fill='seagreen') +
  xlab('')+
  ylab('Enrichment Ratio')+
  scale_y_continuous(position = 'right')+
  coord_flip() +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank())

