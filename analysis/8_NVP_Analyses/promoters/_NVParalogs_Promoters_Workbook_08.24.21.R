
library(tidyverse)
library(patchwork)

svmp.nvp.promRes <- read_csv('analysis/8_NVP_Analyses/promoters/svmpNVP_prom_ciiiderRes/Enrichment: Text8805771021201769772_MostSigDeficit.csv') %>% 
  janitor::clean_names() %>% 
  mutate(nvp_family = 'SVMP')

svsp.nvp.promRes <- read_csv('analysis/8_NVP_Analyses/promoters/svspNVP_prom_ciiiderRes/Enrichment: Text6595001954856759384_MostSigDeficit.csv') %>% 
  janitor::clean_names()%>% 
  mutate(nvp_family = 'SVSP')

pla2.nvp.promRes <- read_csv('analysis/8_NVP_Analyses/promoters/pla2NVP_prom_ciiiderRes/Enrichment: Text3232898237227925445_MostSigDeficit.csv') %>% 
  janitor::clean_names()%>% 
  mutate(nvp_family = 'PLA2')

svmp.nvp.sigUp <- svmp.nvp.promRes %>% 
  filter(gene_p_value < 0.05 & gene_representation == 'Up')

svsp.nvp.sigUp <- svsp.nvp.promRes %>% 
  filter(gene_p_value < 0.05 & gene_representation == 'Up')

pla2.nvp.sigUp <- pla2.nvp.promRes %>% 
  filter(gene_p_value < 0.05 & gene_representation == 'Up')

all.nvp.sigUp <- svmp.nvp.sigUp %>% 
  bind_rows(pla2.nvp.sigUp,svsp.nvp.sigUp) %>% 
  mutate(percent_nvps = no_transcription_factor_search_genes / total_no_search_genes) %>% 
  mutate(enriched_in_nvps = T)

p1 <- ggplot(all.nvp.sigUp,aes(x=percent_nvps,y=reorder(transcription_factor_name,percent_nvps),fill=nvp_family)) +
  geom_bar(stat='identity',show.legend = F) +
  facet_grid(rows = vars(nvp_family),space = 'free',scales = 'free_y') +
  ylab('Enriched TFBS in Non-Venom Paralog Promoters') +
  xlab('Percent of Non-venom Paralogs with one or more binding site') +
  theme_linedraw()



# Read in bound TFBS from promoter analyses

bound.tfbs <- read_tsv('analysis/5_promoters/_tfbs_analysis/allBoundTFBS_Promoters_08.24.21.tsv')

pla2.nvp.scan <- read_csv('analysis/8_NVP_Analyses/promoters/pla2NVP_prom_ciiiderRes/Promoter Panel.csv') %>% janitor::clean_names() %>% mutate(group = 'PLA2')
svmp.nvp.scan <- read_csv('analysis/8_NVP_Analyses/promoters/svmpNVP_prom_ciiiderRes/Promoter Panel.csv') %>% janitor::clean_names() %>% mutate(group = 'SVMP')
svsp.nvp.scan <- read_csv('analysis/8_NVP_Analyses/promoters/svspNVP_prom_ciiiderRes/Promoter Panel.csv') %>% janitor::clean_names() %>% mutate(group = 'SVSP')

all.nvp.scan <- pla2.nvp.scan %>% 
  bind_rows(svmp.nvp.scan,svsp.nvp.scan) %>% 
  select(group,transcription_factor_name) %>% 
  unique() %>% 
  mutate(seen_in_nvps = T)

all.nvp.sigUp.slim <- all.nvp.sigUp %>% select(transcription_factor_name,nvp_family,enriched_in_nvps)

combined.tfbs <- bound.tfbs %>% 
  left_join(all.nvp.scan,by=c('tfbs_id'='transcription_factor_name','vgfamily'='group')) %>% 
  left_join(all.nvp.sigUp.slim,by=c('tfbs_id'='transcription_factor_name','vgfamily'='nvp_family'))


p2 <- ggplot(combined.tfbs,aes(x=vgfamily,y=reorder(tfbs_id,desc(tfbs_id)),fill=vgfamily)) +
  geom_point(data=subset(combined.tfbs,seen_in_nvps==T),size=4,color='black')+
  geom_point(data=subset(combined.tfbs,enriched_in_nvps==T),size=4,color='red')+
  scale_shape_manual(values = c('Primary'=21,'Secondary'=23)) +
  ylab('TFBS') + xlab('Non-venom paralog group') +
  theme_linedraw()
  

p1 + p2 

