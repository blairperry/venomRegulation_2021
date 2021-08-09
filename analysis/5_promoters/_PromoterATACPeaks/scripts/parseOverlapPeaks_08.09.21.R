
library(tidyverse)


overlaps <- read_tsv('./analysis/5_promoters/_PromoterATACPeaks/Overlap_FocalATACpeaks_PromWindows_08.09.21.txt',col_names = F) %>% 
  select(-X4) %>% 
  group_by(X7) %>% 
  mutate(multis = n())

overlaps.multi <- overlaps %>% filter(multis > 1)

# Number of genes w/ more than 1 ATAC-seq peak in promoter window
length(unique(overlaps.multi$X7))

overlaps.multi.plus <- overlaps.multi %>% 
  mutate(strand = str_split_fixed(X8,'[@]',2)[,1]) %>% 
  filter(strand == '#+') %>% 
  filter(X2 == max(X2))%>% 
  select(-multis,-strand)
  
overlaps.multi.minus <- overlaps.multi %>% 
  mutate(strand = str_split_fixed(X8,'[@]',2)[,1]) %>% 
  filter(strand == '#-') %>% 
  filter(X2 == min(X2)) %>% 
  select(-multis,-strand)


overlaps.filtered <- overlaps %>% 
  filter(!(X8 %in% overlaps.multi$X8)) %>% 
  select(-multis) %>% 
  bind_rows(overlaps.multi.plus) %>% 
  bind_rows(overlaps.multi.minus)

overlaps.filtered.formatted <- overlaps.filtered %>% 
  ungroup() %>% 
  mutate(id = str_split_fixed(X8,';',4)[,3] %>% str_remove_all('Crovir_Transcript_ID=')) %>% 
  mutate(strand = str_split_fixed(X8,'[@]',2)[,1] %>% str_remove_all('[#]')) %>% 
  mutate(score = '.') %>% 
  select(X1,X2,X3,id,score,strand)

head(overlaps.filtered.formatted)

write_tsv(overlaps.filtered.formatted,'./analysis/5_promoters/_PromoterATACPeaks/Overlap_ATACpeaks_PromWindows_08.09.21.filtered.bed',col_names = F)
