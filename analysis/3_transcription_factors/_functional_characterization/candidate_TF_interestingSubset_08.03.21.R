
library(tidyverse)
library(scico)


cand_tfs <- read_tsv('./analysis/3_transcription_factors/allCandidateTFs_08.03.21.tsv') %>% 
  mutate(id = id.y)

## Importing ERK interactions and filtering to 1st degree interactions

erk_interact_all <- read_tsv('./analysis/3_transcription_factors/_functional_characterization/stringdb_ExpCurOnly_08.03.21/ERK_1st2ndDegreeTFs_08.03.21.tsv') %>% 
  filter(degree == 1)

# Import Kegg enrichment analysis/characterization results

kegg_path <- read_tsv('./analysis/3_transcription_factors/_functional_characterization/webgestalt_kegg_08.03.21/Project_wg_result1628014729/enrichment_results_wg_result1628014729.txt')

adren_enrich <- kegg_path %>% 
  filter(str_detect(description,'Adrenergic')) %>% 
  mutate(tfs = str_split(userId,';')) %>% 
  unnest(tfs) %>% 
  select(tfs)

ap1_component <- read_tsv('./analysis/3_transcription_factors/_functional_characterization/stringdb_ExpCurOnly_08.03.21/enrichment.Component.tsv') %>% 
  filter(str_detect(`term description`,'AP-1')) %>% 
  mutate(tfs = str_split(`matching proteins in your network (labels)`,',')) %>% 
  unnest(tfs)

upr <- read_tsv('analysis/3_transcription_factors/_functional_characterization/UPR_GO.0030968_ProtIDs.txt',col_names = 'id')

prev_impl <- read_tsv('./analysis/3_transcription_factors/_functional_characterization/TFs_PrevImplicated_wTigris_02.23.21.txt',col_names = 'id')

cand_tfs.interest <- cand_tfs %>% 
  mutate(prev_implicated = ifelse(id %in% prev_impl$id | id %in% ap1_component$tfs,T,F)) %>% 
  mutate(erk_interact = ifelse(id %in% erk_interact_all$TF,T,F)) %>% 
  mutate(adren_path = ifelse(id %in% adren_enrich$tfs,T,F)) %>% 
  mutate(ap1_complex = ifelse(id %in% ap1_component$tfs,T,F)) %>%
  mutate(upr_pathway = ifelse(id %in% upr$id,T,F)) %>%
  filter(prev_implicated | erk_interact | ap1_complex | adren_path | upr_pathway) %>% 
  select(id, 
         `Upregulated (RNA-seq)` = Upreg.Ven.vs.NonVen, 
         `SE-Associated (ChIP-seq)` = SEassociated, 
         `Previously implicated in venom regulation` = prev_implicated,
         `Interacts with ERK` = erk_interact, 
         `Adrenoceptor Signaling Pathway`=adren_path,
         `Unfolded Protein Response Pathway`=upr_pathway,
         `AP-1 Complex Member` = ap1_complex) %>% 
  pivot_longer(-1,names_to = 'feature',values_to='value') %>% 
  mutate(feature = factor(feature,levels=c("Upregulated (RNA-seq)","SE-Associated (ChIP-seq)","Differentially Bound (ATAC-seq)","Previously implicated in venom regulation","Adrenoceptor Signaling Pathway","AP-1 Complex Member","Interacts with ERK","Unfolded Protein Response Pathway")))

# write_tsv(cand_tfs.interest,'CandTFs_InterestingSubset_FULL_03.20.21.tsv')

cand_tfs.interest %>% 
  filter(str_detect(feature,'ERK'),value==TRUE)


ggplot(cand_tfs.interest,aes(y=reorder(id,desc(id)),x=feature,fill=feature,alpha=value)) +
  geom_point(size=4,show.legend = F,pch=21) +
  scale_alpha_manual(values=c('TRUE'=1,'FALSE'=0)) +
  scale_fill_scico_d() +
  theme_linedraw() +
  xlab('') +
  ylab('Transcription Factor') +
  ggtitle('Subset of Candidate TFs') +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),plot.title.position = 'plot', plot.title = element_text(face='bold'),axis.title = element_text(face = 'bold'))


