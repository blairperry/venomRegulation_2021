
library(tidyverse)

stringNet <- read_tsv('./analysis/3_transcription_factors/_functional_characterization/stringdb_ExpCurOnly_08.03.21/string_interactions.tsv') %>% janitor::clean_names()

stringNet.ERKInteractors <- stringNet %>% filter(number_node1 == 'MAPK1' | node2 == 'MAPK1')

ERK.1stDegreeTFs <- union(stringNet.ERKInteractors$number_node1,stringNet.ERKInteractors$node2)

stringNet.AdditionalInteractors <- stringNet %>% filter(number_node1 %in% ERK.1stDegreeTFs | node2 %in% ERK.1stDegreeTFs)

ERK.AllInteractors <- as_tibble(union(union(ERK.1stDegreeTFs,stringNet.AdditionalInteractors$number_node1),stringNet.AdditionalInteractors$node2)) %>% 
  select(TF=1) %>% 
  mutate(degree = ifelse(TF %in% ERK.1stDegreeTFs,1,2))

# write_tsv(ERK.AllInteractors,'./analysis/3_transcription_factors/_functional_characterization/stringdb_ExpCurOnly_08.03.21/ERK_1st2ndDegreeTFs_08.03.21.tsv')


stringNet.1st2ndFiltered <- stringNet %>% filter(number_node1 %in% ERK.AllInteractors$TF & node2 %in% ERK.AllInteractors$TF)

write_tsv(stringNet.1st2ndFiltered,'./analysis/3_transcription_factors/_functional_characterization/stringdb_ExpCurOnly_08.03.21/string_interactions.1st2ndDegree.tsv')
