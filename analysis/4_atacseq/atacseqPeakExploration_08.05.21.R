
library(tidyverse)
library(ggvenn)

# a <- list(`Set 1` = c(1, 3, 5, 7, 9),
#           `Set 2` = c(1, 5, 9, 13),
#           `Set 3` = c(1, 2, 8, 9),
#           `Set 4` = c(6, 7, 10, 12))
# ggvenn(a, c("Set 1", "Set 2"))            # draw two-set venn
# ggvenn(a, c("Set 1", "Set 2", "Set 3"))   # draw three-set venn
# ggvenn(a)   # without set names, the first 4 elements in list will be chose to draw four-set venn

mergedPeaks <- read_tsv('analysis/4_atacseq/peak_regions/merged_peaks/mergedPeakIntersects.txt',col_names = F) %>% 
  mutate(peakID = paste(X1,X2,X3,sep = '_')) %>% 
  select(peakID, sample = 8) %>% 
  mutate(sample = str_split_fixed(sample,'[_]',2)[,1])

rvg1.peaks <- mergedPeaks %>% 
  filter(sample == 'RVG1') %>% 
  select(peakID)

rvg4.peaks <- mergedPeaks %>% 
  filter(sample == 'RVG4') %>% 
  select(peakID)

rvg11.peaks <- mergedPeaks %>% 
  filter(sample == 'RVG11') %>% 
  select(peakID)

peak.overlap <- list('RVG1'=rvg1.peaks$peakID,'RVG4'=rvg4.peaks$peakID,'RVG11'=rvg11.peaks$peakID)

ggvenn(peak.overlap)


focalPeaks <- mergedPeaks %>% 
  group_by(peakID) %>% 
  add_tally() %>% 
  nest(sample_nest=sample) %>% 
  mutate(samples = map_chr(sample_nest, ~ .[[1]] %>% str_c(collapse = "_"))) %>% 
  select(peakID,samples) %>% 
  filter(str_detect(samples,'[_]')) %>% 
  separate(peakID,into = c('chr','start','end'),sep = '[_]')

head(focalPeaks)

write_tsv(focalPeaks,'analysis/4_atacseq/peak_regions/_atacPeaks_2orMoreSamples_08.16.21.tsv')
