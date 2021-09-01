
library(tidyverse)

atac.peaks <- read_tsv('./analysis/4_atacseq/peak_blacklist/_atacPeaks_2orMoreSamples_08.16.21.simple.SCORED.max.bed',col_names = c('chr','start','end','RVG1_density','RVG4_density','RVG11_density')) %>% 
  pivot_longer(4:6,names_to = 'Sample',values_to = 'Density')

head(atac.peaks)

ggplot(atac.peaks,aes(x=Density+1)) +
  geom_histogram(bins = 100) +
  facet_wrap(~Sample) +
  scale_x_log10() +
  theme_classic()

atac.peaks.max1000 <- atac.peaks %>% filter(Density > 1000) %>% select(chr,start,end) %>% unique()

write_tsv(atac.peaks.max1000,'./analysis/4_atacseq/peak_blacklist/Blacklist_MaxDensityOver1k.bed',col_names = F)
