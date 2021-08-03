
#install.packages('ggalluvial')
#install.packages('networkD3')

library(ggalluvial)
library(dplyr)
library(networkD3)
library(stringr)
library(stringi)
library(ggrepel)

SE.Gene <- read.table('./SEtoGene_Association/Cvv_SE_NearestTSS_simpleInfo.bed',sep='\t',stringsAsFactors = F)
withinTAD <- read.table('SEtoGene_Association/Cvv_SEassocPairs_ContainedWithinTAD_12.10.18.bed',sep='\t',stringsAsFactors = F)
withinTAD[,5:6] <- str_split_fixed(withinTAD$V5, "@", 2)


new.df <- SE.Gene
new.df$Arrangement <- ifelse(new.df$V9 == 0, 'Overlapping','Adjacent')
new.df$withinTAD <- ifelse(new.df$V8 %in% withinTAD$V6, 'Within','Outside')
new.df$annotation <- ifelse(grepl('trna',new.df$V8),'tRNA',ifelse(grepl('fgenesh',new.df$V8),'Venom','Non-venom\nProtein Coding'))

venomSE <- new.df[which(new.df$annotation == 'Venom'),]
venomSE <- unique(venomSE$V4)

write.table(venomSE,'VenomSE_IDs.txt',quote = F, col.names=F,row.names=F)


new.df$joined <- paste(new.df$Arrangement,'_',new.df$withinTAD,'_',new.df$annotation)

#Count frequency of each combination (i.e. # Up-Up, # Up-Down, etc.)
sigcount.data <- count(new.df, joined)

#Make alluvial plot input dataframe
se.alluv.data <- as.data.frame(str_split_fixed(sigcount.data$joined, "_", 3))
se.alluv.data$freq <- sigcount.data$n


se.plot <- ggplot(as.data.frame(se.alluv.data),
                  aes(weight = freq, axis1 = V1, axis2 = V2, axis3 = V3)) +
  geom_alluvium(aes(fill = V3),width = 1/6,reverse=F) +
  geom_stratum(width = 1/6, fill = "grey", color = "black",reverse=F) +
  geom_label(stat = "stratum", label.strata = TRUE,reverse=F) +
  scale_x_discrete(limits = c('SE-Gene Arrangement','TAD','Annotation')) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ylab(label='# of SE-Gene Associated Pairs')+
  theme(legend.position = 'none',axis.line = element_line(colour = "black"),panel.background = element_blank())

se.plot

### Sankey Diagram

