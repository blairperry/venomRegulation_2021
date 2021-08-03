

library(tidyverse)
library(stringr)

setwd("~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/ChIPseq/SuperEnhancers/SEtoGene_Association/Characterization")

SE.Gene <- read.table('../Cvv_SE_NearestTSS_simpleInfo.bed',sep='\t',stringsAsFactors = F)
withinTAD <- read.table('../Cvv_SEassocPairs_ContainedWithinTAD_12.10.18.bed',sep='\t',stringsAsFactors = F)
withinTAD[,5:6] <- str_split_fixed(withinTAD$V5, "@", 2)


new.df <- SE.Gene
new.df$Arrangement <- ifelse(new.df$V9 == 0, 'Overlapping','Adjacent')
new.df$withinTAD <- ifelse(new.df$V8 %in% withinTAD$V6, 'Within','Outside')
new.df$annotation <- ifelse(grepl('trna',new.df$V8),'tRNA',ifelse(grepl('fgenesh',new.df$V8),'Venom','Non-venom\nProtein Coding'))

# Counting manually form ChIPseq/SuperEnhancers/updatedVenomSE_DirectOverlap.txt
num.venomSE <- 18
num.venomNotSE <- 36 - 18

num.GenomeSE <- 946
num.GenomeNotSE <- 18538 - 946

venSE_FisherMatrix <- matrix(c(num.venomSE,num.GenomeSE,num.venomNotSE,num.GenomeNotSE),nrow=2)
fisher.test(venSE_FisherMatrix)


prop.df <- tribble(
  ~Type, ~InSE,  ~Prop,
  "Venom", 'Within SE',  num.venomSE/(num.venomSE+num.venomNotSE)*100,
  "Venom", 'Outside SE',  num.venomNotSE/(num.venomSE+num.venomNotSE)*100,
  "All",'Within SE',num.GenomeSE/(num.GenomeSE+num.GenomeNotSE)*100,
  "All",'Outside SE',num.GenomeNotSE/(num.GenomeSE+num.GenomeNotSE)*100,
)


ggplot(prop.df,aes(x=Type,fill=InSE,y=Prop))+
  geom_bar(stat='identity',show.legend = T) +
  scale_fill_manual(values = c('grey70','#FF8B00'))+
  scale_y_continuous(expand=c(0,0))+
  ylab('Percent of Genes')+
  xlab('')+
  theme_half_open()+
  theme(legend.position = 'bottom',legend.title = element_blank())
