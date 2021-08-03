
# BiocManager::install("DESeq2")
# BiocManager::install("IHW")

library(DESeq2)
library(pheatmap)
library(viridis)
library(IHW)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(tidyverse)

tx2gene <- read.csv('./analysis/1_gene_expression/tx2gene_v4.csv',header = F,stringsAsFactors = F)

rawCounts <- read.table('./data/rnaseq/raw_counts/cvv_VenomRegulationRNAseq_wNonVenom_rawCounts_08.02.21.txt',stringsAsFactors = F,skip=1,header=T)
rawCounts <- rawCounts[,c(-2,-3,-4,-5)]

rawCounts.simple <- rawCounts[which(!is.na(rawCounts$Crovir_Transcript_ID)),]

rawCounts.simple <- rawCounts.simple[,c(-1,-2)]

rawCounts.simple <- merge(rawCounts.simple,tx2gene,by=1,all.x=T)

row.names(rawCounts.simple) <- rawCounts.simple$V2
rawCounts.simple <- rawCounts.simple[,c(-1,-30)]

rawCounts.simple <- rawCounts.simple %>% select(-8,-13,-14,-15,-16,-18,-4,-5,-6,-19,-20,-21,-22,-23,-24,-25,-28,-27,-26)
colnames(rawCounts.simple) <- colnames(rawCounts.simple) %>% str_remove_all('..STAR_mapped.|Aligned.sortedByCoord.out.bam')

rawCounts.simple <- rawCounts.simple[rowSums( rawCounts.simple != 0 ) >= 3,]

test <- rawCounts.simple %>% rownames_to_column()


colnames(rawCounts.simple)

#################################################################################### INCLUDE JUST 1 TECH REP


condition <- factor(c('NonVen','NonVen','Ven','NonVen','NonVen','NonVen','Ven','Ven','NonVen'))
batch <- factor(c('2','2','1','1','2','2','2','2','1'))

colData <- DataFrame(condition = condition,batch = batch)


dds <- DESeqDataSetFromMatrix(rawCounts.simple,colData,formula(~batch + condition))

dds <- DESeq(dds)

normcounts <- counts(dds,normalized=TRUE)
#write.csv(normcounts,'./normalized_counts/CvvVenomReg_RepRNAseq_NormCounts_07.28.21.csv',)


vsd.normCounts <- as.data.frame(assay(vst(dds, blind=FALSE)))

# write.csv(vsd.normCounts,'./analysis/1_gene_expression/norm_counts/CvvVenomReg_RepRNAseq_wNonVen_VSTNormCounts_08.02.21.csv')

### Vst PCA
vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c("condition","batch"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, label=name)) +
  geom_point(size=3) +
  geom_text_repel() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + theme_minimal()





###
### Venom vs. NonVenom
###

deResVenomVsNonVen <- as.data.frame(results(dds, contrast=c('condition','Ven','NonVen')))

test <- as.data.frame(assays(dds)[["cooks"]]) %>% rownames_to_column()
summary(results(dds, cooksCutoff = 1, contrast=c('condition','Ven','NonVen')))


ihwRes <- ihw(pvalue ~ baseMean,  data = deResVenomVsNonVen, alpha = 0.05)

rejections(ihwRes)
deRes <- deResVenomVsNonVen
deRes <- na.omit(deRes)

plot(ihwRes)

gg <- ggplot(as.data.frame(ihwRes), aes(x = pvalue, y = adj_pvalue, col = group)) +
  geom_point(size = 0.25) + scale_colour_hue(l = 70, c = 150, drop = FALSE)
gg
gg %+% subset(as.data.frame(ihwRes), adj_pvalue <= 0.2)

ggplot(deRes, aes(x = pvalue)) + geom_histogram(binwidth = 0.025, boundary = 0)

deRes$baseMeanGroup <- groups_by_filter(deRes$baseMean, 10)

ggplot(deRes, aes(x=pvalue)) +
  geom_histogram(binwidth = 0.025, boundary = 0) +
  facet_wrap( ~ baseMeanGroup, nrow = 2)

ggplot(deRes, aes(x = pvalue, col = baseMeanGroup)) + stat_ecdf(geom = "step")

rbind(data.frame(pvalue = deRes$pvalue, covariate = rank(deRes$baseMean)/nrow(deRes),
                 covariate_type="base mean"),
      data.frame(pvalue = deRes$pvalue, covariate = rank(deRes$log2FoldChange)/nrow(deRes),
                 covariate_type="log2 fc")) %>%
  ggplot(aes(x = covariate, y = -log10(pvalue))) + geom_hex(bins = 100) +
  facet_grid( . ~ covariate_type) + ylab(expression(-log[10]~p))

# 
deResVenomVsNonVen$IHW_pvalue <- ihwRes@df$adj_pvalue

deResVenomVsNonVen <- deResVenomVsNonVen[order(deResVenomVsNonVen$IHW_pvalue),]

write.csv(as.data.frame(deResVenomVsNonVen),file='./analysis/1_gene_expression/pairwise_results/cvv_Venom.vs.NonVenom_VG_PairwiseResult_08.02.21.csv',row.names = T)
