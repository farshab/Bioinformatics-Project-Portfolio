read_count <- read.table("Coundata.csv", header = T)
head(read_count)
dim(read_count)
metadata<-read.table("ColData.csv", header=T)
head(metadata)
rownames(metadata)
colnames(read_count)
all(rownames(metadata) == colnames(read_count))
match(rownames(metadata),colnames(read_count))
idx<-match(rownames(metadata),colnames(read_count))
reordered_metadata<-metadata[idx,]
View(reordered_metadata)
all(rownames(reordered_metadata) == colnames(read_count))
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
dds<-DESeqDataSetFromMatrix(countData = read_count, colData = reordered_metadata, design = ~1)
dds<-estimateSizeFactors(dds)
sizeFactors(dds)
plot(sizeFactors(dds),colSums(counts(dds)))
all(rownames(reordered_metadata) == colnames(read_count))
abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))
normlzd_dds <- counts(dds, normalized=T)
head(normlzd_dds)
plot(hclust(dist(t(normlzd_dds))), labels=colData(dds)$protocol)
plot(hclust(dist(t(normlzd_dds))), labels=colData(dds)$Time)
plot(log(normlzd_dds[,1])+1, log(normlzd_dds[,2])+1, cex =.1)
vsd <- vst(dds, blind = T)
vsd_mat <- assay(vsd)
vsd_cor <- cor(vsd_mat)
vsd_cor
pheatmap(vsd_cor)
plotPCA(vsd, intgroup="protocol")
plotPCA(vsd, intgroup="Time")
dds <- DESeqDataSetFromMatrix(countData = read_count, colData = reordered_metadata, design = ~ Time + protocol)
dds<DESeq(dds)
res <-results(dds)
library(DESeq2)
res <-results(dds)
summary(res)
res
plotMA(res, ylim=c(-5,5) )
resBigFC <- results(dds, lfcThreshold=1, altHypothesis="greaterAbs")
resBigFC
plotMA(resBigFC, ylim=c(-5,5))
abline(h=c(-1,1),lwd=5)
plotDispEsts(dds)
resSort <- res[order(res$pvalue),]
head(resSort)
library(org.Rn.eg.db)
keytypes(org.Rn.eg.db)
head(rownames(dds))
head(resSort,n=10)
geneinfo
geneinfo <- select(org.Rn.eg.db, keys=rownames(resSort)[1:20], columns=c("ENSEMBL", "SYMBOL","GENENAME"), keytype="ENSEMBL")
geneinfo
res_sig <- data.frame(normlzd_dds[geneinfo$ENSEMBL, ])
pheatmap(res_sig,color = heat_colors, cluster_rows = T, show_rownames = T, annotation = dplyr::select(metadata, protocol), scale = "row")
heat_colors <- brewer.pal(6, "YlOrRd")
pheatmap(res_sig,color = heat_colors, cluster_rows = T, show_rownames = T, annotation = dplyr::select(metadata, protocol), scale = "row")
top_20 <- gather(top_20, key = "samplename", value = "normalized_counts", 2:8)
top_20 <- inner_join(top_20, rownames_to_column(reordered_metadata, var  = "samplename"), by = "samplename")
ggplot(top_20) + geom_point(aes(x = ensgene, y = normalized_counts, color = protocol)) + scale_y_log10() + xlab("Genes") + ylab("Normalized Counts") + ggtitle("Top 20 Significant DE Genes") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(plot.title = element_text(hjust = 0.5))
ggplot(top_20) + geom_point(aes(x = ensgene, y = normalized_counts, color = protocol)) + scale_y_log10() + xlab("Genes") + ylab("Normalized Counts") + ggtitle("Top 20 Significant DE Genes") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(plot.title = element_text(hjust = 0.5))

vsd <- vst(dds, blind = T, show_rownames=T)
vsd_mat <- assay(vsd)
vsd_cor <- cor(vsd_mat)
vsd_cor
pheatmap(vsd_cor)

plot(sizeFactors(dds),colSums(counts(dds)))

padj.cutoff <- 0.05
lfc.cutoff <- 0.58
threshold <- res$padj < padj.cutoff & abs(res$log2FoldChange) > lfc.cutoff
res$threshold <- threshold
sigOE <- data.frame(subset(res, threshold==TRUE))
sigOE_ordered <- sigOE[order(sigOE$padj), ]
top20_sigOE <- rownames(sigOE_ordered[1:20, ])
normalized_counts <- counts(dds, normalized=T)
top_20_sig0E_norm <- normlzd_dds[top20_sigOE, ]
sigOE_ordered
top20_sigOE
top



view(top_20)
view(top_20_sig0E_norm)

?melt

install.packages("reshape")
melted_top20_sigOE <- data.frame(melt(top_20_sig0E_norm))
view(melted_top20_sigOE)
colnames(melted_top20_sigOE) <- c("gene", "samplename", "normalized_counts", "protocol")
view(melted_top20_sigOE)
meta$samplename <- rownames(meta)
melted_top20_sigOE <- merge(melted_top20_sigOE, meta)
ggplot(melted_top20_sigOE) +
  geom_point(aes(x = gene, y = normalized_counts, color = samplename)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title=element_text(hjust=0.5))
