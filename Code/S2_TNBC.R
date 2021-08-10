library(ggplot2)
library(fgsea)
library(ggrepel)
library(rPanglaoDB)
library(Seurat)
library(harmony)
library(Matrix)
library(limma)
library(DESeq2)
library(patchwork)


### PCS
load('~/breastCancer/preData.RData')
iData <- iData[,iData$CellTypes %in% c('Epithelial cells I', 'Epithelial cells II')]
newIdent <- paste0(iData$orig.ident, '_', iData$CellTypes)
pseudoCounts <- sapply(unique(newIdent), function(X){
  rowSums(iData@assays$RNA@counts[,newIdent %in% X])
})
pseudoCounts <- pseudoCounts[rowSums(pseudoCounts) > 0,]

tnbc <- which(grepl('GSM', colnames(pseudoCounts)) & grepl('cells II$', colnames(pseudoCounts)))
h <- which(grepl('sc5r', colnames(pseudoCounts)) & grepl('cells I$', colnames(pseudoCounts)))
eData <- pseudoCounts[,c(tnbc, h)]
dData <- data.frame(status = ifelse(grepl('GSM', colnames(eData)), 'C', 'H'))
tData <- DESeqDataSetFromMatrix(eData, dData, design = ~ status)
tData <- DESeq(tData)
tData <- as.data.frame(results(tData))
tData <- tData[complete.cases(tData),]
tData$G <- rownames(tData)
tData <- tData[order(abs(tData$log2FoldChange), decreasing = TRUE),]
tData$G[101:nrow(tData)] <- NA
tData$color <- 'black'
tData$color[((tData$log2FoldChange > 1) & (tData$padj < 0.05))] <- 'red'
tData$color[((tData$log2FoldChange < -1) & (tData$padj < 0.05))] <- 'blue'

pcFC <- tData$log2FoldChange
names(pcFC) <- rownames(tData)

# SC
DE <- read.csv('~/L1000-TNBC/Data/de_EC_TNBC-H.csv')
DE$p_val_adj[DE$p_val_adj == 0] <- min(DE$p_val_adj[DE$p_val_adj != 0])
DE$G <- DE$X
DE <- DE[order(abs(DE$avg_log2FC), decreasing = TRUE),]
DE$G[101:nrow(DE)] <- NA
DE$color <- 'black'
DE$color[(DE$avg_log2FC > 1) & (DE$p_val_adj < 0.05)] <- 'red'
DE$color[(DE$avg_log2FC < -1) & (DE$p_val_adj < 0.05)] <- 'blue'

scFC <- DE$avg_log2FC
names(scFC) <- DE$X

iGenes <- intersect(names(scFC), names(pcFC))

# Bulk
load('~/BRCA.RData')
selectedSamples <- data@colData$paper_BRCA_Subtype_PAM50 %in% c('Normal', 'Basal')
eData <- data@assays@data$`HTSeq - FPKM-UQ`[,selectedSamples]
eData <- as.matrix(eData)
eData <- apply(eData,2,as.integer)
rownames(eData) <- data@rowRanges$external_gene_name
colnames(eData) <- data@colData$barcode[selectedSamples]
eData <- eData[complete.cases(eData),]
eData <- eData[intersect(rownames(eData), iGenes),]
lData <- data@colData$paper_BRCA_Subtype_PAM50[selectedSamples]

design <- cbind(h_vs_tnbc = lData)
rownames(design) = colnames(eData)
eData <- DESeqDataSetFromMatrix(round(eData), colData = design, design = ~ h_vs_tnbc)
eData <- DESeq(eData)

oFit <- as.data.frame(results(eData))
oFit$log2FoldChange <- -1 * oFit$log2FoldChange
oFit <- oFit[complete.cases(oFit),]
oFit$G <- rownames(oFit)
oFit <- oFit[order(abs(oFit$log2FoldChange), decreasing = TRUE),]
oFit$G[101:nrow(oFit)] <- NA
oFit$color <- 'black'
oFit$color[((oFit$log2FoldChange > 1) & (oFit$padj < 0.05))] <- 'red'
oFit$color[((oFit$log2FoldChange < -1) & (oFit$padj < 0.05))] <- 'blue'


bulkFC <- oFit$log2FoldChange
names(bulkFC) <- rownames(oFit)

# Plots
iGenes <- intersect(intersect(names(pcFC),names(scFC)), names(bulkFC))
DE <- DE[DE$X %in% iGenes,]
tData <- tData[rownames(tData) %in% iGenes,]
oFit <- oFit[rownames(oFit) %in% iGenes,]
tData$G <- rownames(tData)
tData$G[101:nrow(tData)] <- NA

P1A <- ggplot(DE, aes(avg_log2FC, -log10(p_val_adj), label = G)) + 
  geom_point(pch = 16, alpha = 0.5, color = DE$color) + 
  theme_bw() +
  geom_text_repel(min.segment.length = 0, fontface = 'italic', size = 3.5, max.overlaps = 10) + 
  xlab(log[2]~(Fold-Change~Single-Cell~RNA-seq)) +
  ylab(-log[10]~(P-value)) +
  labs(tag = 'A', title = 'TNBC - Healthy', subtitle = 'Counts from Single-Cell RNA-seq\n11219 TNBC vs. 11553 Healthy Cells') +
  theme(plot.title = element_text(face = 2))

P1B <- ggplot(tData, aes(log2FoldChange, -log10(pvalue), label = G)) + 
  geom_point(pch = 16, alpha = 0.5, color = tData$color) + 
  xlim(-10,10) + 
  theme_bw() + 
  geom_text_repel(min.segment.length = 0, fontface = 'italic', size = 3.5, max.overlaps = 25) +
  xlab(log[2]~(Fold-Change~PseudoCounts)) +
  ylab(-log[10]~(P-value)) +
  labs(title = 'TNBC - Healthy', subtitle = 'PseudoCounts from Single-Cell RNA-seq\n9 TNBC vs. 13 Healthy Samples') +
  theme(plot.title = element_text(face = 2))

P1C <- ggplot(oFit, aes(log2FoldChange, -log10(pvalue), label = G)) + 
  geom_point(pch = 16, alpha = 0.5, color = oFit$color) + 
  xlim(-10,10) + 
  theme_bw() + 
  geom_text_repel(min.segment.length = 0, fontface = 'italic', size = 3.5, max.overlaps = 20) +
  xlab(log[2]~(Fold-Change~Bulk~RNA-seq)) +
  ylab(-log[10]~(P-value)) +
  labs(title = 'TNBC - Healthy', subtitle = 'Counts from BRCA-TCGA RNA-seq\n194 TNBC vs. 40 Healthy Samples') +
  theme(plot.title = element_text(face = 2))

df <- data.frame(sc = scFC[iGenes], pc = pcFC[iGenes], bulk = bulkFC[iGenes])
df <- df[complete.cases(df),]

o <- predict(lm(pc~sc, df), newdata = data.frame(sc=df$sc), interval = 'prediction', level = 0.95)
df$pc_sc_lwr <- o[,2]
df$pc_sc_upr <- o[,3]

o <- predict(lm(bulk~sc, df), newdata = data.frame(sc=df$sc), interval = 'prediction', level = 0.95)
df$bulk_sc_lwr <- o[,2]
df$bulk_sc_upr <- o[,3]

corValue <- round(cor(df$sc, df$pc, method = 'sp'),3)
df$g <- rownames(df)
df$g[!(((df$pc < -1) & (df$sc < -1)) | ((df$pc > 1) & (df$sc > 1)))] <- NA
df$color <- 'black'
df$color[df$sc > 1 & df$pc > 1] <- 'red'
df$color[df$sc < -1 & df$pc < -1] <- 'blue'
P1D <- ggplot(df, aes(sc, pc, label = g)) + 
  geom_point(pch = 16, alpha = 0.5, color = df$color) + 
  geom_density2d() + 
  geom_smooth(method = 'lm', se = FALSE, col = 'red') +
  geom_line(aes(sc, pc_sc_lwr), data = df, col = 'red', lty = 2) +
  geom_line(aes(sc, pc_sc_upr), data = df, col = 'red', lty = 2) +
  labs(tag = 'B', title = 'PseudoCounts vs. Single-Cell', subtitle = parse(text = paste0('rho == ', corValue))) +
  ylim(-10,10) +
  geom_text_repel(min.segment.length = 0, fontface = 'italic', size = 3.5) + 
  theme_bw() +
  xlab(log[2]~(Fold-Change~Single-Cell~RNA-seq)) +
  ylab(log[2]~(Fold-Change~PseudoCounts)) +
  theme(plot.title = element_text(face = 2))

corValue <- round(cor(df$sc, df$bulk, method = 'sp'),3)
df$g <- rownames(df)
df$g[!(((df$bulk < -1) & (df$sc < -1)) | ((df$bulk > 1) & (df$sc > 1)))] <- NA
df$color <- 'black'
df$color[df$sc > 1 & df$bulk > 1] <- 'red'
df$color[df$sc < -1 & df$bulk < -1] <- 'blue'
P1E <- ggplot(df, aes(sc, bulk, label = g)) + 
  geom_point(pch = 16, alpha = 0.5, color = df$color) + 
  geom_density2d() + 
  geom_smooth(method = 'lm', se = FALSE, col = 'red') +
  geom_line(aes(sc, bulk_sc_lwr), data = df, col = 'red', lty = 2) +
  geom_line(aes(sc, bulk_sc_upr), data = df, col = 'red', lty = 2) +
  labs(title = 'Bulk vs. Single-Cell', subtitle = parse(text = paste0('rho == ', corValue))) +
  geom_text_repel(min.segment.length = 0, fontface = 'italic', size = 3.5) + 
  theme_bw() +
  xlab(log[2]~(Fold-Change~Single-Cell~RNA-seq)) +
  ylab(log[2]~(Fold-Change~BRCA-TCGA~RNA-seq)) +
  theme(plot.title = element_text(face = 2))

png('../Figures/F2.png', width = 4800, height = 2400, res = 300)
(P1A | P1B | P1C)/(P1D | P1E)
dev.off()
