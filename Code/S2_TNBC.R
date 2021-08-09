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

# Sample types
tnbc <- which(grepl('GSM', colnames(pseudoCounts)) & grepl('cells II$', colnames(pseudoCounts)))
h <- which(grepl('sc5r', colnames(pseudoCounts)) & grepl('cells I$', colnames(pseudoCounts)))
eData <- pseudoCounts[,c(tnbc, h)]
dData <- data.frame(status = ifelse(grepl('GSM', colnames(eData)), 'C', 'H'))
tData <- DESeqDataSetFromMatrix(eData, dData, design = ~ status)
tData <- DESeq(tData)
tData <- as.data.frame(results(tData))
tData <- tData[complete.cases(tData),]
tData$G <- rownames(tData)
tData <- tData[order(-log10(tData$pvalue), abs(tData$log2FoldChange), decreasing = TRUE),]
tData$G[101:nrow(tData)] <- NA

P1A <- ggplot(tData, aes(log2FoldChange, -log10(pvalue), label = G)) + 
  geom_point(pch = 16, alpha = 0.5) + 
  xlim(-10,10) + 
  theme_bw() + 
  geom_text_repel() +
  xlab(log[2]~(Fold-Change~PseudoCounts)) +
  ylab(-log[10]~(P-value))

pcFC <- tData$log2FoldChange
names(pcFC) <- rownames(tData)

# SC
DE <- read.csv('~/L1000-TNBC/Data/de_EC_TNBC-H.csv')
scFC <- DE$avg_log2FC
names(scFC) <- DE$X
DE$G <- DE$X
DE <- DE[order(abs(DE$avg_log2FC), decreasing = TRUE),]
DE$G[101:nrow(DE)] <- NA

ggplot(DE, aes(avg_log2FC, -log10(p_val), label = G)) + 
  geom_point(pch = 16, alpha = 0.5) + 
  xlim(-10,10) + 
  theme_bw() +
  geom_text_repel() + 
  xlab(log[2]~(Fold-Change~PseudoCounts)) +
  ylab(-log[10]~(P-value))

# Bulk
load('~/BRCA.RData')
selectedSamples <- data@colData$paper_BRCA_Subtype_PAM50 %in% c('Normal', 'Basal')
eData <- data@assays@data$`HTSeq - FPKM-UQ`[,selectedSamples]
eData <- as.matrix(eData)
eData <- apply(eData,2,as.integer)
rownames(eData) <- data@rowRanges$external_gene_name
colnames(eData) <- data@colData$barcode[selectedSamples]
eData <- eData[complete.cases(eData),]
lData <- data@colData$paper_BRCA_Subtype_PAM50[selectedSamples]

design <- cbind(h_vs_tnbc = lData)
rownames(design) = colnames(eData)
eData <- DESeqDataSetFromMatrix(round(eData), colData = design, design = ~ h_vs_tnbc)
eData <- DESeq(eData)

oFit <- as.data.frame(results(eData))

bulkFC <- -1 * oFit$log2FoldChange
names(bulkFC) <- rownames(oFit)

# Plots
iGenes <- intersect(intersect(names(pcFC),names(scFC)), names(bulkFC))
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
PA <- ggplot(df, aes(sc, pc, label = g)) + 
  geom_point(pch = 16, alpha = 0.5) + 
  geom_density2d() + 
  geom_smooth(method = 'lm', se = FALSE, col = 'red') +
  geom_line(aes(sc, pc_sc_lwr), data = df, col = 'red', lty = 2) +
  geom_line(aes(sc, pc_sc_upr), data = df, col = 'red', lty = 2) +
  labs(subtitle = parse(text = paste0('rho == ', corValue))) +
  ylim(-10,10) +
  geom_text_repel(min.segment.length = 0, fontface = 'italic') + 
  theme_bw() +
  xlab(log[2]~('Fold-Change Single-cell RNA-seq')~'by MAST') +
  ylab(log[2]~('Fold-Change Pseudocounts')~'by DESeq2')

corValue <- round(cor(df$sc, df$bulk, method = 'sp'),3)
df$g <- rownames(df)
df$g[!(((df$bulk < -1) & (df$sc < -1)) | ((df$bulk > 1) & (df$sc > 1)))] <- NA
PB <- ggplot(df, aes(sc, bulk, label = g)) + 
  geom_point(pch = 16, alpha = 0.5) + 
  geom_density2d() + 
  geom_smooth(method = 'lm', se = FALSE, col = 'red') +
  geom_line(aes(sc, bulk_sc_lwr), data = df, col = 'red', lty = 2) +
  geom_line(aes(sc, bulk_sc_upr), data = df, col = 'red', lty = 2) +
  labs(subtitle = parse(text = paste0('rho == ', corValue))) +
  geom_text_repel(min.segment.length = 0, fontface = 'italic') + 
  theme_bw() +
  xlab(log[2]~('Fold-Change Single-cell RNA-seq')~'by MAST') +
  ylab(log[2]~('Fold-Change Bulk TCGA RNA-seq ')~'by DESeq2')

PA | PB
