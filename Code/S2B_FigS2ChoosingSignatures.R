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
library(statsExpressions)
library(GSVA)
library(ComplexHeatmap)
source('https://raw.githubusercontent.com/dosorio/utilities/master/cheatmap2ggplot.R')

MSigDB_Hallmarks <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020')

### PCS
load('../Data/EpithelialCells.RData')
breastData <- breastData[,Idents(breastData) %in% c('Cancer_TN', 'Healthy_TP')]
newIdent <- paste0(breastData$orig.ident, '_', Idents(breastData))

pseudoCounts <- sapply(unique(newIdent), function(X){
  rowSums(breastData@assays$RNA@counts[,newIdent %in% X])
})
pseudoCounts <- pseudoCounts[rowSums(pseudoCounts) > 0,]

eData <- pseudoCounts
dData <- data.frame(status = ifelse(grepl('Cancer', colnames(eData)), 'C', 'H'))
dData$status <- factor(dData$status, levels = c('C', 'H'))
tData <- DESeqDataSetFromMatrix(eData, dData, design = ~ status)
tData <- DESeq(tData)
tData <- as.data.frame(results(tData))
tData <- tData[complete.cases(tData),]
tData$log2FoldChange <- -1 * tData$log2FoldChange
tData$G <- rownames(tData)
tData <- tData[order(abs(tData$log2FoldChange), decreasing = TRUE),]
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
DE$color <- 'black'
DE$color[(DE$avg_log2FC > 1) & (DE$p_val_adj < 0.05)] <- 'red'
DE$color[(DE$avg_log2FC < -1) & (DE$p_val_adj < 0.05)] <- 'blue'

scFC <- DE$avg_log2FC
names(scFC) <- DE$X

DE2 <- read.csv('~/L1000-TNBC/Data/de_EC_TNBC-allH.csv')
DE2$p_val_adj[DE2$p_val_adj == 0] <- min(DE2$p_val_adj[DE2$p_val_adj != 0])
DE2$G <- DE2$X
DE2 <- DE2[order(abs(DE2$avg_log2FC), decreasing = TRUE),]
DE2$color <- 'black'
DE2$color[(DE2$avg_log2FC > 1) & (DE2$p_val_adj < 0.05)] <- 'red'
DE2$color[(DE2$avg_log2FC < -1) & (DE2$p_val_adj < 0.05)] <- 'blue'

scFC2 <- DE2$avg_log2FC
names(scFC2) <- DE2$X

# Bulk
load('~/BRCA.RData')
selectedSamples <- data@colData$paper_BRCA_Subtype_PAM50 %in% c('Normal', 'Basal')
eData <- data@assays@data$`HTSeq - FPKM-UQ`[,selectedSamples]
eData <- as.matrix(eData)
eData <- apply(eData,2,as.integer)
rownames(eData) <- data@rowRanges$external_gene_name
colnames(eData) <- data@colData$barcode[selectedSamples]
eData <- eData[complete.cases(eData),]
eData <- eData[rownames(eData) %in% intersect(names(pcFC), names(scFC)),]
lData <- data@colData$paper_BRCA_Subtype_PAM50[selectedSamples]
lData[lData == 'Basal'] <- 'C'
lData[lData == 'Normal'] <- 'H'
lData <- factor(lData, levels = c('C', 'H'))
design <- data.frame(h_vs_tnbc = lData)
rownames(design) = colnames(eData)

eData <- DESeqDataSetFromMatrix(eData, colData = design, design = ~ h_vs_tnbc)
eData <- DESeq(eData)

oFit <- as.data.frame(results(eData))
oFit <- oFit[complete.cases(oFit),]
oFit$log2FoldChange <- -1* oFit$log2FoldChange
oFit$G <- rownames(oFit)
oFit <- oFit[order(abs(oFit$log2FoldChange), decreasing = TRUE),]
oFit$color <- 'black'
oFit$color[((oFit$log2FoldChange > 1) & (oFit$padj < 0.05))] <- 'red'
oFit$color[((oFit$log2FoldChange < -1) & (oFit$padj < 0.05))] <- 'blue'

bulkFC <- oFit$log2FoldChange
names(bulkFC) <- rownames(oFit)

g1 <- DE$X[((abs(DE$avg_log2FC) > 1) & (DE$avg_log2FC < 0.05))]
g2 <- rownames(tData)[((abs(tData$log2FoldChange) > 1) & (tData$padj < 0.05))]
g3 <- rownames(oFit)[((abs(oFit$log2FoldChange) > 1) & (oFit$padj < 0.05))]
g4 <- DE2$X[((abs(DE2$avg_log2FC) > 1) & (DE2$avg_log2FC < 0.05))]

sGenes <- c(g1,g2,g3,g4)
sGenes <- table(sGenes)
sGenes <- names(sGenes[sGenes >= 2])

iGenes <- intersect(intersect(intersect(names(pcFC),names(scFC)), names(bulkFC)), names(scFC2))

df <- data.frame(sc = scFC[iGenes], sc2 = scFC2[iGenes], pc = pcFC[iGenes], bulk = bulkFC[iGenes])
df <- df[complete.cases(df),]

o <- predict(lm(bulk~sc, df), newdata = data.frame(sc=df$sc), interval = 'prediction', level = 0.95)
df$bulk_sc_lwr <- o[,2]
df$bulk_sc_upr <- o[,3]

corValue <- statsExpressions::corr_test(df,sc, bulk, type = 'non')$expression[[1]]
df$g <- rownames(df)
df$g[!(((df$bulk < -1) & (df$sc < -1)) | ((df$bulk > 1) & (df$sc > 1)))] <- NA
df$color <- 'black'
df$color[df$sc > 1 & df$bulk > 1] <- 'red'
df$color[df$sc < -1 & df$bulk < -1] <- 'blue'
P1A <- ggplot(df, aes(sc, bulk, label = g)) + 
  geom_point(pch = 16, alpha = ifelse(df$color == 'black', 0.25,1), color = df$color) + 
  geom_abline(slope = 1, intercept = 0, color = 'red', lty = 2) + 
  geom_density2d() + 
  labs(title = 'TN Breast Cancer - TP Healthy', subtitle = corValue) +
  geom_text_repel(min.segment.length = 0, fontface = 'italic', size = 3.5, bg.color = 'white', bg.r = 0.1) + 
  theme_bw() +
  xlab(log[2]~(Fold-Change~Single-Cell~RNA-seq)) +
  ylab(log[2]~(Fold-Change~BRCA-TCGA~RNA-seq)) +
  theme(plot.title = element_text(face = 2), plot.subtitle = element_text(size = 10))
P1A 


# SC
DE <- read.csv('~/L1000-TNBC/Data/de_EC_C-H.csv')
DE$p_val_adj[DE$p_val_adj == 0] <- min(DE$p_val_adj[DE$p_val_adj != 0])
DE$G <- DE$X
DE <- DE[order(abs(DE$avg_log2FC), decreasing = TRUE),]
DE$color <- 'black'
DE$color[(DE$avg_log2FC > 1) & (DE$p_val_adj < 0.05)] <- 'red'
DE$color[(DE$avg_log2FC < -1) & (DE$p_val_adj < 0.05)] <- 'blue'

scFC <- DE$avg_log2FC
names(scFC) <- DE$X

g1 <- DE$X[((abs(DE$avg_log2FC) > 1) & (DE$avg_log2FC < 0.05))]
g2 <- rownames(tData)[((abs(tData$log2FoldChange) > 1) & (tData$padj < 0.05))]
g3 <- rownames(oFit)[((abs(oFit$log2FoldChange) > 1) & (oFit$padj < 0.05))]
g4 <- DE2$X[((abs(DE2$avg_log2FC) > 1) & (DE2$avg_log2FC < 0.05))]

sGenes <- c(g1,g2,g3,g4)
sGenes <- table(sGenes)
sGenes <- names(sGenes[sGenes >= 2])

iGenes <- intersect(intersect(intersect(names(pcFC),names(scFC)), names(bulkFC)), names(scFC2))

df <- data.frame(sc = scFC[iGenes], sc2 = scFC2[iGenes], pc = pcFC[iGenes], bulk = bulkFC[iGenes])
df <- df[complete.cases(df),]

o <- predict(lm(bulk~sc, df), newdata = data.frame(sc=df$sc), interval = 'prediction', level = 0.95)
df$bulk_sc_lwr <- o[,2]
df$bulk_sc_upr <- o[,3]

corValue <- statsExpressions::corr_test(df,sc, bulk, type = 'non')$expression[[1]]
df$g <- rownames(df)
df$g[!(((df$bulk < -1) & (df$sc < -1)) | ((df$bulk > 1) & (df$sc > 1)))] <- NA
df$color <- 'black'
df$color[df$sc > 1 & df$bulk > 1] <- 'red'
df$color[df$sc < -1 & df$bulk < -1] <- 'blue'
P1B <- ggplot(df, aes(sc, bulk, label = g)) + 
  geom_point(pch = 16, alpha = ifelse(df$color == 'black', 0.25,1), color = df$color) + 
  geom_abline(slope = 1, intercept = 0, color = 'red', lty = 2) + 
  geom_density2d() + 
  labs(title = 'Breast Cancer - Healthy', subtitle = corValue) +
  geom_text_repel(min.segment.length = 0, fontface = 'italic', size = 3.5, bg.color = 'white', bg.r = 0.1) + 
  theme_bw() +
  xlab(log[2]~(Fold-Change~Single-Cell~RNA-seq)) +
  ylab(log[2]~(Fold-Change~BRCA-TCGA~RNA-seq)) +
  theme(plot.title = element_text(face = 2), plot.subtitle = element_text(size = 10))
P1B

# SC
DE <- read.csv('~/L1000-TNBC/Data/de_EC_TNBC-TNH.csv')
DE$p_val_adj[DE$p_val_adj == 0] <- min(DE$p_val_adj[DE$p_val_adj != 0])
DE$G <- DE$X
DE <- DE[order(abs(DE$avg_log2FC), decreasing = TRUE),]
DE$color <- 'black'
DE$color[(DE$avg_log2FC > 1) & (DE$p_val_adj < 0.05)] <- 'red'
DE$color[(DE$avg_log2FC < -1) & (DE$p_val_adj < 0.05)] <- 'blue'

scFC <- DE$avg_log2FC
names(scFC) <- DE$X

# Label Selection
g1 <- DE$X[((abs(DE$avg_log2FC) > 1) & (DE$avg_log2FC < 0.05))]
g2 <- rownames(tData)[((abs(tData$log2FoldChange) > 1) & (tData$padj < 0.05))]
g3 <- rownames(oFit)[((abs(oFit$log2FoldChange) > 1) & (oFit$padj < 0.05))]
g4 <- DE2$X[((abs(DE2$avg_log2FC) > 1) & (DE2$avg_log2FC < 0.05))]

sGenes <- c(g1,g2,g3,g4)
sGenes <- table(sGenes)
sGenes <- names(sGenes[sGenes >= 2])

iGenes <- intersect(intersect(intersect(names(pcFC),names(scFC)), names(bulkFC)), names(scFC2))

df <- data.frame(sc = scFC[iGenes], sc2 = scFC2[iGenes], pc = pcFC[iGenes], bulk = bulkFC[iGenes])
df <- df[complete.cases(df),]

o <- predict(lm(bulk~sc, df), newdata = data.frame(sc=df$sc), interval = 'prediction', level = 0.95)
df$bulk_sc_lwr <- o[,2]
df$bulk_sc_upr <- o[,3]


corValue <- statsExpressions::corr_test(df,sc, bulk, type = 'non')$expression[[1]]
df$g <- rownames(df)
df$g[!(((df$bulk < -1) & (df$sc < -1)) | ((df$bulk > 1) & (df$sc > 1)))] <- NA
df$color <- 'black'
df$color[df$sc > 1 & df$bulk > 1] <- 'red'
df$color[df$sc < -1 & df$bulk < -1] <- 'blue'
P1C <- ggplot(df, aes(sc, bulk, label = g)) + 
  geom_point(pch = 16, alpha = ifelse(df$color == 'black', 0.25,1), color = df$color) + 
  geom_abline(slope = 1, intercept = 0, color = 'red', lty = 2) + 
  geom_density2d() + 
  labs(title = 'TN Breast Cancer - TN Healthy', subtitle = corValue) +
  geom_text_repel(min.segment.length = 0, fontface = 'italic', size = 3.5, bg.color = 'white', bg.r = 0.1) + 
  theme_bw() +
  xlab(log[2]~(Fold-Change~Single-Cell~RNA-seq)) +
  ylab(log[2]~(Fold-Change~BRCA-TCGA~RNA-seq)) +
  theme(plot.title = element_text(face = 2), plot.subtitle = element_text(size = 10))
P1C

attach(df)
o <- predict(lm(bulk~sc2, df), newdata = data.frame(sc=df$sc2), interval = 'prediction', level = 0.95)
df$bulk_sc_lwr <- o[,2]
df$bulk_sc_upr <- o[,3]


corValue <- statsExpressions::corr_test(df,sc2, bulk, type = 'non')$expression[[1]]
df$g <- rownames(df)
df$g[!(((df$bulk < -1) & (df$sc2 < -1)) | ((df$bulk > 1) & (df$sc2 > 1)))] <- NA
df$color <- 'black'
df$color[df$sc2 > 1 & df$bulk > 1] <- 'red'
df$color[df$sc2 < -1 & df$bulk < -1] <- 'blue'
P1D <- ggplot(df, aes(sc, bulk, label = g)) + 
  geom_point(pch = 16, alpha = ifelse(df$color == 'black', 0.25,1), color = df$color) + 
  geom_abline(slope = 1, intercept = 0, color = 'red', lty = 2) + 
  geom_density2d() + 
  labs(title = 'TN Breast Cancer - Healthy Epithelial', subtitle = corValue) +
  geom_text_repel(min.segment.length = 0, fontface = 'italic', size = 3.5, bg.color = 'white', bg.r = 0.1) + 
  theme_bw() +
  xlab(log[2]~(Fold-Change~Single-Cell~RNA-seq)) +
  ylab(log[2]~(Fold-Change~BRCA-TCGA~RNA-seq)) +
  theme(plot.title = element_text(face = 2), plot.subtitle = element_text(size = 10))
P1D


png('../Figures/S2.png', width = 3400, height = 3400, res = 300)
P1B + P1C + P1D +  P1A
dev.off()
