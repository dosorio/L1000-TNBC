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

# Bulk
load('~/BRCA.RData')
selectedSamples <- data@colData$paper_BRCA_Subtype_PAM50 %in% c('Normal', 'Basal')
eData <- data@assays@data$`HTSeq - FPKM-UQ`[,selectedSamples]
eData <- as.matrix(eData)
eData <- apply(eData,2,as.integer)
rownames(eData) <- data@rowRanges$external_gene_name
colnames(eData) <- data@colData$barcode[selectedSamples]
eData <- eData[complete.cases(eData),]
eData <- eData[rownames(eData) %in% intersect(names(scFC), names(pcFC)),]
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

# Label Selection
g1 <- DE$X[((abs(DE$avg_log2FC) > 1) & (DE$avg_log2FC < 0.05))]
g2 <- rownames(tData)[((abs(tData$log2FoldChange) > 1) & (tData$padj < 0.05))]
g3 <- rownames(oFit)[((abs(oFit$log2FoldChange) > 1) & (oFit$padj < 0.05))]

sGenes <- c(g1,g2,g3)
sGenes <- table(sGenes)
sGenes <- names(sGenes[sGenes >= 2])

# Plots
iGenes <- intersect(intersect(names(pcFC),names(scFC)), names(bulkFC))

DE <- DE[DE$X %in% iGenes,]
tData <- tData[rownames(tData) %in% iGenes,]
oFit <- oFit[rownames(oFit) %in% iGenes,]

tData$G <- rownames(tData)
DE$G <- DE$X
oFit$G <- rownames(oFit)

tData$G[!tData$G %in% sGenes] <- NA
DE$G[!DE$G %in% sGenes] <- NA
oFit$G[!oFit$G %in% sGenes] <- NA

DE$G[DE$color == 'black'] <- NA
tData$G[tData$color == 'black'] <- NA
oFit$G[oFit$color == 'black'] <- NA

P1A <- ggplot(DE, aes(avg_log2FC, -log10(p_val_adj), label = G)) + 
  geom_point(pch = 16, alpha = ifelse(DE$color == 'black', 0.25,1), color = DE$color) + 
  theme_bw() +
  geom_text_repel(min.segment.length = 0, fontface = 'italic', size = 3.5, max.overlaps = 5, force = 5, segment.color = 'gray60', bg.color = 'white', bg.r = 0.1) +
  xlab(log[2]~(Fold-Change~Single-Cell~RNA-seq)) +
  ylab(-log[10]~(P-value)) +
  labs(tag = 'A', title = 'Single-cell RNA-seq', subtitle = '2998 TNBC vs. 6206 Healthy Epithelial Cells') +
  theme(plot.title = element_text(face = 2), plot.subtitle = element_text(size = 10))
P1A

P1B <- ggplot(tData, aes(log2FoldChange, -log10(pvalue), label = G)) + 
  geom_point(pch = 16, alpha = ifelse(tData$color == 'black', 0.25,1), color = tData$color) + 
  xlim(-10,10) + 
  ylim(0, 25) +
  theme_bw() + 
  geom_text_repel(min.segment.length = 0, fontface = 'italic', size = 3.5, max.overlaps = 10, force = 5, segment.color = 'gray60', bg.color = 'white', bg.r = 0.1) +
  xlab(log[2]~(Fold-Change~PseudoCounts)) +
  ylab(-log[10]~(P-value)) +
  labs(title = 'Pseudobulk', subtitle = '24 TNBC vs. 10 Healthy Samples') +
  theme(plot.title = element_text(face = 2), plot.subtitle = element_text(size = 10))
P1B

P1C <- ggplot(oFit, aes(log2FoldChange, -log10(pvalue), label = G)) + 
  geom_point(pch = 16, alpha = ifelse(oFit$color == 'black', 0.25, 1), color = oFit$color) + 
  #xlim(-6,6) + 
  #ylim(0, 40) +
  theme_bw() + 
  geom_text_repel(min.segment.length = 0, fontface = 'italic', size = 3.5, max.overlaps = 10, force = 5, segment.color = 'gray60', bg.color = 'white', bg.r = 0.1) +
  xlab(log[2]~(Fold-Change~Bulk~RNA-seq)) +
  ylab(-log[10]~(P-value)) +
  labs(title = 'Bulk BRCA-TCGA RNA-seq', subtitle = '194 TNBC vs. 40 Healthy Samples') +
  theme(plot.title = element_text(face = 2), plot.subtitle = element_text(size = 10))
P1C

df <- data.frame(sc = scFC[iGenes], pc = pcFC[iGenes], bulk = bulkFC[iGenes])
df <- df[complete.cases(df),]

# Heatmap
MSigDB_Hallmarks <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020')
ssgsea <- gsva(as.matrix(df), MSigDB_Hallmarks, method='ssgsea')
ssgsea <- round(ssgsea,2)
colnames(ssgsea) <- c('Single-Cell', 'Pseudo-Bulk', 'Bulk')
rownames(ssgsea) <- gsub('Pathway','P.',rownames(ssgsea))
sPath <- rowSums(ssgsea < 0) %in% c(0,3)
HM <- Heatmap(ssgsea, show_row_names = TRUE, name = 'ES', column_order = 1:3, show_row_dend = FALSE, show_column_dend = FALSE,
              heatmap_legend_param = list(at = c(-1, 0, 1), grid_width = unit(2, "mm"))) +
  rowAnnotation(link = anno_mark(at = which(sPath),
                                 labels = rownames(ssgsea)[sPath]))
HM <- ggplotify::as.ggplot(grid::grid.grabExpr(ComplexHeatmap::draw(HM)))
HM <- HM + labs(tag = 'C')
HM

o <- predict(lm(pc~sc, df), newdata = data.frame(sc=df$sc), interval = 'prediction', level = 0.95)
df$pc_sc_lwr <- o[,2]
df$pc_sc_upr <- o[,3]

o <- predict(lm(bulk~sc, df), newdata = data.frame(sc=df$sc), interval = 'prediction', level = 0.95)
df$bulk_sc_lwr <- o[,2]
df$bulk_sc_upr <- o[,3]

corValue <- statsExpressions::corr_test(df,sc, pc, type = 'non')$expression[[1]]
df$g <- rownames(df)
df$g[!(((df$pc < -1) & (df$sc < -1)) | ((df$pc > 1) & (df$sc > 1)))] <- NA
df$color <- 'black'
df$color[df$sc > 1 & df$pc > 1] <- 'red'
df$color[df$sc < -1 & df$pc < -1] <- 'blue'
P1D <- ggplot(df, aes(sc, pc, label = g)) + 
  geom_point(pch = 16, alpha = ifelse(df$color == 'black', 0.25,1), color = df$color) +
  geom_abline(slope = 1, intercept = 0, color = 'red', lty = 2) + 
  geom_density2d() + 
  #geom_smooth(method = 'lm', se = FALSE, col = 'red') +
  #geom_line(aes(sc, pc_sc_lwr), data = df, col = 'red', lty = 2) +
  #geom_line(aes(sc, pc_sc_upr), data = df, col = 'red', lty = 2) +
  labs(tag = 'B', title = 'Pseudobulk vs. Single-Cell', subtitle = corValue) +
  ylim(-10,10) +
  geom_text_repel(min.segment.length = 0, fontface = 'italic', size = 3.5, bg.color = 'white', bg.r = 0.1) + 
  theme_bw() +
  xlab(log[2]~(Fold-Change~Single-Cell~RNA-seq)) +
  ylab(log[2]~(Fold-Change~PseudoCounts)) +
  theme(plot.title = element_text(face = 2), plot.subtitle = element_text(size = 10))
P1D

corValue <- statsExpressions::corr_test(df,sc, bulk, type = 'non')$expression[[1]]
df$g <- rownames(df)
df$g[!(((df$bulk < -1) & (df$sc < -1)) | ((df$bulk > 1) & (df$sc > 1)))] <- NA
df$color <- 'black'
df$color[df$sc > 1 & df$bulk > 1] <- 'red'
df$color[df$sc < -1 & df$bulk < -1] <- 'blue'
P1E <- ggplot(df, aes(sc, bulk, label = g)) + 
  geom_point(pch = 16, alpha = ifelse(df$color == 'black', 0.25,1), color = df$color) + 
  geom_abline(slope = 1, intercept = 0, color = 'red', lty = 2) + 
  geom_density2d() + 
  labs(title = 'Bulk vs. Single-Cell', subtitle = corValue) +
  geom_text_repel(min.segment.length = 0, fontface = 'italic', size = 3.5, bg.color = 'white', bg.r = 0.1) + 
  theme_bw() +
  xlab(log[2]~(Fold-Change~Single-Cell~RNA-seq)) +
  ylab(log[2]~(Fold-Change~BRCA-TCGA~RNA-seq)) +
  theme(plot.title = element_text(face = 2), plot.subtitle = element_text(size = 10))
P1E 

# Personalized medicine
TN <- breastData[,Idents(breastData) %in% c('Cancer_TN', 'Healthy_TP')]
S <- unique(data.frame(TN$orig.ident, TN$diseaseStatus))
E <- sapply(S$TN.orig.ident, function(OI){
  rowSums(TN@assays$RNA@data[,TN$orig.ident %in% OI])
})
E <- gsva(E, MSigDB_Hallmarks)
PS <- ComplexHeatmap::Heatmap(E, column_title_gp = gpar(fill = c("#F24405", "#348888", "#22BABB", "#FA7F08", "#9EF8EE"), font = 2, color = 'white'),
                              column_split = S$Subtype, 
                              name = 'ES', 
                              show_row_dend = FALSE, 
                              show_column_names = TRUE, 
                              heatmap_legend_param = list(at = c(-1, 0, 1), grid_width = unit(2, "mm"))) +
  rowAnnotation(link = anno_mark(at = which(sPath),
                                 labels = rownames(E)[sPath]))
PS
#PS <- ggplotify::as.ggplot(grid::grid.grabExpr(ComplexHeatmap::draw(PS)))

dP <- '
AABBCCFFF
AABBCCFFF
DDDEEEFFF
DDDEEEFFF'

png('../Figures/F2.png', width = 4650, height = 3000, res = 300)
P1A + P1B + P1C + P1D + P1E + HM +  plot_layout(design = dP)
dev.off()

png('../Figures/F5.png', width = 4000, height = 2000, res = 300)
PS
dev.off()




### Levels of expression
load('../Data/EpithelialCells.RData')
Idents(breastData) <- gsub('_', ' ', Idents(breastData))
Idents(breastData) <- factor(Idents(breastData), levels = c('Cancer TP', 'Healthy TP','Cancer TN', 'Healthy TN'))
PS2 <- DotPlot(breastData, features = c('ESR1', 'PGR', 'ERBB2'), dot.scale = 10, cols = c('white', 'darkblue'), col.min = -0.2, col.max = 1)
PS2 <- PS2 + theme_bw() + xlab('Genes') +
  theme(panel.grid = element_blank()) + ylab('Epithelial Cell')

png('../Figures/S4.png', width = 1250, height = 1000, res = 300)
print(PS2)
dev.off()
