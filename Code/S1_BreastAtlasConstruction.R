library(ggplot2)
library(fgsea)
library(Seurat)
library(harmony)
library(Matrix)
library(patchwork)
library(Nebulosa)
library(ggrepel)
library(circlize)
library(shadowtext)
hsaPanglaoDB <- gmtPathways('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/markerGenes/hsaPanglaoDB.gmt')


# # Gene dictionary to current
# ENSEMBLid <- read.csv('../Data/ENSEMBL_GeneAlias.txt', sep = '\t')
# geneDictionary <- c(ENSEMBLid$Gene.name, ENSEMBLid$Gene.name)
# names(geneDictionary) <- c(ENSEMBLid$Gene.Synonym, ENSEMBLid$Gene.name)
# 
# # FileNames
# matricesList <- list.files(path = '../Data/UMI/', pattern = '.mtx', full.names = TRUE)
# 
# # Collecting data
# # breastData <- lapply(matricesList, function(fileName){
# #   X <- readMM(fileName)
# #   geneList <- geneDictionary[readLines(gsub('matrix','genes',gsub('mtx', 'txt', fileName)))]
# #   X <- X[!is.na(geneList),]
# #   rownames(X) <- geneList[!is.na(geneList)]
# #   colnames(X) <- readLines(gsub('matrix','barcodes',gsub('mtx', 'txt', fileName)))
# #   sampleName <- gsub('.mtx', '', unlist(strsplit(fileName, '_'))[2])
# #   colnames(X) <- paste0(sampleName, '_C', seq_len(ncol(X)))
# #   X <- CreateSeuratObject(X, project = sampleName)
# #   Idents(X) <- sampleName
# #   X$orig.ident <- Idents(X)
# #   return(X)  
# # })
# readSample <- function(fileName){
#   X <- readMM(fileName)
#   geneList <- geneDictionary[readLines(gsub('matrix','genes',gsub('mtx', 'txt', fileName)))]
#   X <- X[!is.na(geneList),]
#   rownames(X) <- geneList[!is.na(geneList)]
#   colnames(X) <- readLines(gsub('matrix','barcodes',gsub('mtx', 'txt', fileName)))
#   sampleName <- gsub('.mtx', '', unlist(strsplit(fileName, '_'))[2])
#   colnames(X) <- paste0(sampleName, '_C', seq_len(ncol(X)))
#   X <- CreateSeuratObject(X, project = sampleName)
#   Idents(X) <- sampleName
#   X$orig.ident <- Idents(X)
#   return(X)
# }
# 
# # Final Object
# breastData <- readSample(matricesList[1])
# for(i in seq_along(matricesList)[-1]){
#   breastData <- merge(breastData, readSample(matricesList[i]))
# }
# 
# # QC
# source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/scQC.R')
# breastData <- scQC(breastData)
# 
# # Processing
# breastData <- NormalizeData(breastData)
# breastData <- FindVariableFeatures(breastData)
# breastData <- ScaleData(breastData)
# breastData <- RunPCA(breastData)
# breastData <- RunHarmony(breastData, group.by.vars = 'orig.ident', max.iter.harmony = 100)
# 
# # TSNE
# # breastData <- RunTSNE(breastData, reduction = 'harmony', dims = 1:20)
# # TSNEPlot(breastData)
# 
# # UMAP
# breastData <- RunUMAP(breastData, reduction = 'harmony', dims = 1:20)
# UMAPPlot(breastData)

# Loading pre-build Seurat object
load('../Data/breastData.RData')

# Clusters
set.seed(1)
breastData <- FindNeighbors(breastData, reduction = 'umap', dims = 1:2)
set.seed(1)
breastData <- FindClusters(breastData, resolution = 0.005)

# Cell Type Assignation
# cellTypeAssignation <- lapply((seq_along(unique(Idents(breastData))) - 1), function(cID){
#   selectedCells <- Idents(breastData) %in% cID
#   otherAvg <- rowMeans(breastData@assays$RNA@counts[,!selectedCells] != 0)
#   cellAvg <- rowMeans(breastData@assays$RNA@counts[,selectedCells] != 0)
#   log2Freq <- log2(cellAvg/otherAvg)
#   log2Freq <- log2Freq[is.finite(log2Freq)]
#   E <- fgseaMultilevel(hsaPanglaoDB, log2Freq, eps = 0)
#   E <- E[E$NES > 0,]
#   E <- E[order(E$padj),]
# })

levels(Idents(breastData)) <- c('T cells', 
  'Epithelial cells', 
  'Fibroblasts',
  'Smooth muscle cells',
  'Macrophages',
  'B cells',
  'Gamma delta T cells',
  'Endothelial cells',
  'Myoepithelial cells',
  'Memory\nB cells', 
  'Macrophages'
  )

# Cell-types plot
cellTypesPlot <- UMAPPlot(breastData, label = TRUE) + 
  theme_bw() +
  xlim(-15,15) + 
  labs(tag = 'A', title = 'Breast', subtitle = expression(italic(n)==77384~Cells)) + 
  theme(legend.position = 'none', plot.title = element_text(face = 2)) +
  xlab('UMAP 1') +
  ylab('UMAP 2')

# Disease Status
diseaseStatus <- as.vector(breastData$orig.ident)
diseaseStatus[grepl('CID', diseaseStatus)] <- 'Cancer'
diseaseStatus[grepl('sc', diseaseStatus)] <- 'Cancer'
diseaseStatus[grepl('GSM', diseaseStatus)] <- 'Healthy'
diseaseStatus[grepl('TS', diseaseStatus)] <- 'Healthy'
breastData$diseaseStatus <- diseaseStatus

cancerPlot <- UMAPPlot(breastData[,diseaseStatus == 'Cancer'])
cancerPlot <- cancerPlot + 
  labs(title = 'Cancer Samples', subtitle = expression(italic(n)==46594~Cells)) + 
  theme_bw() +
  theme(legend.position = 'none', plot.title = element_text(face = 2)) +
  xlab('UMAP 1') +
  ylab('UMAP 2')

healthyPlot <- UMAPPlot(breastData[,diseaseStatus == 'Healthy'])
healthyPlot <- healthyPlot + 
  labs(tag = 'B', title = 'Healthy Samples', subtitle = expression(italic(n)==30790~Cells)) + 
  theme_bw() +
  theme(legend.position = 'none', plot.title = element_text(face = 2)) +
  xlab('UMAP 1') +
  ylab('UMAP 2')

# Hormone receptors plot
densityPlot <- plot_density(breastData, c('ESR1', 'PGR', 'ERBB2'), joint = TRUE, size = 0.05)
densityPlot <- densityPlot[[4]]
densityPlot <- densityPlot + 
  theme_bw() +
  labs(subtitle = expression(italic(n)==8938~Cells)) +
  theme(legend.position = 'none', plot.title = element_text(face = 2))

# MarkersPlot
# markerGenes <- unlist(hsaPanglaoDB[c('T cells', 
#                       'Epithelial cells', 
#                       'Fibroblasts',
#                       'Smooth muscle cells',
#                       'Macrophages',
#                       'B cells',
#                       'Gamma delta T cells',
#                       'Endothelial cells',
#                       'Myoepithelial cells',
#                       'Memory B cells', 
#                       'Macrophages')])
# markerGenes <- table(markerGenes)
# markerGenes <- names(markerGenes[markerGenes == 1])
# markerGenes <- markerGenes[markerGenes %in% rownames(breastData)]
# 
# markerDE <- FindAllMarkers(breastData, features = markerGenes)
# write.csv(markerDE, file = '../Results/cellTypeMarkerDE.csv')
markerDE <- read.csv('../Results/cellTypeMarkerDE.csv')
markerDE <- markerDE[markerDE$p_val_adj < 0.05,]
markerDE <- lapply(unique(markerDE$cluster),function(X){
  tDExpression <- markerDE[markerDE$cluster %in% X,]
  tDExpression <- tDExpression[order(tDExpression$avg_log2FC, decreasing = TRUE),]
  FC <- tDExpression$avg_log2FC
  names(FC) <- tDExpression$gene
  return(FC[1:5])
})
gLabels <- unique(unlist(lapply(markerDE, names)))

levels(Idents(breastData)) <- c('T cells', 
                                'Epithelial cells', 
                                'Fibroblasts',
                                'Smooth muscle cells',
                                'Macrophages',
                                'B cells',
                                'Gamma delta T cells',
                                'Endothelial cells',
                                'Myoepithelial cells',
                                'Memory B cells', 
                                'Macrophages')

dotPlot <- DotPlot(breastData, features = gLabels, dot.min = 0.20, scale = TRUE, dot.scale = 5) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab('Genes') +
  ylab(NULL) +
  labs(tag = 'C') +
  scale_y_discrete(position = "right") +
  theme(legend.position="bottom", axis.text.x = element_text(face = "italic"))
  
dotPlot


# Plot
pLayout <- '
AAAAAA
AAAAAA
AAAAAA
AAAAAA
AAAAAA
AAAAAA
AAAAAA
AAAAAA
AAAAAA
AAAAAA
BBCCDD
BBCCDD
BBCCDD
BBCCDD
EEEEEE
EEEEEE
EEEEEE
'

png('../Figures/F1.png', width = 4800*.6, height = 4800*.6, res = 300)
cellTypesPlot + healthyPlot + cancerPlot + densityPlot + dotPlot + 
  plot_layout(design = pLayout)
dev.off()

# TNBC Single cell Transcriptional Signatures
breastData <- breastData[,Idents(breastData) %in% 'Epithelial cells']
breastData <- FindNeighbors(breastData, reduction = 'umap', dims = 1:2)
breastData <- FindClusters(breastData, resolution = 0.01)
Idents(breastData) <- ifelse(Idents(breastData) == 0, yes = 'TP', no = 'TN')
Idents(breastData) <- paste0(breastData$diseaseStatus, '_', Idents(breastData))
UMAPPlot(breastData)
DE <- FindMarkers(breastData, ident.1 = 'Cancer_TN', ident.2 = 'Healthy_TP', logfc.threshold = 0)
write.csv(DE, '../Data/de_EC_TNBC-H.csv')
save(breastData, file = '../Data/EpithelialCells.RData')

# Comparing to healthy negative hormone receptors
oDE <- FindMarkers(breastData, ident.1 = 'Cancer_TN', ident.2 = 'Healthy_TN', logfc.threshold = 0)
write.csv(oDE, '../Data/de_EC_TNBC-TNH.csv')

# Comparing H vs C
Idents(breastData) <- breastData$diseaseStatus
oDE <- FindMarkers(breastData, ident.1 = 'Cancer', ident.2 = 'Healthy', logfc.threshold = 0)
write.csv(oDE, '../Data/de_EC_C-H.csv')

# table(Idents(breastData))
# Cancer_TN  Cancer_TP Healthy_TN Healthy_TP 
# 2998       2732       6117       6206 