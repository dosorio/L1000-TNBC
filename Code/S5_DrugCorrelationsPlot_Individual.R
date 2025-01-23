library(ggplot2)
library(dplyr)
library(readr)
library(ComplexHeatmap)
library(reshape2)

dPotential <- read.csv('../Results/S5_drugPotential.csv', row.names = 1)

topSingle <- sort(table(dPotential$compound[dPotential$rank == 1]), decreasing = TRUE)

dPotentialMatrix <- dcast(dPotential, compound ~ sampleID, value.var = 'drugPotential') 
rownames(dPotentialMatrix) <- dPotentialMatrix$compound
dPotentialMatrix <- dPotentialMatrix[,-1]
dPotentialMatrix <- dPotentialMatrix[names(topSingle),]
topSingleLabels <- gsub('\\.', '-', paste0(toupper(names(topSingle)), ' (', topSingle, '/', sum(topSingle), ')'))
sampleID <- unique(dPotential$sampleID)
sampleID <- paste0(sampleID, ifelse(sampleID %in% unique(dPotential$sampleID[dPotential$FDR < 0.05]), '*', ''))
sampleID <- gsub('_Cancer_TN', '', sampleID)
dPotentialMatrix <- as.matrix(dPotentialMatrix)

png('../Figures/S5_dPotential.png', width = 2500, height = 750, res = 300)
HM1 <- Heatmap(dPotentialMatrix, name = 'Drug Potential', 
        row_labels = topSingleLabels, 
        column_labels = sampleID, cluster_rows = FALSE, show_column_dend = FALSE)
draw(HM1, heatmap_legend_side = "right")
dev.off()


cPotential <- read.csv('../Results/S5_combinationPotential.csv', row.names = 1)
topCombination <- sort(table(cPotential$combination[cPotential$rank == 1]), decreasing = TRUE)

cPotentialMatrix <- dcast(cPotential, combination ~ sampleID, value.var = 'combinationPotential') 
rownames(cPotentialMatrix) <- cPotentialMatrix$combination
cPotentialMatrix <- cPotentialMatrix[,-1]
cPotentialMatrix <- cPotentialMatrix[names(topCombination),]
topCombinationLabels <- gsub('\\.', '-', paste0(toupper(names(topCombination)), ' (', topCombination, '/', sum(topCombination), ')'))
sampleID <- unique(cPotential$sampleID)
sampleID <- paste0(sampleID, ifelse(sampleID %in% unique(cPotential$sampleID[cPotential$FDR < 0.05]), '*', ''))
sampleID <- gsub('_Cancer_TN', '', sampleID)
cPotentialMatrix <- as.matrix(cPotentialMatrix)

png('../Figures/S5_cPotential.png', width = 2500, height = 1000, res = 300)
HM2 <- Heatmap(cPotentialMatrix, name = 'Combination Potential', 
               row_labels = topCombinationLabels)
draw(HM2, heatmap_legend_side = "right")
dev.off()

library(circlize)
col_fun = colorRamp2(c(-0.3, 0, 0.3), c("red", "white", "blue"))
rowSplit <- factor(c(rep('Compound', length(topSingleLabels)), rep('Combination', length(topCombinationLabels))),
                   c('Compound', 'Combination'))
fontType <- rep(1, length(rowSplit))
fontType[1] <- 2
fontType[length(topSingleLabels)+1] <- 2
png('../Figures/F5.png', width = 3000, height = 850, res = 300)
HM3 <- Heatmap(rbind(dPotentialMatrix, cPotentialMatrix), 
               name = 'Potential',
               col = col_fun,
               row_labels = c(topSingleLabels, topCombinationLabels),
               column_labels = sampleID, 
               cluster_rows = FALSE, 
               show_column_dend = FALSE,
               row_split = rowSplit,
               show_column_names = FALSE,
               row_gap = unit(5, "mm"),
               heatmap_legend_param = list(title = expression(rho), direction = "vertical")) +
  rowAnnotation(link = anno_mark(at = seq_along(c(topSingleLabels, topCombinationLabels)), 
                                 labels =c(topSingleLabels, topCombinationLabels), 
                                 labels_gp = gpar(fontsize = 8, font = fontType), padding = unit(0.5, "mm")))
draw(HM3, heatmap_legend_side = "right")
dev.off()
