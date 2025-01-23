library(pbapply)
library(fgsea)
library(ggplot2)
library(ggrepel)
library(statsExpressions)
library(patchwork)
library(dplyr)

expressionProfile <- read.csv('../Results/S5_IndividualFC.csv', row.names = 1)
drugProfiles <- read.csv('../Results/S1_Profiles.csv', row.names = 1)

iGenes <- intersect(rownames(expressionProfile), rownames(drugProfiles))

expressionProfile <- expressionProfile[iGenes, ]
drugProfiles <- drugProfiles[iGenes, ]

dPotential <- lapply(colnames(expressionProfile), function(sampleID){
  drugPotential <- t(pbsapply(seq_len(ncol(drugProfiles)), function(X){
    df <- data.frame(A = expressionProfile[,sampleID], B = drugProfiles[,X])
    cValue <- statsExpressions::corr_test(df,A,B, type = 'nonparametric')
    c(drugPotential = cValue$estimate, lowCI = cValue$conf.low, highCI = cValue$conf.high, P= cValue$p.value)
  }))
  drugPotential <- as.data.frame(drugPotential)
  rownames(drugPotential) <- colnames(drugProfiles)
  drugPotential$FDR <- p.adjust(drugPotential$P, method = 'fdr')
  drugPotential <- drugPotential[order(drugPotential$drugPotential),]
  drugPotential <- data.frame(sampleID = sampleID, rank = seq_len(nrow(drugPotential)), compound = rownames(drugPotential), drugPotential)
  return(drugPotential)
})

dPotential <- do.call(rbind.data.frame, dPotential)
dPotential <- dPotential %>% group_by(compound) %>% mutate(compound_avgRank = mean(rank))
write.csv(dPotential, '../Results/S5_drugPotential.csv')

nameCombinations <- t(combn(colnames(drugProfiles),2))

combPotential <- pblapply(colnames(expressionProfile), function(sampleID){
  combinationPotential <- t(apply(nameCombinations,1,function(X){
    df <- data.frame(A = expressionProfile[,sampleID], B = rowMeans(drugProfiles[,X]))
    cValue <- statsExpressions::corr_test(df,A,B, type = 'nonparametric')
    c(combinationPotential = cValue$estimate, lowCI = cValue$conf.low, highCI = cValue$conf.high, P= cValue$p.value)
  }))
  combinationPotential <- data.frame(combinationPotential)
  combinationPotential$FDR <- p.adjust(combinationPotential$P, method = 'fdr')
  rownames(combinationPotential) <- apply(nameCombinations,1,function(X){paste0(X[1],' + ', X[2])})
  combinationPotential <- combinationPotential[order(combinationPotential$combinationPotential),]
  combinationPotential <- data.frame(sampleID = sampleID, rank = seq_len(nrow(combinationPotential)), combination = rownames(combinationPotential), combinationPotential)
  return(combinationPotential)
})

combPotential <- do.call(rbind.data.frame, combPotential)
combPotential <- combPotential %>% group_by(combination) %>% mutate(combination_avgRank = mean(rank))
write.csv(combPotential, '../Results/S5_combinationPotential.csv')


