library(pbapply)
library(fgsea)
library(ggplot2)
library(ggrepel)
library(statsExpressions)
library(patchwork)

expressionProfile <- read.csv('../Data/de_EC_TNBC-H.csv', row.names = 1)
drugProfiles <- read.csv('../Results/S1_Profiles.csv', row.names = 1)

iGenes <- intersect(rownames(expressionProfile), rownames(drugProfiles))

expressionProfile <- expressionProfile[iGenes, ]
drugProfiles <- drugProfiles[iGenes, ]

drugPotential <- t(pbsapply(seq_len(ncol(drugProfiles)), function(X){
  df <- data.frame(A = expressionProfile$avg_log2FC, B = drugProfiles[,X])
  cValue <- statsExpressions::corr_test(df,A,B, type = 'nonparametric')
  c(drugPotential = cValue$estimate, lowCI = cValue$conf.low, highCI = cValue$conf.high, P= cValue$p.value)
}))
drugPotential <- as.data.frame(drugPotential)
rownames(drugPotential) <- colnames(drugProfiles)
drugPotential$FDR <- p.adjust(drugPotential$P, method = 'fdr')
drugPotential <- drugPotential[order(drugPotential$drugPotential),]
write.csv(drugPotential, '../Results/S3_drugPotential.csv')

nameCombinations <- t(combn(colnames(drugProfiles),2))

combinationPotential <- t(pbapply(nameCombinations,1,function(X){
  df <- data.frame(A = expressionProfile$avg_log2FC, B = rowMeans(drugProfiles[,X]))
  cValue <- statsExpressions::corr_test(df,A,B, type = 'nonparametric')
  c(combinationPotential = cValue$estimate, lowCI = cValue$conf.low, highCI = cValue$conf.high, P= cValue$p.value)
}))
combinationPotential <- data.frame(combinationPotential)
combinationPotential$FDR <- p.adjust(combinationPotential$P, method = 'fdr')
rownames(combinationPotential) <- apply(nameCombinations,1,function(X){paste0(X[1],' + ', X[2])})
combinationPotential <- combinationPotential[order(combinationPotential$combinationPotential),]
write.csv(combinationPotential, '../Results/S3_combinationPotential.csv')


