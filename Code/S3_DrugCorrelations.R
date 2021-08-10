library(pbapply)
library(fgsea)
library(ggplot2)
library(ggrepel)

expressionProfile <- read.csv('../Data/de_EC_TNBC-H.csv', row.names = 1)
drugProfiles <- read.csv('../Results/S1_Profiles.csv', row.names = 1)

iGenes <- intersect(rownames(expressionProfile), rownames(drugProfiles))

expressionProfile <- expressionProfile[iGenes, ]
drugProfiles <- drugProfiles[iGenes, ]

drugPotential <- data.frame(TNBC = expressionProfile$avg_log2FC, drugProfiles)
drugPotential <- sort(cor(drugPotential, method = 'sp')[,1])
drugPotential <- data.frame(drugPotential)
write.csv(drugPotential, '../Results/S3_drugPotential.csv')

nameCombinations <- t(combn(colnames(drugProfiles),2))
combinationPotential <- pbapply(nameCombinations,1,function(X){
  cor(expressionProfile$avg_log2FC, rowMeans(drugProfiles[,X]), method = 'sp')
})
names(combinationPotential) <- apply(nameCombinations,1,function(X){paste0(X[1],' + ', X[2])})
combinationPotential <- sort(combinationPotential)
combinationPotential <- data.frame(combinationPotential)
write.csv(combinationPotential, '../Results/S3_combinationPotential.csv')


# Improved Performance
(combinationPotential[1,1]/drugPotential[1,1])-1



# Pharmacological effects
drugProfiles <- read.csv('../Results/S1_Profiles.csv', row.names = 1)
drugProfile <- drugProfiles[,rownames(drugPotential)[1]]
names(drugProfile) <- rownames(drugProfiles)

# Plot 
df <- data.frame(SC = expressionProfile[iGenes,2], DP = drugProfile[iGenes])
df <- df[abs(df$SC) > 1,]
df$G <- rownames(df)
df$G[!((abs(df$SC) > 0.6) & (abs(df$DP) > 0.6))] <- NA
ggplot(df, aes(SC, DP, label = G)) + 
  geom_smooth(method = 'lm', color = 'red', alpha = 0.25) + 
  geom_point() + 
  geom_text_repel(min.segment.length = 0, fontface = 3) +
  theme_bw() +
  labs(title = 'QL-XII-47') +
  theme(plot.title = element_text(face = 2))

Hallmarks <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020')
E <- fgseaMultilevel(Hallmarks, drugProfile)
E <- E[E$padj < 0.05,]


