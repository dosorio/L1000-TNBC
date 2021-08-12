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
eLabel <- corr_test(df, SC,DP, type = 'nonp')$expression[[1]]
df$G <- rownames(df)
fcLimit <- 0.25
df$G[(df$SC > -fcLimit & df$DP > -fcLimit)] <- NA
df$G[(df$SC < fcLimit & df$DP < fcLimit)] <- NA
df$alpha <- 0.5
df$alpha[is.na(df$G)] <- 0.1

F3A <- ggplot(df, aes(SC, DP, label = G)) + 
  geom_smooth(method = 'lm', color = 'red', alpha = 0.25) + 
  geom_point(pch = 16, alpha = df$alpha) + 
  geom_text_repel(min.segment.length = 0, fontface = 3, size = 4) +
  theme_bw() +
  labs(tag = 'A', title = 'QL-XII-47', subtitle = eLabel) +
  theme(plot.title = element_text(face = 2)) +
  xlab(expression(log[2]~(Fold-Change~Single-Cell~RNA-seq))) +
  ylab(expression(Effect~Sizes~LINC~L1000~Project)) +
  theme(plot.subtitle = element_text())

MSigDB_Hallmarks <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020')
set.seed(1)
E <- fgseaMultilevel(MSigDB_Hallmarks, drugProfile, eps = 0)
E <- E[E$padj < 0.05,]
E <- E[order(E$padj),]

eRank <- 1
F3B1 <- plotEnrichment(MSigDB_Hallmarks[[E$pathway[eRank]]], drugProfile) +
  xlab('Gene Rank') +
  ylab('Enrichment Score') +
  labs(title = E$pathway[eRank], 
       subtitle = parse(text = paste0('NES == ',round(E$NES[eRank],2),'~~P-adj==',formatC(E$padj[eRank],digits = 2, format = 'g')))) +
  theme_bw() +
  theme(plot.title = element_text(face = 2))

eRank <- 2
F3B2 <- plotEnrichment(MSigDB_Hallmarks[[E$pathway[eRank]]], drugProfile) +
  xlab('Gene Rank') +
  ylab('Enrichment Score') +
  labs(title = E$pathway[eRank], 
       subtitle = parse(text = paste0('NES == ',round(E$NES[eRank],2),'~~P-adj==',formatC(E$padj[eRank],digits = 2, format = 'g')))) +
  theme_bw() +
  theme(plot.title = element_text(face = 2))

eRank <- 3
F3B3 <- plotEnrichment(MSigDB_Hallmarks[[E$pathway[eRank]]], drugProfile) +
  xlab('Gene Rank') +
  ylab('Enrichment Score') +
  labs(title = E$pathway[eRank], 
       subtitle = parse(text = paste0('NES == ',round(E$NES[eRank],2),'~~P-adj==',formatC(E$padj[eRank],digits = 2, format = 'g')))) +
  theme_bw() +
  theme(plot.title = element_text(face = 2))

eRank <- 4
F3B4 <- plotEnrichment(MSigDB_Hallmarks[[E$pathway[eRank]]], drugProfile) +
  xlab('Gene Rank') +
  ylab('Enrichment Score') +
  labs(title = E$pathway[eRank], 
       subtitle = parse(text = paste0('NES == ',round(E$NES[eRank],2),'~~P-adj==',formatC(E$padj[eRank],digits = 2, format = 'g')))) +
  theme_bw() +
  theme(plot.title = element_text(face = 2))

eRank <- 5
F3B5 <- plotEnrichment(MSigDB_Hallmarks[[E$pathway[eRank]]], drugProfile) +
  xlab('Gene Rank') +
  ylab('Enrichment Score') +
  labs(title = E$pathway[eRank], 
       subtitle = parse(text = paste0('NES == ',round(E$NES[eRank],2),'~~P-adj==',formatC(E$padj[eRank],digits = 2, format = 'g')))) +
  theme_bw() +
  theme(plot.title = element_text(face = 2))

drugCombination <- rowMeans(drugProfiles[,c('nilotinib', 'QL.XII.47')])
df <- data.frame(SC = expressionProfile[iGenes,2], DP = drugCombination[iGenes])
eLabel <- corr_test(df, SC,DP, type = 'nonp')$expression[[1]]
df$G <- rownames(df)
fcLimit <- 0.25
df$G[(df$SC > -fcLimit & df$DP > -fcLimit)] <- NA
df$G[(df$SC < fcLimit & df$DP < fcLimit)] <- NA
df$alpha <- 0.5
df$alpha[is.na(df$G)] <- 0.1

F3C <- ggplot(df, aes(SC, DP, label = G)) + 
  geom_smooth(method = 'lm', color = 'red', alpha = 0.25) + 
  geom_point(pch = 16, alpha = df$alpha) + 
  geom_text_repel(min.segment.length = 0, fontface = 3, size = 4) +
  theme_bw() +
  labs(tag = 'B', title = 'QL-XII-47 + Nilotinib', subtitle = eLabel) +
  theme(plot.title = element_text(face = 2)) +
  xlab(expression(log[2]~(Fold-Change~Single-Cell~RNA-seq))) +
  ylab(expression(Effect~Sizes~LINC~L1000~Combinations)) +
  theme(plot.subtitle = element_text())

set.seed(1)
E <- fgseaMultilevel(MSigDB_Hallmarks, drugCombination)
E <- E[E$padj < 0.05,]
E <- E[order(E$padj),]

eRank <- 1
F3D1 <- plotEnrichment(MSigDB_Hallmarks[[E$pathway[eRank]]], drugCombination) +
  xlab('Gene Rank') +
  ylab('Enrichment Score') +
  labs(title = E$pathway[eRank], 
       subtitle = parse(text = paste0('NES == ',round(E$NES[eRank],2),'~~P-adj==',formatC(E$padj[eRank],digits = 2, format = 'g')))) +
  theme_bw() +
  theme(plot.title = element_text(face = 2))

eRank <- 2
F3D2 <- plotEnrichment(MSigDB_Hallmarks[[E$pathway[eRank]]], drugCombination) +
  xlab('Gene Rank') +
  ylab('Enrichment Score') +
  labs(title = E$pathway[eRank], 
       subtitle = parse(text = paste0('NES == ',round(E$NES[eRank],2),'~~P-adj==',formatC(E$padj[eRank],digits = 2, format = 'g')))) +
  theme_bw() +
  theme(plot.title = element_text(face = 2))

eRank <- 3
F3D3 <- plotEnrichment(MSigDB_Hallmarks[[E$pathway[eRank]]], drugCombination) +
  xlab('Gene Rank') +
  ylab('Enrichment Score') +
  labs(title = E$pathway[eRank], 
       subtitle = parse(text = paste0('NES == ',round(E$NES[eRank],2),'~~P-adj==',formatC(E$padj[eRank],digits = 2, format = 'g')))) +
  theme_bw() +
  theme(plot.title = element_text(face = 2))

eRank <- 4
F3D4 <- plotEnrichment(MSigDB_Hallmarks[[E$pathway[eRank]]], drugCombination) +
  xlab('Gene Rank') +
  ylab('Enrichment Score') +
  labs(title = E$pathway[eRank], 
       subtitle = parse(text = paste0('NES == ',round(E$NES[eRank],2),'~~P-adj==',formatC(E$padj[eRank],digits = 2, format = 'g')))) +
  theme_bw() +
  theme(plot.title = element_text(face = 2))

eRank <- 5
F3D5 <- plotEnrichment(MSigDB_Hallmarks[[E$pathway[eRank]]], drugCombination) +
  xlab('Gene Rank') +
  ylab('Enrichment Score') +
  labs(title = E$pathway[eRank], 
       subtitle = parse(text = paste0('NES == ',round(E$NES[eRank],2),'~~P-adj==',formatC(E$padj[eRank],digits = 2, format = 'g')))) +
  theme_bw() +
  theme(plot.title = element_text(face = 2))

eRank <- 6
F3D6 <- plotEnrichment(MSigDB_Hallmarks[[E$pathway[eRank]]], drugCombination) +
  xlab('Gene Rank') +
  ylab('Enrichment Score') +
  labs(title = E$pathway[eRank], 
       subtitle = parse(text = paste0('NES == ',round(E$NES[eRank],2),'~~P-adj==',formatC(E$padj[eRank],digits = 2, format = 'g')))) +
  theme_bw() +
  theme(plot.title = element_text(face = 2))

eRank <- 7
F3D7 <- plotEnrichment(MSigDB_Hallmarks[[E$pathway[eRank]]], drugCombination) +
  xlab('Gene Rank') +
  ylab('Enrichment Score') +
  labs(title = E$pathway[eRank], 
       subtitle = parse(text = paste0('NES == ',round(E$NES[eRank],2),'~~P-adj==',formatC(E$padj[eRank],digits = 2, format = 'g')))) +
  theme_bw() +
  theme(plot.title = element_text(face = 2))

eRank <- 8
F3D8 <- plotEnrichment(MSigDB_Hallmarks[[E$pathway[eRank]]], drugCombination) +
  xlab('Gene Rank') +
  ylab('Enrichment Score') +
  labs(title = E$pathway[eRank], 
       subtitle = parse(text = paste0('NES == ',round(E$NES[eRank],2),'~~P-adj==',formatC(E$padj[eRank],digits = 2, format = 'g')))) +
  theme_bw() +
  theme(plot.title = element_text(face = 2))

plotLayout <- '
AAAABC
AAAADE
AAAAF#
AAAA##
GGGGHI
GGGGJK
GGGGLM
GGGGNO'

png('../Figures/F3.png', width = 4800 * 0.9, height = 4800  * 0.9, res = 300)
F3A + F3B1 + F3B2 + F3B3 + F3B4 + F3B5 + F3C + F3D1 + F3D2 + F3D3 + F3D4 + F3D5 + F3D6 + F3D7 + F3D8 + plot_layout(design = plotLayout)
dev.off()
