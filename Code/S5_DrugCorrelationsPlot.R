source('S4_DrugCorrelations.R')

# Improved Performance
(combinationPotential[1,1]/drugPotential[1,1])-1
# 0.1114373

# Pharmacological effects
drugProfiles <- read.csv('../Results/S1_Profiles.csv', row.names = 1)
drugProfile <- drugProfiles[,rownames(drugPotential)[1]]
names(drugProfile) <- rownames(drugProfiles)

# Plot 
df <- data.frame(SC = expressionProfile[iGenes,2], DP = drugProfile[iGenes])
eLabel <- corr_test(df, SC,DP, type = 'nonp')$expression[[1]]
df$G <- rownames(df)
fcLimit <- 0.25
df$C <- 'black'
df$G[(df$SC > -fcLimit & df$DP > -fcLimit)] <- NA
df$C[(df$SC < -fcLimit & df$DP > -fcLimit)] <- 'red'
df$G[(df$SC < fcLimit & df$DP < fcLimit)] <- NA
df$C[(df$SC > -fcLimit & df$DP < -fcLimit)] <- 'blue'
df$C[is.na(df$G)] <- 'black'
df$alpha <- 1
df$alpha[is.na(df$G)] <- 0.25

F3A <- ggplot(df, aes(SC, DP, label = G)) + 
  geom_abline(slope = -1, intercept = 0, lty = 2, color = 'red') +
  geom_point(pch = 16, alpha = df$alpha, color = df$C) + 
  geom_density_2d() + 
  geom_text_repel(min.segment.length = 0, fontface = 3, size = 4, bg.color = 'white') +
  theme_bw() +
  labs(tag = 'A', title = 'QL-XII-47', subtitle = eLabel) +
  theme(plot.title = element_text(face = 2)) +
  xlab(expression(log[2]~(Fold-Change~Single-Cell~RNA-seq))) +
  ylab(expression(Effect~Sizes~LINC~L1000~Project)) +
  theme(plot.subtitle = element_text())
F3A

MSigDB_Hallmarks <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020')
KEGG <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human')
BIOP <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')

set.seed(1)
#E <- fgseaMultilevel(BIOP, drugProfile, eps = 0)
E <- fgseaMultilevel(MSigDB_Hallmarks, drugProfile, eps = 0)
E <- E[E$padj < 0.05,]
E <- E[order(E$padj),]

eRank <- 1
F3B1 <- plotEnrichment(MSigDB_Hallmarks[[E$pathway[eRank]]], drugProfile) +
  xlab('Gene Rank') +
  ylab('Enrichment Score') +
  labs(tag = 'B', title = E$pathway[eRank], 
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

drugCombination <- rowMeans(drugProfiles[,c('GSK.690693', 'QL.XII.47')])
df <- data.frame(SC = expressionProfile[iGenes,2], DP = drugCombination[iGenes])
eLabel <- corr_test(df, SC,DP, type = 'nonp')$expression[[1]]
df$G <- rownames(df)
fcLimit <- 0.25
df$C <- 'black'
df$G[(df$SC > -fcLimit & df$DP > -fcLimit)] <- NA
df$C[(df$SC < -fcLimit & df$DP > -fcLimit)] <- 'red'
df$G[(df$SC < fcLimit & df$DP < fcLimit)] <- NA
df$C[(df$SC > -fcLimit & df$DP < -fcLimit)] <- 'blue'
df$C[is.na(df$G)] <- 'black'
df$alpha <- 1
df$alpha[is.na(df$G)] <- 0.25

F3C <- ggplot(df, aes(SC, DP, label = G)) + 
  geom_abline(slope = -1, intercept = 0, lty = 2, color = 'red') +
  geom_point(pch = 16, alpha = df$alpha, color = df$C) + 
  geom_density2d() + 
  geom_text_repel(min.segment.length = 0, fontface = 3, size = 4, bg.color = 'white') +
  theme_bw() +
  labs(tag = 'D', title = 'QL-XII-47 + GSK-690693', subtitle = eLabel) +
  theme(plot.title = element_text(face = 2)) +
  xlab(expression(log[2]~(Fold-Change~Single-Cell~RNA-seq))) +
  ylab(expression(Effect~Sizes~LINC~L1000~Combinations)) +
  theme(plot.subtitle = element_text())
F3C

set.seed(1)
E <- fgseaMultilevel(MSigDB_Hallmarks, drugCombination)
E <- E[E$padj < 0.05,]
E <- E[order(E$padj),]

eRank <- 1
F3D1 <- plotEnrichment(MSigDB_Hallmarks[[E$pathway[eRank]]], drugCombination) +
  xlab('Gene Rank') +
  ylab('Enrichment Score') +
  labs(tag = 'E', title = E$pathway[eRank], 
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

eRank <- 7
F3D6 <- plotEnrichment(MSigDB_Hallmarks[[E$pathway[eRank]]], drugCombination) +
  xlab('Gene Rank') +
  ylab('Enrichment Score') +
  labs(title = E$pathway[eRank], 
       subtitle = parse(text = paste0('NES == ',round(E$NES[eRank],2),'~~P-adj==',formatC(E$padj[eRank],digits = 2, format = 'g')))) +
  theme_bw() +
  theme(plot.title = element_text(face = 2))


EXP <- read.csv('../Data/dataset_20367_20210819011500.csv')
EXP <- EXP[EXP$Drug.name.Name %in% c('QL-XII-47'),]
F3A1 <- ggplot(EXP, aes(EXP$GRmax, EXP$Cell.line.Name)) + 
  geom_boxplot(fill = NA) + 
  geom_jitter(height = 0) + 
  theme_bw() + 
  xlim(c(-1,1)) +
  xlab(parse(text = 'GR[max]~value' )) + ylab('TNBC Cell Lines') +
  labs(tag = 'C', title = 'QL-XII-47 Sensitivity', subtitle = 'TNBC Cell Lines') +
  theme(plot.title = element_text(face = 2)) +
  geom_vline(xintercept = 0, lty = 2, col = 'red')

plotLayout <- '
AABC
AADE
AAFN
GGHI
GGJK
GGLM'

png('../Figures/F4.png', width = 4800 * 0.85, height = 4800 * 0.7, res = 300)
F3A + F3B1 + F3B2 + F3B3 + F3B4 + F3B5 + F3C + F3D1 + F3D2 + F3D3 + F3D4 + F3D5 + F3D6 + F3A1 + plot_layout(design = plotLayout)
dev.off()

