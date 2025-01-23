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

library(dplyr)
library(ggpubr)
library(patchwork)
t1Data <- data.frame(A = c(52.81, 54.48, 14.56, NA),  B= c(88.29, 87.45, 92.12, 90.04), C = c(46.07, 13.75, 48.48, 40.83), D = c(100.1,100,99.9, NA))
t1Data <- reshape2::melt(t1Data)
t1Data <- na.omit(t1Data)
levels(t1Data$variable) <- c('QL-XII-47', 'GSK-690693', 'Combination', 'Control')
t1Data$variable <- factor(t1Data$variable, c('Combination', 'QL-XII-47', 'GSK-690693', 'Control') )
t1Data$cellLine <- 'DU4475'
t1Summary <- t1Data %>% group_by(variable) %>% summarize(M = mean(value), L = M - sd(value), U = M + sd(value))

T1 <- ggplot(t1Data, aes(value, variable)) + 
  geom_errorbarh(data = t1Summary, aes(y = variable, x = M, xmin = L, xmax = U), height = 0.2) + 
  geom_bar(position = "dodge", stat = "summary",fun = "mean", fill = c('gray30','gray50','gray70','gray90')) +
  geom_jitter(height = 0.1, alpha = 0.5, pch = 16) + theme_light() +
  ylab('Treatment') +
  xlab('Viability') +
  stat_compare_means(aes(label = ..p.signif..), 
                     comparisons = list(c(1,2), c(1,3), c(1,4)), 
                     method = 't.test', 
                     method.args = list(alternative = 'less'), 
                     hide.ns = TRUE, tip.length = 0.01) +
  theme(plot.title = element_text(face = 2)) +
  labs(title = 'DU4475') +
  geom_vline(xintercept = prod(t1Summary$M[c(2,3)])/100, color = 'red', lty = 2)


t1Data <- data.frame(A = c(41.28, 35.8, 23.06, NA),  B= c(36.67, 31.37, 36.76, 41.46), C = c(14.35, 20.59, 27.47, 15.58), D = c(100.1,100,99.9, NA))
t1Data <- reshape2::melt(t1Data)
t1Data <- na.omit(t1Data)
levels(t1Data$variable) <- c('QL-XII-47', 'GSK-690693', 'Combination', 'Control')
t1Data$variable <- factor(t1Data$variable, c('Combination', 'QL-XII-47', 'GSK-690693', 'Control') )
t1Data$cellLine <- 'BT20'
t1Summary <- t1Data %>% group_by(variable) %>% summarize(M = mean(value), L = M - sd(value), U = M + sd(value))

T2 <- ggplot(t1Data, aes(value, variable)) + 
  geom_errorbarh(data = t1Summary, aes(y = variable, x = M, xmin = L, xmax = U), height = 0.2) + 
  geom_bar(position = "dodge", stat = "summary",fun = "mean", fill = c('gray30','gray50','gray70','gray90')) +
  geom_jitter(height = 0.1, alpha = 0.5, pch = 16) + theme_light() +
  ylab('Treatment') +
  xlab('Viability') +
  stat_compare_means(aes(label = ..p.signif..), 
                     comparisons = list(c(1,2), c(1,3), c(1,4)), 
                     method = 't.test', 
                     method.args = list(alternative = 'less'), 
                     hide.ns = TRUE, tip.length = 0.01) +
  theme(plot.title = element_text(face = 2)) +
  labs(title = 'BT20') +
  geom_vline(xintercept = prod(t1Summary$M[c(2,3)])/100, color = 'red', lty = 2)


t1Data <- data.frame(A = c(72.65, 90.84, 78.16),  B= c(76.94, 97.23, 79.75), C = c(58.59, 71.96, 61.70), D = c(100.1,100,99.9))
t1Data <- reshape2::melt(t1Data)
t1Data <- na.omit(t1Data)
levels(t1Data$variable) <- c('QL-XII-47', 'GSK-690693', 'Combination', 'Control')
t1Data$variable <- factor(t1Data$variable, c('Combination', 'QL-XII-47', 'GSK-690693', 'Control') )
t1Data$cellLine <- 'CAL120'
t1Summary <- t1Data %>% group_by(variable) %>% summarize(M = mean(value), L = M - sd(value), U = M + sd(value))

T3 <- ggplot(t1Data, aes(value, variable)) + 
  geom_errorbarh(data = t1Summary, aes(y = variable, x = M, xmin = L, xmax = U), height = 0.2) + 
  geom_bar(position = "dodge", stat = "summary",fun = "mean", fill = c('gray30','gray50','gray70','gray90')) +
  geom_jitter(height = 0.1, alpha = 0.5, pch = 16) + theme_light() +
  ylab('Treatment') +
  xlab('Viability') +
  stat_compare_means(aes(label = ..p.signif..), 
                     comparisons = list(c(1,2), c(1,3), c(1,4)), 
                     method = 't.test', 
                     method.args = list(alternative = 'less'), 
                     hide.ns = TRUE, tip.length = 0.01) +
  labs(title = 'CAL120')  +
  theme(plot.title = element_text(face = 2)) +
  labs(tag = 'F') +
  geom_vline(xintercept = prod(t1Summary$M[c(2,3)])/100, color = 'red', lty = 2)


png('../Figures/F4A.png', width = 4800 * 0.85, height = 4800 * 0.12, res = 300)
T3 | T2 | T1
dev.off()

viabilityData <- bind_rows(T1$data, T2$data, T3$data)
colnames(viabilityData) <- c('Treatment', 'Viability', 'Cell Line')
readr::write_csv(viabilityData, '../Results/S6_Viability.csv')
