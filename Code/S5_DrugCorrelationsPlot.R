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
# Label only the 10 strongest reversed genes on each side (|fold-change x effect size|) to avoid overlaps
labelScore <- abs(df$SC * df$DP)
for (cl in c('red', 'blue')) { i <- which(df$C == cl); df$G[i[rank(-labelScore[i]) > 10]] <- NA }

F3A <- ggplot(df, aes(SC, DP, label = G)) + 
  geom_abline(slope = -1, intercept = 0, lty = 2, color = 'red') +
  geom_point(pch = 16, alpha = df$alpha, color = df$C) + 
  geom_density_2d() + 
  geom_text_repel(min.segment.length = 0, fontface = 3, size = 3.5, bg.color = 'white', box.padding = 0.8, point.padding = 0.3, force = 6, max.overlaps = Inf, seed = 42) +
  theme_bw() +
  scale_x_continuous(expand = expansion(mult = 0.08)) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  labs(tag = 'A', title = 'QL-XII-47', subtitle = eLabel) +
  theme(plot.title = element_text(face = 2)) +
  xlab(expression(log[2]~(Fold-Change~Single-Cell~RNA-seq))) +
  ylab(expression(Effect~Sizes~LINCS~L1000~Project)) +
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
# Label only the 10 strongest reversed genes on each side (|fold-change x effect size|) to avoid overlaps
labelScore <- abs(df$SC * df$DP)
for (cl in c('red', 'blue')) { i <- which(df$C == cl); df$G[i[rank(-labelScore[i]) > 10]] <- NA }

F3C <- ggplot(df, aes(SC, DP, label = G)) + 
  geom_abline(slope = -1, intercept = 0, lty = 2, color = 'red') +
  geom_point(pch = 16, alpha = df$alpha, color = df$C) + 
  geom_density2d() + 
  geom_text_repel(min.segment.length = 0, fontface = 3, size = 3.5, bg.color = 'white', box.padding = 0.8, point.padding = 0.3, force = 6, max.overlaps = Inf, seed = 42) +
  theme_bw() +
  scale_x_continuous(expand = expansion(mult = 0.08)) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  labs(tag = 'D', title = 'QL-XII-47 + GSK-690693', subtitle = eLabel) +
  theme(plot.title = element_text(face = 2)) +
  xlab(expression(log[2]~(Fold-Change~Single-Cell~RNA-seq))) +
  ylab(expression(Effect~Sizes~LINCS~L1000~Combinations)) +
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

# Cell viability (Figure 5F): Bliss independence expectation and
# two-sided Welch t-tests of the combination against each single agent, Holm-adjusted across cell lines.
library(dplyr)
library(ggpubr)
library(patchwork)
viability <- read.csv('../Results/S6_CellularViability.csv', check.names = FALSE)
colnames(viability) <- c('Replicate', 'Treatment', 'Viability', 'CellLine')
viability$Treatment <- factor(viability$Treatment, c('Combination', 'QL-XII-47', 'GSK-690693', 'Control'))
tests <- do.call(rbind, lapply(split(viability, viability$CellLine), function(d){
  cmb <- d$Viability[d$Treatment == 'Combination']
  do.call(rbind, lapply(c('QL-XII-47', 'GSK-690693'), function(s)
    data.frame(CellLine = d$CellLine[1], group1 = 'Combination', group2 = s, p = t.test(cmb, d$Viability[d$Treatment == s])$p.value)))
}))
tests <- tests %>% group_by(group2) %>% mutate(p.adj = p.adjust(p, 'holm')) %>% ungroup() %>%
  mutate(p.adj.signif = cut(p.adj, c(0, 1e-4, 1e-3, 1e-2, 0.05, 1), c('****', '***', '**', '*', 'ns')))
viabilityPlot <- function(cl, tag = NULL){
  d <- viability[viability$CellLine == cl, ]
  s <- d %>% group_by(Treatment) %>% summarize(M = mean(Viability), L = M - sd(Viability), U = M + sd(Viability))
  bliss <- prod(s$M[s$Treatment %in% c('QL-XII-47', 'GSK-690693')]) / 100
  tt <- as.data.frame(tests[tests$CellLine == cl, ]); tt$y <- c(106, 116)
  tt$x1 <- match(tt$group1, levels(d$Treatment)); tt$x2 <- match(tt$group2, levels(d$Treatment))
  ggplot(d, aes(Treatment, Viability)) +
    geom_col(data = s, aes(Treatment, M), fill = c('gray30', 'gray50', 'gray70', 'gray90'), inherit.aes = FALSE) +
    geom_errorbar(data = s, aes(x = Treatment, ymin = L, ymax = U), width = 0.2, inherit.aes = FALSE) +
    geom_jitter(width = 0.1, height = 0, alpha = 0.6, pch = 16) +
    geom_hline(yintercept = bliss, color = 'red', lty = 2) +
    # significance brackets drawn manually so that labels stay horizontal
    geom_segment(data = tt, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE) +
    geom_segment(data = tt, aes(x = x1, xend = x1, y = y, yend = y - 2), inherit.aes = FALSE) +
    geom_segment(data = tt, aes(x = x2, xend = x2, y = y, yend = y - 2), inherit.aes = FALSE) +
    geom_text(data = tt, aes(x = (x1 + x2) / 2, y = y + 1.5, label = p.adj.signif), hjust = 0, size = 3.2, inherit.aes = FALSE) +
    scale_y_continuous(limits = c(0, 128), breaks = seq(0, 100, 25), expand = expansion(mult = c(0, 0.02))) +
    coord_flip(clip = 'off') + theme_light() +
    labs(title = cl, tag = tag, x = 'Treatment', y = 'Viability (% of control)') +
    theme(plot.title = element_text(face = 2), plot.margin = margin(5.5, 20, 5.5, 5.5))
}
png('../Figures/F4A.png', width = 4800 * 0.85, height = 4800 * 0.17, res = 300)
viabilityPlot('CAL120', tag = 'F') | viabilityPlot('BT20') | viabilityPlot('DU4475')
dev.off()
write.csv(as.data.frame(tests), '../Results/S6_ViabilityTests.csv', row.names = FALSE)
