# Figures: benchmark and GO BP heatmap.
suppressMessages({library(ggplot2); library(patchwork)})
dir.create('figures', showWarnings = FALSE)
lab <- c(retriever = 'retriever', naive_average = 'Naive average', PRL = 'PRL (Iorio 2010)', metaLINCS = 'metaLINCS')
b <- read.csv('results/benchmark_results.csv', row.names = 1)[names(lab), ]
b$method <- factor(lab, levels = rev(lab))
s <- read.csv('results/stability_results.csv', row.names = 1)[names(lab), ]
s$method <- factor(lab, levels = rev(lab))
n <- readRDS('results/benchmark_scores.rds')$common
nc <- as.integer(sub('.*: ', '', grep('^Compounds scored by every method', readLines('results/04_stability.txt'), value = TRUE)))
A <- ggplot(b, aes(AUC_rho, method)) + geom_vline(xintercept = 0, lty = 2) +
  geom_errorbarh(aes(xmin = AUC_lo, xmax = AUC_hi), height = 0.2) + geom_point(size = 2.5) + theme_bw() +
  labs(tag = 'A', x = expression(Spearman~rho~'(predicted reversal vs GDSC1 sensitivity)'), y = NULL, title = 'Accuracy',
       subtitle = paste(length(n), 'compounds, 21 held-out ER-/HER2- breast cancer lines'))
B <- ggplot(s, aes(mean_rank_correlation, method)) + geom_col(fill = 'gray60') + coord_cartesian(xlim = c(0.5, 1)) + theme_bw() +
  labs(tag = 'B', x = expression('Mean'~rho~'between leave-one-cell-line-out rankings'), y = NULL, title = 'Robustness',
       subtitle = paste(nc, 'compounds scored by all methods'))
ggsave('figures/Fig_R3_benchmark.png', A + B, width = 12, height = 3.8, dpi = 300, bg = 'white')

R <- read.csv('results/GO_BP_GSEA.csv', check.names = FALSE)
inputs <- c('TNBC signature', 'QL-XII-47', 'GSK-690693', 'QL-XII-47 + GSK-690693')
sel <- unique(R$pathway[R$input != 'TNBC signature' & R$padj < 0.05])
L <- R[R$pathway %in% sel, ]
ord <- L[L$input == 'QL-XII-47 + GSK-690693', ]; L$pathway <- factor(L$pathway, levels = ord$pathway[order(ord$NES)])
L$input <- factor(L$input, levels = inputs); L$star <- ifelse(L$padj < 0.05, '*', '')
D <- ggplot(L, aes(input, pathway, fill = NES)) + geom_tile(color = 'white') + geom_text(aes(label = star), vjust = 0.75) +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(), plot.title.position = 'plot') +
  labs(x = NULL, y = NULL, fill = 'NES', title = 'GO Biological Process (GSEA)', subtitle = '610 genes shared by the TNBC signature and LINCS-L1000; * FDR < 0.05')
ggsave('figures/Fig_R3_GO_BP_heatmap.png', D, width = 9, height = max(6, 0.22 * length(sel) + 2), dpi = 300, bg = 'white')
