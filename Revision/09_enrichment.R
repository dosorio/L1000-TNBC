# Point 5. GO Biological Process and Hallmark GSEA for QL-XII-47, GSK-690693, their combination
# and the TNBC signature, all on the genes shared by the signature and LINCS-L1000.
suppressMessages(library(fgsea)); source('functions.R')
sink('results/09_enrichment.txt', split = TRUE)
P <- as.matrix(read.csv('../Results/S1_Profiles.csv', row.names = 1, check.names = FALSE))
fc <- loadDisease(rownames(P)); g <- names(fc); cat('Shared genes:', length(g), '\n')
combo <- (P[, 'QL-XII-47'] + P[, 'GSK-690693']) / 2
st <- list(`TNBC signature` = fc, `QL-XII-47` = P[g, 'QL-XII-47'], `GSK-690693` = P[g, 'GSK-690693'], `QL-XII-47 + GSK-690693` = combo[g])
run <- function(sets) do.call(rbind, lapply(names(st), function(n){
  set.seed(1); r <- suppressWarnings(fgseaMultilevel(sets, st[[n]], minSize = 10, maxSize = 500))
  data.frame(pathway = r$pathway, NES = r$NES, padj = r$padj, size = r$size, input = n,
             leadingEdge = sapply(r$leadingEdge, paste, collapse = ','))
}))
wide <- function(R){
  W <- reshape(R[, c('pathway', 'input', 'NES', 'padj')], idvar = 'pathway', timevar = 'input', direction = 'wide')
  colnames(W) <- sub('^NES\\.(.*)$', '\\1 NES', sub('^padj\\.(.*)$', '\\1 FDR', colnames(W))); W
}
BP <- suppressWarnings(gmtPathways('public/GO_BP_2023.gmt'))
R <- run(BP); write.csv(R, 'results/GO_BP_GSEA.csv', row.names = FALSE)
drugs <- names(st)[-1]
cat('GO BP terms with FDR < 0.05 per input:\n'); print(table(factor(R$input[R$padj < 0.05], levels = names(st))))
sel <- unique(R$pathway[R$input %in% drugs & R$padj < 0.05])
W <- wide(R[R$pathway %in% sel, ])
cs <- W[W$`QL-XII-47 + GSK-690693 FDR` < 0.05, ]
cat('Combination terms:', nrow(cs), '| opposite sign to TNBC signature:', sum(sign(cs$`TNBC signature NES`) != sign(cs$`QL-XII-47 + GSK-690693 NES`)),
    '| minimum TNBC-signature FDR among them:', signif(min(cs$`TNBC signature FDR`), 3), '\n')
write.csv(W, 'results/GO_BP_heatmap_table.csv', row.names = FALSE)
cat('\nLeading edge of combination GO BP terms:\n')
lc <- R[R$input == 'QL-XII-47 + GSK-690693' & R$padj < 0.05, ]; lc <- lc[order(-lc$NES), ]
for (i in seq_len(nrow(lc))) cat(sprintf('%s | NES %.2f | FDR %.3g | %s\n', lc$pathway[i], lc$NES[i], lc$padj[i], paste(head(strsplit(lc$leadingEdge[i], ',')[[1]], 15), collapse = ',')))

H <- suppressWarnings(gmtPathways('public/Hallmark_2020.gmt'))
RH <- run(H); write.csv(RH, 'results/Hallmark_GSEA.csv', row.names = FALSE)
WH <- wide(RH)
cat('\nHallmarks, TNBC signature, ordered by FDR:\n')
print(format(WH[order(WH$`TNBC signature FDR`), c('pathway', 'TNBC signature NES', 'TNBC signature FDR', 'QL-XII-47 + GSK-690693 NES', 'QL-XII-47 + GSK-690693 FDR')][1:8, ], digits = 2), row.names = FALSE)
cat('\nHallmarks significant for the combination:\n')
print(format(WH[WH$`QL-XII-47 + GSK-690693 FDR` < 0.05, c('pathway', 'TNBC signature NES', 'TNBC signature FDR', 'QL-XII-47 + GSK-690693 NES', 'QL-XII-47 + GSK-690693 FDR')], digits = 2), row.names = FALSE)
set.seed(1); r1001 <- suppressWarnings(fgseaMultilevel(H, combo))    # default sizes, as in the published Fig. 5E code (S9_PRLComparison.R)
cat('\nCombination Hallmarks on all', length(combo), 'LINCS genes (FDR < 0.05, as in Fig. 5E):\n')
print(as.data.frame(r1001[r1001$padj < 0.05, c('pathway', 'NES', 'padj')])[order(r1001$NES[r1001$padj < 0.05]), ], row.names = FALSE)

genes <- c('GADD45B', 'GADD45A', 'HMOX1', 'CEBPD', 'CDKN1A', 'CDKN1B', 'FOXO3', 'ERBB3', 'IGF1R', 'EGFR', 'MYC', 'JUN', 'EGR1', 'FOSL1', 'NR3C1', 'HIF1A', 'SMAD3')
cat('\nGenes discussed in the response:\n')
dd <- read.csv('../Data/de_EC_TNBC-H.csv', row.names = 1)
print(data.frame(TNBC_log2FC = round(dd[genes, 'avg_log2FC'], 2), TNBC_FDR = signif(dd[genes, 'p_val_adj'], 2),
                 QL = round(P[genes, 'QL-XII-47'], 2), GSK = round(P[genes, 'GSK-690693'], 2), combination = round(combo[genes], 2), row.names = genes))
sink()
