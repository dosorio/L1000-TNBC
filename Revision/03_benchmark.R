# Point 1. Accuracy against GDSC1 sensitivity in TNBC lines not profiled by LINCS-L1000.
suppressMessages(library(readxl)); source('functions.R')
set.seed(1)
sink('results/03_benchmark.txt', split = TRUE)
d <- loadTNBC(); E <- d$E; md <- d$md
sig <- readRDS('results/signatures.rds')
fc <- loadDisease(rownames(E)); cat('Genes shared by the disease signature and LINCS:', length(fc), '\n')
cat('Minimum experiments per compound:', min(table(md$compound)), '\n')
scores <- list(retriever = connectivity(sig$Rfix, fc),
               retriever_published = connectivity(sig$Rpub, fc),
               naive_average = connectivity(sig$Naive, fc),
               PRL = connectivity(sig$PRL, fc),
               metaLINCS = metaLincsScore(E, md, fc))
cat('Compounds scored:\n'); print(sapply(scores, length))
cat('\nRank (1 = strongest predicted reversal):\n')
print(sapply(scores, function(s) { r <- rank(s); c(QL_XII_47 = unname(r['QL-XII-47']), GSK_690693 = unname(r['GSK-690693'])) }))
cat('\nTop 5 per method:\n'); print(sapply(scores, function(s) names(sort(s))[1:5]))

g1 <- read_excel('public/GDSC1_fitted_dose_response_27Oct23.xlsx')
si <- read.csv('public/sample_info.csv'); si <- si[si$lineage == 'breast' & si$Sanger_Model_ID != '', ]
held <- si$Sanger_Model_ID[si$lineage_sub_subtype == 'ERneg_HER2neg' & !si$stripped_cell_line_name %in% TN]
tn <- g1[g1$SANGER_MODEL_ID %in% held, ]; tn$drug <- normName(tn$DRUG_NAME)
cat('\nHeld-out ER-/HER2- breast lines in GDSC1:', length(unique(tn$SANGER_MODEL_ID)), '\n')
a <- aggregate(cbind(AUC, Z_SCORE) ~ drug + SANGER_MODEL_ID, tn, mean)          # duplicate GDSC drug IDs averaged
gt <- aggregate(cbind(AUC, Z_SCORE) ~ drug, a, median)                           # median across held-out lines
nl <- table(a$drug); gt <- gt[gt$drug %in% names(nl)[nl >= 5], ]
common <- intersect(Reduce(intersect, lapply(scores, function(s) normName(names(s)))), gt$drug)
cat('Compounds scored by all methods and in GDSC1:', length(common), '\n')
res <- t(sapply(names(scores), function(m){
  s <- scores[[m]]; set.seed(1)    # same bootstrap draws for every method, independent of which methods are included
  s <- s[match(common, normName(names(s)))]; y <- gt[match(common, gt$drug), ]
  out <- c()
  for (v in c('AUC', 'Z_SCORE')) {
    ct <- cor.test(s, y[[v]], method = 'spearman', exact = FALSE)
    b <- replicate(2000, { i <- sample(length(s), replace = TRUE); cor(s[i], y[[v]][i], method = 'spearman') })
    out <- c(out, setNames(c(ct$estimate, quantile(b, c(0.025, 0.975)), ct$p.value), paste0(v, c('_rho', '_lo', '_hi', '_p'))))
  }
  q <- y$AUC <= quantile(y$AUC, 0.25)
  c(out, precision_at_10 = mean(q[order(s)[1:10]]), random_expectation = mean(q))
}))
cat('\nSpearman correlation of predicted reversal with GDSC1 sensitivity (positive = more reversal, more sensitive):\n')
print(signif(res, 3))
write.csv(res, 'results/benchmark_results.csv')
cm <- sapply(scores, function(s) s[match(common, normName(names(s)))])
cat('\nAgreement between methods (Spearman):\n'); print(round(cor(cm, method = 'spearman'), 2))
gt$potency_rank <- rank(gt$AUC); g <- gt[gt$drug %in% common, ]; g$potency_rank <- rank(g$AUC)
cat('\nPotency rank among the', length(common), 'compounds (1 = lowest median AUC):\n'); print(g[g$drug %in% c('QLXII47', 'GSK690693'), ])
# Compounds removed by the retriever consistency filter
kept <- normName(names(scores$retriever)); allc <- normName(names(scores$naive_average))
cat('\nCompounds removed by retriever:', length(allc) - length(kept), '\n')
for (grp in c('kept', 'removed')) {
  ids <- if (grp == 'kept') intersect(kept, gt$drug) else intersect(setdiff(allc, kept), gt$drug)
  for (m in c('naive_average', 'PRL', 'metaLINCS')) {
    ct <- cor.test(scores[[m]][match(ids, normName(names(scores[[m]])))], gt$AUC[match(ids, gt$drug)], method = 'spearman', exact = FALSE)
    cat(sprintf('%-8s n = %3d  %-24s rho = %6.3f  P = %.3g\n', grp, length(ids), m, ct$estimate, ct$p.value))
  }
}
saveRDS(list(scores = scores, gt = gt, common = common), 'results/benchmark_scores.rds')
sink()
