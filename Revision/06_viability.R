# Point 2. Reanalysis of Supplementary Table 6 (Figure 5F) with Bliss independence.
set.seed(1)
sink('results/06_viability.txt', split = TRUE)
v <- read.csv('../Results/S6_CellularViability.csv', check.names = FALSE)
colnames(v) <- c('rep', 'trt', 'via', 'cell')
res <- do.call(rbind, lapply(split(v, v$cell), function(d){
  q <- d$via[d$trt == 'QL-XII-47'] / 100; k <- d$via[d$trt == 'GSK-690693'] / 100; cb <- d$via[d$trt == 'Combination'] / 100
  boot <- replicate(10000, mean(sample(cb, replace = TRUE)) - mean(sample(q, replace = TRUE)) * mean(sample(k, replace = TRUE)))
  data.frame(cell = d$cell[1], n_QL = length(q), n_GSK = length(k), n_combo = length(cb),
             QL_mean = 100 * mean(q), QL_sd = 100 * sd(q), GSK_mean = 100 * mean(k), GSK_sd = 100 * sd(k),
             combo_mean = 100 * mean(cb), combo_sd = 100 * sd(cb), bliss = 100 * mean(q) * mean(k),
             combo_minus_bliss = 100 * (mean(cb) - mean(q) * mean(k)),
             ci_lo = 100 * quantile(boot, 0.025), ci_hi = 100 * quantile(boot, 0.975),
             p_vs_QL = t.test(cb, q)$p.value, p_vs_GSK = t.test(cb, k)$p.value,
             CV_QL = 100 * sd(q) / mean(q), CV_combo = 100 * sd(cb) / mean(cb))
}))
res$p_vs_QL_holm <- p.adjust(res$p_vs_QL, 'holm'); res$p_vs_GSK_holm <- p.adjust(res$p_vs_GSK, 'holm')
print(format(res, digits = 3), row.names = FALSE)
write.csv(res, 'results/viability_reanalysis.csv', row.names = FALSE)
sink()
