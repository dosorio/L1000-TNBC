# Point 3. Sensitivity of retriever outputs to its correlation threshold.
source('functions.R')
sink('results/07_threshold.txt', split = TRUE)
d <- loadTNBC(); E <- d$E; md <- d$md; fc <- loadDisease(rownames(E)); rfc <- rank(fc)
out <- do.call(rbind, lapply(c(0.4, 0.5, 0.6, 0.7, 0.8), function(thr){
  R <- retrieverExact(E, md, thr = thr)[names(fc), , drop = FALSE]
  s <- connectivity(R, fc)
  cb <- combn(colnames(R), 2)
  cs <- apply(cb, 2, function(p) cor(rank((R[, p[1]] + R[, p[2]]) / 2), rfc))
  i <- which(apply(cb, 2, function(p) setequal(p, c('QL-XII-47', 'GSK-690693'))))
  data.frame(threshold = thr, compounds = ncol(R), combinations = ncol(cb), top_compound = names(which.min(s)),
             QL_XII_47_rank = if ('QL-XII-47' %in% names(s)) unname(rank(s)['QL-XII-47']) else NA,
             GSK_690693_present = 'GSK-690693' %in% names(s),
             top_combination = paste(cb[, which.min(cs)], collapse = ' + '),
             QL_GSK_rank = if (length(i)) rank(cs)[i] else NA)
}))
print(out, row.names = FALSE)
write.csv(out, 'results/threshold_sensitivity.csv', row.names = FALSE)
sink()
