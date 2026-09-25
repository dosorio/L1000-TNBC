# Point 5. Direction of the 205 differentially expressed genes (Fig. 3).
sink('results/08_deg_direction.txt', split = TRUE)
d <- read.csv('../Data/de_EC_TNBC-H.csv', row.names = 1)
s <- d[abs(d$avg_log2FC) > 1 & d$p_val_adj < 0.05, ]
cat('DEGs:', nrow(s), '| up in TNBC:', sum(s$avg_log2FC > 0), '| down in TNBC:', sum(s$avg_log2FC < 0), '\n')
cat('Sign check (positive = higher in TNBC):\n'); print(round(d[c('TFF1', 'AGR2', 'ESR1', 'KRT14', 'S100A8', 'VIM'), 'avg_log2FC'], 2))
up <- s[order(-s$avg_log2FC), ]; up <- up[up$avg_log2FC > 0, ]
dn <- s[order(s$avg_log2FC), ]; dn <- dn[dn$avg_log2FC < 0, ]
cat('Top 25 up:', paste(rownames(up)[1:25], collapse = ', '), '\n')
cat('Top 25 down:', paste(rownames(dn)[1:25], collapse = ', '), '\n')
out <- data.frame(gene = rownames(s), s, direction = ifelse(s$avg_log2FC > 0, 'Up in TNBC', 'Down in TNBC'))
write.csv(out[order(-out$avg_log2FC), ], 'results/DEG_205_direction.csv', row.names = FALSE)
sink()
