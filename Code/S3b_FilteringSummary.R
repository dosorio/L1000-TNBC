# Fraction of profiles removed by the correlation threshold at steps 1 and 2 (exact matching),
# and a check that Results/S1_Profiles.csv equals the output of the corrected retriever package.
suppressMessages({library(ccdata); library(preprocessCore)})
data('l1000_es')
n <- colnames(l1000_es); s <- strsplit(n, '_')
md <- data.frame(compound = sapply(s, function(X) paste0(X[1:(length(X)-3)], collapse = '_')),
                 cellLine = sapply(s, function(X) X[length(X)-2]), conc = sapply(s, function(X) X[length(X)-1]), name = n)
md <- md[md$cellLine %in% c('BT20', 'HS578T', 'MDAMB231'),]
pass <- function(X, thr = 0.6){ X <- as.matrix(X); if (ncol(X) > 1) X <- normalize.quantiles(X); cor(data.frame(rowMeans(X), X), method = 'sp')[,1][-1] > thr }
k1 <- paste(md$compound, md$cellLine, md$conc, sep = '_')
p1 <- unlist(lapply(split(md$name, k1), function(nm) pass(l1000_es[, nm, drop = FALSE])))
cat('Step 1: profiles', length(p1), '| removed', sum(!p1), sprintf('(%.2f%%)', 100 * mean(!p1)), '\n')
S1 <- read.csv('../Results/S1_CellLinesConcentrationProfiles.csv', row.names = 1, check.names = FALSE)
p2 <- unlist(lapply(split(colnames(S1), sub('_[^_]+$', '', colnames(S1))), function(nm) pass(S1[, nm, drop = FALSE])))
cat('Step 2: profiles', length(p2), '| removed', sum(!p2), sprintf('(%.2f%%)', 100 * mean(!p2)), '\n')
P <- as.matrix(read.csv('../Results/S1_Profiles.csv', row.names = 1, check.names = FALSE))
R <- suppressMessages(retriever::retriever(c('BT20', 'HS578T', 'MDAMB231')))
cat('Profiles:', ncol(P), '| identical to retriever package:', identical(sort(colnames(P)), sort(colnames(R))) && max(abs(P[rownames(R), colnames(R)] - R)) < 1e-12, '\n')
cat('Compounds retained:', ncol(P), 'of', length(unique(md$compound)), sprintf('(%.2f%%)', 100 * ncol(P) / length(unique(md$compound))), '| removed:', length(unique(md$compound)) - ncol(P), '\n')
