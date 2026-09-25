# Shared helpers for the revision-3 analyses.
suppressMessages(library(preprocessCore))

TN <- c('BT20', 'HS578T', 'MDAMB231')    # LINCS-L1000 TNBC cell lines
normName <- function(s) toupper(gsub('[^A-Za-z0-9]', '', s))    # drug-name matching across databases

# Loads the LINCS-L1000 TNBC signatures (ccdata) and their metadata (cached).
loadTNBC <- function(cache = 'results/tnbc_l1000.rds'){
  if (file.exists(cache)) return(readRDS(cache))
  suppressMessages(library(ccdata))
  data('l1000_es', package = 'ccdata', envir = environment())
  colnames(l1000_es) <- gsub('\\[|\\]', '', colnames(l1000_es))
  s <- strsplit(colnames(l1000_es), '_')
  md <- data.frame(compound = sapply(s, function(X) paste0(X[1:(length(X) - 3)], collapse = '_')),
                   cellLine = sapply(s, function(X) X[length(X) - 2]),
                   conc = sapply(s, function(X) X[length(X) - 1]),
                   time = sapply(s, function(X) X[length(X)]),
                   name = colnames(l1000_es), stringsAsFactors = FALSE)
  md <- md[md$cellLine %in% TN, ]
  out <- list(E = l1000_es[, md$name], md = md)
  saveRDS(out, cache)
  out
}

# Single-cell TNBC signature restricted to the genes measured in LINCS-L1000.
loadDisease <- function(genes){
  dis <- read.csv('../Data/de_EC_TNBC-H.csv', row.names = 1)
  g <- intersect(rownames(dis), genes)
  setNames(dis[g, 'avg_log2FC'], g)
}

# retriever: one collapsing step (quantile normalization, keep profiles with rho > thr with the mean).
collapse <- function(X, thr){
  X <- as.matrix(X)
  if (ncol(X) > 1) X <- normalize.quantiles(X)
  keep <- cor(data.frame(rowMeans(X), X), method = 'spearman')[, 1][-1] > thr
  if (sum(keep) > 1) rowMeans(X[, keep, drop = FALSE]) else rep(NA, nrow(X))
}

# retriever (time, then concentration, then cell line) with exact metadata matching.
retrieverExact <- function(E, md, thr = 0.6, cells = TN){
  md <- md[md$cellLine %in% cells, ]
  k1 <- paste(md$compound, md$cellLine, md$conc, sep = '\r')
  S1 <- sapply(split(md$name, k1), function(nm) collapse(E[, nm, drop = FALSE], thr))
  S1 <- S1[, colSums(is.na(S1)) == 0, drop = FALSE]
  p1 <- do.call(rbind, strsplit(colnames(S1), '\r'))
  S2 <- sapply(split(seq_len(ncol(S1)), paste(p1[, 1], p1[, 2], sep = '\r')), function(i) collapse(S1[, i, drop = FALSE], thr))
  S2 <- S2[, colSums(is.na(S2)) == 0, drop = FALSE]
  p2 <- do.call(rbind, strsplit(colnames(S2), '\r'))
  S3 <- sapply(split(seq_len(ncol(S2)), p2[, 1]), function(i) collapse(S2[, i, drop = FALSE], thr))
  S3 <- S3[, colSums(is.na(S3)) == 0, drop = FALSE]
  rownames(S3) <- rownames(E)
  S3
}

# Naive average of all experiments per compound.
naiveMean <- function(E, md){
  M <- sapply(split(md$name, md$compound), function(nm) rowMeans(E[, nm, drop = FALSE]))
  rownames(M) <- rownames(E)
  M
}

# Prototype Ranked List (Iorio et al., PNAS 2010): iterative Borda merging of the closest pair of
# ranked lists (Spearman footrule), first within each cell line, then across cell lines.
bordaMerge <- function(R){
  R <- as.matrix(R)
  while (ncol(R) > 1){
    D <- as.matrix(dist(t(R), method = 'manhattan')); diag(D) <- Inf
    ij <- which(D == min(D), arr.ind = TRUE)[1, ]
    R <- cbind(R[, -ij, drop = FALSE], rank(rowMeans(R[, ij]), ties.method = 'first'))
  }
  R[, 1]
}
prl <- function(E, md){
  P <- sapply(split(md, md$compound), function(d){
    perCell <- sapply(split(d$name, d$cellLine), function(nm) bordaMerge(apply(-E[, nm, drop = FALSE], 2, rank, ties.method = 'first')))
    -bordaMerge(perCell)    # rank 1 = most upregulated; sign flipped so larger = more upregulated
  })
  rownames(P) <- rownames(E)
  P
}

# Spearman correlation of each column with the disease signature (negative = predicted reversal).
connectivity <- function(M, fc) apply(M[names(fc), , drop = FALSE], 2, function(x) cor(x, fc, method = 'spearman'))

# metaLINCS drug-level NES (negative = predicted reversal).
metaLincsScore <- function(E, md, fc, nmin = 10){
  suppressMessages(library(metaLINCS))
  # metaLINCS takes the drug name as the text before the first '@' or '_', so underscores in
  # compound names (e.g. ligands such as 'EGF_lig') are protected and restored afterwards
  safe <- gsub('_', '.', md$compound)
  mD <- E; colnames(mD) <- paste0(safe, '@', md$cellLine, '_', md$conc, '_', md$time)
  # names = is passed as well: metaLINCS 0.9.0 stops on a named vector without it
  ml <- suppressMessages(computeConnectivityEnrichment(fc, names = names(fc), mDrugEnrich = mD, nmin = nmin, nprune = 0))
  setNames(ml$X[, 1], md$compound[match(rownames(ml$X), safe)])
}

# All benchmark scores for a set of experiments.
allScores <- function(E, md, fc, cells = TN){
  m <- md[md$cellLine %in% cells, ]; X <- E[, m$name]
  list(retriever = connectivity(retrieverExact(X, m, cells = cells), fc),
       naive_average = connectivity(naiveMean(X, m), fc),
       PRL = connectivity(prl(X, m), fc),
       metaLINCS = metaLincsScore(X, m, fc))
}
