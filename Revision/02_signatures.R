# Drug-level signatures for the benchmark, and a check of the published retriever profiles.
source('functions.R')
sink('results/02_signatures.txt', split = TRUE)
# Published (revision 2) results are read from the repository history (commit ad5332e), because
# Results/ now holds the profiles regenerated with exact matching.
published <- function(file, ...) read.csv(pipe(sprintf('git -C .. show ad5332e:Results/%s', file)), ...)
d <- loadTNBC(); E <- d$E; md <- d$md
Rpub <- as.matrix(published('S1_Profiles.csv', row.names = 1, check.names = FALSE))[rownames(E), ]
Rfix <- retrieverExact(E, md)
cat('Published retriever profiles:', ncol(Rpub), '| exact-matching profiles:', ncol(Rfix), '\n')
cat('Only in published:', setdiff(colnames(Rpub), colnames(Rfix)), '| only in exact:', setdiff(colnames(Rfix), colnames(Rpub)), '\n')
cc <- intersect(colnames(Rpub), colnames(Rfix))
cat('Shared profiles identical (max abs difference < 1e-12):', sum(apply(abs(Rpub[, cc] - Rfix[, cc]), 2, max) < 1e-12), 'of', length(cc), '\n')
cat('Published TNF equals exact TNF_lig:', max(abs(Rpub[, 'TNF'] - Rfix[, 'TNF_lig'])) < 1e-12, '\n')
# Source of the published EGF profile: substring matching also captured HBEGF_lig signatures
S1pub <- colnames(published('S1_CellLinesConcentrationProfiles.csv', row.names = 1, check.names = FALSE, nrows = 1))
cat('Step-1 profiles matched by grepl("EGF_lig_BT20"):', grep('EGF_lig_BT20', S1pub, fixed = TRUE, value = TRUE), '\n')
cat('Step-1 profiles matched by grepl("EGF_lig_HS578T"):', grep('EGF_lig_HS578T', S1pub, fixed = TRUE, value = TRUE), '\n')
dp <- published('S3_drugPotential.csv', row.names = 1)
cp <- published('S3_combinationPotential.csv', row.names = 1)
cat('Published rank of EGF (single agent):', which(rownames(dp) == 'EGF'), 'of', nrow(dp), '\n')
cat('Published top combinations:\n'); print(head(cp[, 1, drop = FALSE], 3))
saveRDS(list(Rpub = Rpub, Rfix = Rfix, Naive = naiveMean(E, md), PRL = prl(E, md)), 'results/signatures.rds')
sink()
