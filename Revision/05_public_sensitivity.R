# Point 2. Public sensitivity data for QL-XII-47, GSK-690693 and related compounds.
suppressMessages(library(readxl)); source('functions.R')
sink('results/05_public_sensitivity.txt', split = TRUE)
si <- read.csv('public/sample_info.csv')
grp <- function(sub) ifelse(sub == 'ERneg_HER2neg', 'TNBC (ER-/HER2-)',
                     ifelse(sub %in% c('ERpos_HER2neg', 'ERneg_HER2pos', 'ERpos_HER2pos'), 'non-TNBC breast', NA))
br <- si[si$lineage == 'breast' & si$Sanger_Model_ID != '', ]
g1 <- read_excel('public/GDSC1_fitted_dose_response_27Oct23.xlsx')
x <- g1[g1$DRUG_ID %in% c(235, 1155, 326), ]
x$IC50_uM <- exp(x$LN_IC50)
x$cell <- br$stripped_cell_line_name[match(x$SANGER_MODEL_ID, br$Sanger_Model_ID)]
x$group <- grp(br$lineage_sub_subtype[match(x$SANGER_MODEL_ID, br$Sanger_Model_ID)])
cat('== GDSC1 drug entries\n'); print(as.data.frame(unique(x[, c('DRUG_ID', 'DRUG_NAME', 'PUTATIVE_TARGET', 'MAX_CONC')])))
gsum <- list(); gcell <- list()
for (id in c(235, 326)) {
  y <- x[x$DRUG_ID == id & !is.na(x$group), ]
  cat('\n== GDSC1 DRUG_ID', id, unique(y$DRUG_NAME), '\n')
  print(do.call(rbind, lapply(split(y, y$group), function(z) data.frame(n = nrow(z), median_AUC = round(median(z$AUC), 3), median_IC50_uM = signif(median(z$IC50_uM), 3)))))
  wp <- wilcox.test(y$AUC[y$group == 'TNBC (ER-/HER2-)'], y$AUC[y$group == 'non-TNBC breast'])$p.value
  cat('Wilcoxon AUC, TNBC vs non-TNBC breast: P =', signif(wp, 3), '\n')
  thr <- if (id == 235) 0.6 else 0.8
  t <- y[y$group == 'TNBC (ER-/HER2-)', ]
  cat('TNBC lines with IC50 below', thr, 'uM:', sum(t$IC50_uM < thr), 'of', nrow(t), paste0('(', paste(t$cell[t$IC50_uM < thr], collapse = ', '), ')'), '\n')
  gsum[[length(gsum) + 1]] <- data.frame(drug_id = id, drug = unique(y$DRUG_NAME), group = c('TNBC', 'non-TNBC breast'),
    n = c(sum(y$group == 'TNBC (ER-/HER2-)'), sum(y$group == 'non-TNBC breast')),
    median_AUC = c(median(y$AUC[y$group == 'TNBC (ER-/HER2-)']), median(y$AUC[y$group == 'non-TNBC breast'])),
    wilcox_p = wp, test_conc_uM = thr, TNBC_lines_IC50_below_test_conc = sum(t$IC50_uM < thr))
  z <- y[y$cell %in% c('BT20', 'CAL120', 'DU4475', 'HCC70', 'HCC1806', 'MDAMB231', 'HS578T'), c('cell', 'AUC', 'IC50_uM', 'Z_SCORE')]
  print(as.data.frame(z[order(z$cell), ]), row.names = FALSE)
  gcell[[length(gcell) + 1]] <- data.frame(drug = unique(y$DRUG_NAME), as.data.frame(z))
}
write.csv(do.call(rbind, gsum), 'results/gdsc_summary.csv', row.names = FALSE)
write.csv(do.call(rbind, gcell), 'results/gdsc_cells.csv', row.names = FALSE)
cat('\n== HMS LINCS 20344 (Breast Cancer Profiling)\n')
h <- as.data.frame(read_excel('public/Screen20344_DrugSensitivity2.xlsx', sheet = 'Data'))
h <- h[, c('Cell Line', 'Drug Name', 'GRmax', 'GR_AOC')]; colnames(h) <- c('cell', 'drug', 'GRmax', 'GR_AOC')
h$group <- ifelse(h$cell %in% c('MCF10A', 'HME1'), 'non-tumorigenic', grp(si$lineage_sub_subtype[match(h$cell, si$stripped_cell_line_name)]))
cat('Lines excluded (no subtype annotation):', paste(sort(unique(h$cell[is.na(h$group)])), collapse = ', '), '\n')
hsum <- list()
for (dr in c('Ipatasertib/GDC0068', 'Torin2')) {
  y <- h[h$drug == dr & !is.na(h$group), ]
  cat('\n', dr, '\n')
  print(do.call(rbind, lapply(split(y, y$group), function(z) data.frame(n = nrow(z), median_GRmax = round(median(z$GRmax), 3), median_GR_AOC = round(median(z$GR_AOC), 3)))))
  print(y[y$group == 'non-tumorigenic', c('cell', 'GRmax', 'GR_AOC')], row.names = FALSE)
  hsum[[length(hsum) + 1]] <- rbind(aggregate(GR_AOC ~ group, y, median), data.frame(group = y$cell[y$group == 'non-tumorigenic'], GR_AOC = y$GR_AOC[y$group == 'non-tumorigenic']))
  hsum[[length(hsum)]]$n <- c(as.vector(table(y$group)[hsum[[length(hsum)]]$group[1:3]]), 1, 1)
  hsum[[length(hsum)]]$drug <- dr
}
write.csv(do.call(rbind, hsum), 'results/hms20344_summary.csv', row.names = FALSE)
cat('\n== HMS LINCS 20367 (QL-XII-47 and Torin2 in HCC70 and HCC1806, mean of replicates)\n')
q <- read.csv('../Data/dataset_20367_20210819011500.csv', check.names = FALSE)
q <- q[q[['Drug name Name']] %in% c('QL-XII-47', 'Torin2'), c('Cell line Name', 'Drug name Name', 'GRmax', 'GR50', 'GR_AOC')]
q <- aggregate(. ~ `Cell line Name` + `Drug name Name`, q, mean)
print(q)
write.csv(q, 'results/hms20367_summary.csv', row.names = FALSE)
sink()
