# Supplementary Figure: public drug-sensitivity data for the prioritized compounds (Point 2).
suppressMessages({library(readxl); library(ggplot2); library(ggrepel); library(patchwork)})
si <- read.csv('public/sample_info.csv')
grp <- function(sub) ifelse(sub == 'ERneg_HER2neg', 'TNBC', ifelse(sub %in% c('ERpos_HER2neg', 'ERneg_HER2pos', 'ERpos_HER2pos'), 'Non-TNBC breast', NA))
br <- si[si$lineage == 'breast' & si$Sanger_Model_ID != '', ]
g1 <- as.data.frame(read_excel('public/GDSC1_fitted_dose_response_27Oct23.xlsx'))
x <- g1[g1$DRUG_ID %in% c(235, 326), c('DRUG_ID', 'SANGER_MODEL_ID', 'AUC')]
x$drug <- ifelse(x$DRUG_ID == 235, 'QL-XII-47', 'GSK-690693')
x$cell <- br$stripped_cell_line_name[match(x$SANGER_MODEL_ID, br$Sanger_Model_ID)]
x$group <- grp(br$lineage_sub_subtype[match(x$SANGER_MODEL_ID, br$Sanger_Model_ID)])
x$group[is.na(x$cell)] <- 'Other cancer lines'
x <- x[!is.na(x$group), ]
x$group <- factor(x$group, c('TNBC', 'Non-TNBC breast', 'Other cancer lines'))
x$drug <- factor(x$drug, c('QL-XII-47', 'GSK-690693'))
x$label <- ifelse(x$cell %in% c('BT20', 'CAL120', 'DU4475'), x$cell, NA)
A <- ggplot(x, aes(group, AUC)) + geom_boxplot(outlier.shape = NA, fill = 'gray90', width = 0.45) +
  geom_point(data = x[x$group != 'Other cancer lines' & is.na(x$label), ], position = position_jitter(width = 0.1, height = 0, seed = 1), alpha = 0.5, pch = 16) +
  geom_point(data = x[!is.na(x$label), ], color = 'red', size = 2) +
  geom_text_repel(data = x[!is.na(x$label), ], aes(label = label), color = 'red', size = 3, nudge_x = 0.3, direction = 'y', hjust = 0,
                  min.segment.length = 0, segment.color = 'red', bg.color = 'white', bg.r = 0.15, seed = 1) +
  facet_wrap(~drug) + theme_bw() + theme(plot.title = element_text(face = 2)) +
  labs(tag = 'A', title = 'GDSC1 drug sensitivity', subtitle = 'Lower AUC = more sensitive', x = NULL, y = 'AUC')
h <- as.data.frame(read_excel('public/Screen20344_DrugSensitivity2.xlsx', sheet = 'Data'))[, c('Cell Line', 'Drug Name', 'GR_AOC')]
colnames(h) <- c('cell', 'drug', 'GR_AOC')
h <- h[h$drug %in% c('Ipatasertib/GDC0068', 'Torin2'), ]
h$drug <- ifelse(h$drug == 'Torin2', 'Torin2 (mTOR)', 'Ipatasertib (pan-AKT)')
h$group <- ifelse(h$cell %in% c('MCF10A', 'HME1'), 'Non-tumorigenic', grp(si$lineage_sub_subtype[match(h$cell, si$stripped_cell_line_name)]))
h <- h[!is.na(h$group), ]
h$group <- factor(h$group, c('TNBC', 'Non-TNBC breast', 'Non-tumorigenic'))
h$label <- ifelse(h$group == 'Non-tumorigenic', h$cell, NA)
B <- ggplot(h, aes(group, GR_AOC)) + geom_boxplot(outlier.shape = NA, fill = 'gray90', width = 0.45) +
  geom_point(position = position_jitter(width = 0.1, height = 0, seed = 1), alpha = 0.6, pch = 16) +
  geom_text_repel(aes(label = label), size = 3, nudge_x = 0.3, direction = 'y', hjust = 0, min.segment.length = 0,
                  bg.color = 'white', bg.r = 0.15, seed = 1, na.rm = TRUE) +
  facet_wrap(~drug) + theme_bw() + theme(plot.title = element_text(face = 2)) +
  labs(tag = 'B', title = 'HMS LINCS Breast Cancer Profiling (dataset 20344)', subtitle = 'Higher GR_AOC = stronger growth inhibition', x = NULL, y = expression(GR[AOC]))
dir.create('figures', showWarnings = FALSE)
ggsave('figures/Fig_R3_public_sensitivity.png', A / B, width = 9, height = 8, dpi = 300, bg = 'white')
