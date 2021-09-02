load('../Data/EpithelialCells.RData')
library(Nebulosa)

P <- plot_density(breastData[,breastData$diseaseStatus == 'Cancer'], c('BTK','BMX', 'TEC'), joint = TRUE, method = 'ks')
P[[1]] <- P[[1]] + theme_bw() + theme(legend.position = 'none', plot.title = element_text(face = 2))
P[[2]] <- P[[2]] + theme_bw() + theme(legend.position = 'none', plot.title = element_text(face = 2))
P[[3]] <- P[[3]] + theme_bw() + theme(legend.position = 'none', plot.title = element_text(face = 2))
P[[4]] <- P[[4]] + theme_bw() + theme(legend.position = 'none', plot.title = element_text(face = 2))

png('../Figures/S3.png', width = 1500, height = 1500, res = 300)
P
dev.off()
