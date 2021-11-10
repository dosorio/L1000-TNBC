library(enrichR)
library(fgsea)
HM <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020')

X <- read.csv('Data/de_EC_TNBC-H.csv')
table(FC = (X$avg_log2FC) < -1, FDR = X$p_val_adj < 0.05)

upR <- X$X[X$avg_log2FC > 1 & X$p_val_adj < 0.05]
dnR <- X$X[X$avg_log2FC < -1 & X$p_val_adj < 0.05]

E1 <- do.call(rbind.data.frame, enrichr(upR, databases = 'MSigDB_Hallmark_2020'))
E1 <- E1[E1$Adjusted.P.value < 0.05,]

E2 <- do.call(rbind.data.frame, enrichr(dnR, databases = 'MSigDB_Hallmark_2020'))
E2 <- E2[E2$Adjusted.P.value < 0.05,]

FC <- X$avg_log2FC
names(FC) <- X$X

GSEA <- fgseaMultilevel(HM,FC, eps = 0)
GSEA <- GSEA[GSEA$padj < 0.05,]

E1 <- E1[E1$Term %in%GSEA$pathway[GSEA$NES > 0],]
E2 <- E2[E2$Term %in%GSEA$pathway[GSEA$NES < 0],]

writeLines(paste0(E1$Term, ' (FDR = ',gsub('e-0',' times 10^{-',formatC(E1$Adjusted.P.value, format = 'e', digits = 2)),'}; ', gsub(';',', ',E1$Genes), ')'), sep = ', ')
writeLines(paste0(E2$Term, ' (FDR = ',gsub('e-0',' times 10^{-',formatC(E2$Adjusted.P.value, format = 'e', digits = 2)),'}; ', gsub(';',', ',E2$Genes), ')'), sep = ', ')
