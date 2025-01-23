library(Seurat)
library(fgsea)
library(GSVA)
library(Matrix)

H <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020')

load('../Data/EpithelialCells.RData')
breastData <- breastData[,breastData$diseaseStatus == 'Cancer']

drugProfiles <- read.csv('../Results/S1_Profiles.csv', row.names = 1)
dE <-  gsva(as.matrix(drugProfiles), H, method = 'ssgsea')


id <- unique(breastData$orig.ident)[1]



X <- breastData@assays$RNA@data[,breastData$orig.ident %in% id]
E <- gsva(X, H, method = 'ssgsea')
E <- as.matrix(E)
C <- cor(E, dE)

O <- cbind(sort(colMeans(C)))
sort(apply(C,2,mean)/apply(C,2,sd))
