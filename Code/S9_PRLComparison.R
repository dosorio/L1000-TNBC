library(ggplot2)
library(ggrepel)
library(fgsea)
library(dplyr)
library(patchwork)

# Loading Hallmarks
H <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020')

# Loading drug profiles
DrugProfiles <- read.csv('Results/S1_Profiles.csv', row.names = 1, check.names = FALSE)

# Extracting QL-XII-47 and GSK-690693 profiles
P1 <- DrugProfiles[,'QL-XII-47']
P2 <- DrugProfiles[,'GSK-690693']

# Computing ranks for QL-XII-47 and GSK-690693 profiles
R1 <- rank(P1)
R2 <- rank(P2)

# Computing averages of profiles and ranks
P <- (P1 + P2)/2
R <- (R1 + R2)/2

# Assigning names
names(P) <- rownames(DrugProfiles)
names(R) <- rownames(DrugProfiles)

# Creating a data.frame
df <- data.frame(Ranks = R, ES = P, Label = rownames(DrugProfiles))
df <- df[order(abs(df$Ranks - rank(df$ES)), decreasing = TRUE),]
df$Label[50:nrow(df)] <- NA

# Plotting
A <- ggplot(df, aes(ES, Ranks, label = Label)) +
  geom_point() +
  theme_minimal() +
  ylab('Average Ranks') +
  xlab('Average Effect Sizes') +
  labs(title = 'QL-XII-47 + GSK-690693', subtitle = parse(text = paste0('rho == ', round(cor(R, P, method = 'sp'), digits = 2)))) +
  geom_text_repel(segment.size = 0.1, min.segment.length = 0, fontface = 3, size = 2, bg.color = 'white')

# Computing GSEA scores
set.seed(1)
EP <- fgseaMultilevel(pathways = H, stats = P)
EP$input <- 'Average Effect Sizes'
set.seed(1)
ER <- fgseaMultilevel(pathways = H, stats = R)
ER$input <- 'Average Ranks'

# Binding results
E <- bind_rows(EP, ER)
# Ordering by effect and significance
E <- E[order(-log10(E$padj)*E$NES, decreasing = FALSE),]
# Filtering non-significant
E <- E[E$padj < 0.05,]
# Creating an order for display purposes
E$pathway <- factor(E$pathway, levels = unique(E$pathway))

# Plotting
B <- ggplot(E, aes(x = NES*-log10(padj), y = pathway, group = input, color = input)) +
  geom_linerange(aes(y = pathway, xmin = 0, xmax =NES*-log10(padj)), position = position_dodge(width = 0.2)) +
  geom_point(position = position_dodge(width = 0.2)) +
  theme_minimal() +
  labs(title = 'GSEA', color = 'Input') +
  ylab('Hallmarks') +
  xlab(parse(text = 'NES%*%-log[10]~(P-adj)')) +
  theme(legend.position = 'bottom')

# Saving plot
png('Figures/S6_PRLComparison.png', width = 4000, height = 1500, res = 300)
A + B
dev.off()       
