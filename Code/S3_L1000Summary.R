# BiocManager::install('ccmap')
# BiocManager::install('ccdata')

library(ccdata)
library(ccmap)
library(pbapply)
library(ggplot2)
library(statsExpressions)
library(ggrepel)
library(patchwork)
data("l1000_es")

corMin <- 0.6
bCL <- c('BT20', 'HS578T', 'MDAMB231')

l1000_names <- colnames(l1000_es)
l1000_compounds <- unlist(lapply(strsplit(l1000_names, '_'), function(X){paste0(X[1:(length(X)-3)], collapse = '_')}))
l1000_cellLine <- unlist(lapply(strsplit(l1000_names, '_'), function(X){X[(length(X)-2)]}))
l1000_concentration <- unlist(lapply(strsplit(l1000_names, '_'), function(X){X[(length(X)-1)]}))
l1000_time <- unlist(lapply(strsplit(l1000_names, '_'), function(X){X[length(X)]}))
sList <- cbind(l1000_names)

S1 <- data.frame(compound = l1000_compounds, cellLine = l1000_cellLine, concentration = l1000_concentration)
S1 <- S1[S1$cellLine %in% bCL,]
nrow(S1)
# Profiles
S1 <- unique(S1)
# Compounds
length(unique(S1$compound))
# Concentrations
table(S1$concentration)

# 4899 profiles 205 compounds, 4 concentrations, 2 time points in 3 cell TNBC cell lines
S1 <- apply(S1, 1, function(X){paste0(X, collapse = '_')})

S1Profiles <- pbsapply(S1, function(N){
  X <- l1000_es[,l1000_names[grepl(N, l1000_names)], drop = FALSE]
  if(ncol(X)>1){
    X <- preprocessCore::normalize.quantiles(as.matrix(X))
  }
  corValue <- cor(data.frame(rowMeans(X),X),method = 'sp')[,1]
  corValue <- (corValue[-1] > corMin)
  if(sum(corValue) > 1){
    return(rowMeans(X[,corValue]))
  } else {
    return(rep(NA, 1001))
  }
})
colnames(S1Profiles) <- S1
S1Profiles <- S1Profiles[,complete.cases(t(S1Profiles))]  
S1_names <- colnames(S1Profiles)
rownames(S1Profiles) <- rownames(l1000_es)
# dim(S1Profiles) : 1001 2111
write.csv(S1Profiles, '../Results/S1_CellLinesConcentrationProfiles.csv')

S2 <- unique(data.frame(compound = l1000_compounds, cellLine = l1000_cellLine))
S2 <- S2[S2$cellLine %in% bCL,]
S2 <- apply(S2, 1, function(X){paste0(X, collapse = '_')})

S2Profiles <- pbsapply(S2, function(N){
  X <- S1Profiles[,grepl(N, S1_names, fixed = TRUE), drop = FALSE]
  if(ncol(X)>1){
    X <- preprocessCore::normalize.quantiles(as.matrix(X))
  }
  corValue <- cor(data.frame(rowMeans(X),X),method = 'sp')[,1]
  corValue <- (corValue[-1] > corMin)
  if(sum(corValue) > 1){
    return(rowMeans(X[,corValue]))
  } else {
    return(rep(NA, 1001))
  }
})
colnames(S2Profiles) <- S2
S2Profiles <- S2Profiles[,complete.cases(t(S2Profiles))]  
S2_names <- colnames(S2Profiles)
rownames(S2Profiles) <- rownames(l1000_es)
# dim(S2Profiles): 1001  528
write.csv(S2Profiles, '../Results/S1_CellLinesProfiles.csv')


S3 <- unique(unlist(lapply(strsplit(S2_names, '_'), function(X){X[1]})))
S3Profiles <- pbsapply(S3, function(N){
  X <- S2Profiles[,grepl(N, S2_names, fixed = TRUE), drop = FALSE]
  if(ncol(X)>1){
    X <- preprocessCore::normalize.quantiles(as.matrix(X))
  }
  corValue <- cor(data.frame(rowMeans(X),X),method = 'sp')[,1]
  corValue <- (corValue[-1] > corMin)
  if(sum(corValue) > 1){
    return(rowMeans(X[,corValue]))
  } else {
    return(rep(NA, 1001))
  }
})
colnames(S3Profiles) <- S3
S3Profiles <- S3Profiles[,complete.cases(t(S3Profiles))]  
S3_names <- colnames(S3Profiles)
rownames(S3Profiles) <- rownames(l1000_es)
# dim(S3Profiles): 1001  176
write.csv(S3Profiles, '../Results/S1_Profiles.csv')


# QL-XII-47 Example
T1 <- l1000_es[,colnames(l1000_es)[grepl('QL-XII-47_MDAMB231_0.08um',colnames(l1000_es))]]
dnT1 <- dimnames(T1)
T1 <- preprocessCore::normalize.quantiles(T1)
dimnames(T1) <- dnT1
T1 <- as.data.frame(T1)[,c(2,1)]
T1[,3] <- rowMeans(T1)
colnames(T1) <- c('X6h', 'X24h', 'Avg')
T1$G <- rownames(T1)
T1 <- T1[order(abs(T1$Avg), decreasing = TRUE),]
T1$G[11:nrow(T1)] <- NA

corValue <- cor(T1$X24h,T1$X6h, method = 'sp')
A1 <- ggplot(T1, aes(X24h,X6h, label = G)) + 
  geom_point(alpha = 0.3, pch = 16) + 
  geom_density2d() + 
  theme_bw() + 
  geom_abline(intercept = 0, slope = ifelse(corValue > 0, 1,-1), lty = 2, col = 'red') +
  geom_text_repel(fontface=3, min.segment.length = 0, bg.color = 'white', bg.r = 0.1) + 
  labs(tag = 'A', title = parse(text = 'QL-XII-47 - MDAMB231 - 0.08~mu*M'), subtitle = parse(text = paste0('rho == ', round(corValue,3)))) +
  xlab('24 h') +
  ylab('6 h') + 
  theme(panel.border = element_rect(colour = ifelse(corValue >= 0.6, "black", "gray70"), fill=NA, size=2), plot.title = element_text(face = 'bold'))

corValue <- cor(T1$Avg,T1$X6h, method = 'sp')
A2 <- ggplot(T1, aes(Avg,X6h, label = G)) + 
  geom_point(alpha = 0.3, pch = 16) + 
  geom_density2d() + 
  theme_bw() + 
  geom_abline(intercept = 0, slope = ifelse(corValue > 0, 1,-1), lty = 2, col = 'red') +
  geom_text_repel(fontface=3, size=3, min.segment.length = 0, bg.color = 'white', bg.r = 0.1) + 
  labs(subtitle = parse(text = paste0('rho == ', round(corValue,3)))) +
  xlab(parse(text = 'Average~(QL-XII-47~0.08~mu*M~Samples-MDAMB231)')) +
  ylab('6 h') + 
  theme(panel.border = element_rect(colour = ifelse(corValue >= 0.6,  "black", "gray70"), fill=NA, size=2), plot.title = element_text(face = 2))

corValue <- cor(T1$Avg,T1$X24h, method = 'sp')
A3 <- ggplot(T1, aes(Avg,X24h , label = G)) + 
  geom_point(alpha = 0.3, pch = 16) + 
  geom_density2d() + 
  theme_bw() + 
  geom_abline(intercept = 0, slope = ifelse(corValue > 0, 1,-1), lty = 2, col = 'red') +
  geom_text_repel(fontface=3, size=3, min.segment.length = 0, bg.color = 'white', bg.r = 0.1) + 
  labs(subtitle = parse(text = paste0('rho == ', round(corValue,3)))) +
  xlab(parse(text = 'Average~(QL-XII-47~0.08~mu*M~Samples-MDAMB231)')) +
  ylab('24 h') + 
  theme(panel.border = element_rect(colour = ifelse(corValue >= 0.6, "black", "gray70"), fill=NA, size=2), plot.title = element_text(face = 2))

T1 <- l1000_es[,colnames(l1000_es)[grepl('QL-XII-47_MDAMB231_0.4um',colnames(l1000_es))]]
dnT1 <- dimnames(T1)
T1 <- preprocessCore::normalize.quantiles(T1)
dimnames(T1) <- dnT1
T1 <- as.data.frame(T1)[,c(2,1)]
T1[,3] <- rowMeans(T1)
colnames(T1) <- c('X6h', 'X24h', 'Avg')
T1$G <- rownames(T1)
T1 <- T1[order(abs(T1$Avg), decreasing = TRUE),]
T1$G[11:nrow(T1)] <- NA
corValue <- cor(T1$X24h,T1$X6h, method = 'sp')
A4 <- ggplot(T1, aes(X24h,X6h, label = G)) + 
  geom_point(alpha = 0.3, pch = 16) + 
  geom_density2d() + 
  theme_bw() + 
  geom_abline(intercept = 0, slope = ifelse(corValue > 0, 1,-1), lty = 2, col = 'red') +
  geom_text_repel(fontface=3, min.segment.length = 0, bg.color = 'white', bg.r = 0.1) + 
  labs(title = parse(text = 'QL-XII-47 - MDAMB231 - 0.4~mu*M'), subtitle = parse(text = paste0('rho == ', round(corValue,3)))) +
  xlab('24 h') +
  ylab('6 h') + 
  theme(panel.border = element_rect(colour = ifelse(corValue >= 0.6, "black", "gray70"), fill=NA, size=2), plot.title = element_text(face = 2))

corValue <- cor(T1$Avg,T1$X6h, method = 'sp')
A5 <- ggplot(T1, aes(Avg,X6h, label = G)) + 
  geom_point(alpha = 0.3, pch = 16) + 
  geom_density2d() + 
  theme_bw() + 
  geom_abline(intercept = 0, slope = ifelse(corValue > 0, 1,-1), lty = 2, col = 'red') +
  geom_text_repel(fontface=3, size=3, min.segment.length = 0, bg.color = 'white', bg.r = 0.1) + 
  labs(subtitle = parse(text = paste0('rho == ', round(corValue,3)))) +
  xlab(parse(text = 'Average~(QL-XII-47~0.4~mu*M~Samples-MDAMB231)')) +
  ylab('6 h') + 
  theme(panel.border = element_rect(colour = ifelse(corValue >= 0.6, "black", "gray70"), fill=NA, size=2), plot.title = element_text(face = 2))

corValue <- cor(T1$Avg,T1$X24h, method = 'sp')
A6 <- ggplot(T1, aes(Avg,X24h , label = G)) + 
  geom_point(alpha = 0.3, pch = 16) + 
  geom_density2d() + 
  theme_bw() + 
  geom_abline(intercept = 0, slope = ifelse(corValue > 0, 1,-1), lty = 2, col = 'red') +
  geom_text_repel(fontface=3, size=3, min.segment.length = 0, bg.color = 'white', bg.r = 0.1) + 
  labs(subtitle = parse(text = paste0('rho == ', round(corValue,3)))) +
  xlab(parse(text = 'Average~(QL-XII-47~0.4~mu*M~Samples-MDAMB231)')) +
  ylab('24 h') + 
  theme(panel.border = element_rect(colour = ifelse(corValue >= 0.6, "black", "gray70"), fill=NA, size=2), plot.title = element_text(face = 2))


T1 <- l1000_es[,colnames(l1000_es)[grepl('QL-XII-47_MDAMB231_2um',colnames(l1000_es))]]
dnT1 <- dimnames(T1)
T1 <- preprocessCore::normalize.quantiles(T1)
dimnames(T1) <- dnT1
T1 <- as.data.frame(T1)[,c(2,1)]
T1[,3] <- rowMeans(T1)
colnames(T1) <- c('X6h', 'X24h', 'Avg')
T1$G <- rownames(T1)
T1 <- T1[order(abs(T1$Avg), decreasing = TRUE),]
T1$G[11:nrow(T1)] <- NA

corValue <- cor(T1$X24h,T1$X6h, method = 'sp')
A7 <- ggplot(T1, aes(X24h,X6h, label = G)) + 
  geom_point(alpha = 0.3, pch = 16) + 
  geom_density2d() + 
  theme_bw() + 
  geom_abline(intercept = 0, slope = ifelse(corValue > 0, 1,-1), lty = 2, col = 'red') +
  geom_text_repel(fontface=3, min.segment.length = 0, bg.color = 'white', bg.r = 0.1) + 
  labs(title = parse(text = 'QL-XII-47 - MDAMB231 - 2~mu*M'), subtitle = parse(text = paste0('rho == ', round(corValue,3)))) +
  xlab('24 h') +
  ylab('6 h') + 
  theme(plot.title = element_text(face = 2), panel.border = element_rect(colour = ifelse(corValue >= 0.6, "black", "gray70"), fill=NA, size=2))

corValue <- cor(T1$Avg,T1$X6h, method = 'sp')
A8 <- ggplot(T1, aes(Avg,X6h, label = G)) + 
  geom_point(alpha = 0.3, pch = 16) + 
  geom_density2d() + 
  theme_bw() + 
  geom_abline(intercept = 0, slope = ifelse(corValue > 0, 1,-1), lty = 2, col = 'red') +
  geom_text_repel(fontface=3, size=3, min.segment.length = 0, bg.color = 'white', bg.r = 0.1) + 
  labs(subtitle = parse(text = paste0('rho == ', round(corValue,3)))) +
  xlab(parse(text = 'Average~(QL-XII-47~2~mu*M~Samples-MDAMB231)')) +
  ylab('6 h') + 
  theme(panel.border = element_rect(colour = ifelse(corValue >= 0.6, "black", "gray70"), fill=NA, size=2), plot.title = element_text(face = 2))

corValue <- cor(T1$Avg,T1$X24h, method = 'sp')
A9 <- ggplot(T1, aes(Avg,X24h , label = G)) + 
  geom_point(alpha = 0.3, pch = 16) + 
  geom_density2d() + 
  theme_bw() + 
  geom_abline(intercept = 0, slope = ifelse(corValue > 0, 1,-1), lty = 2, col = 'red') +
  geom_text_repel(fontface=3, size=3, min.segment.length = 0, bg.color = 'white', bg.r = 0.1) + 
  labs(subtitle = parse(text = paste0('rho == ', round(corValue,3)))) +
  xlab(parse(text = 'Average~(QL-XII-47~2~mu*M~Samples-MDAMB231)')) +
  ylab('24 h') + 
  theme(panel.border = element_rect(colour = ifelse(corValue >= 0.6, "black", "gray70"), fill=NA, size=2), plot.title = element_text(face = 2))

T1 <- l1000_es[,colnames(l1000_es)[grepl('QL-XII-47_MDAMB231_10um',colnames(l1000_es))]]
dnT1 <- dimnames(T1)
T1 <- preprocessCore::normalize.quantiles(T1)
dimnames(T1) <- dnT1
T1 <- as.data.frame(T1)[,c(2,1)]
T1[,3] <- rowMeans(T1)
colnames(T1) <- c('X6h', 'X24h', 'Avg')
T1$G <- rownames(T1)
T1 <- T1[order(abs(T1$Avg), decreasing = TRUE),]
T1$G[11:nrow(T1)] <- NA
corValue <- cor(T1$X24h,T1$X6h, method = 'sp')
A10 <- ggplot(T1, aes(X24h,X6h, label = G)) + 
  geom_point(alpha = 0.3, pch = 16) + 
  geom_density2d() + 
  theme_bw() + 
  geom_abline(intercept = 0, slope = ifelse(corValue > 0, 1,-1), lty = 2, col = 'red') +
  geom_text_repel(fontface=3, min.segment.length = 0, bg.color = 'white', bg.r = 0.1) + 
  labs(title = parse(text = 'QL-XII-47 - MDAMB231 - 10~mu*M'), subtitle = parse(text = paste0('rho == ', round(corValue,3)))) +
  xlab('24 h') +
  ylab('6 h') + 
  theme(panel.border = element_rect(colour = ifelse(corValue >= 0.6, "black", "gray70"), fill=NA, size=2), plot.title = element_text(face = 2))

corValue <- cor(T1$Avg,T1$X6h, method = 'sp')
A11 <- ggplot(T1, aes(Avg,X6h, label = G)) + 
  geom_point(alpha = 0.3, pch = 16) + 
  geom_density2d() + 
  theme_bw() + 
  geom_abline(intercept = 0, slope = ifelse(corValue > 0, 1,-1), lty = 2, col = 'red') +
  geom_text_repel(fontface=3, size=3, min.segment.length = 0, bg.color = 'white', bg.r = 0.1) + 
  labs(subtitle = parse(text = paste0('rho == ', round(corValue,3)))) +
  xlab(parse(text = 'Average~(QL-XII-47~10~mu*M~Samples-MDAMB231)')) +
  ylab('6 h') + 
  theme(panel.border = element_rect(colour = ifelse(corValue >= 0.6, "black", "gray70"), fill=NA, size=2), plot.title = element_text(face = 2))

corValue <- cor(T1$Avg,T1$X24h, method = 'sp')
A12 <- ggplot(T1, aes(Avg,X24h , label = G)) + 
  geom_point(alpha = 0.3, pch = 16) + 
  geom_density2d() + 
  theme_bw() + 
  geom_abline(intercept = 0, slope = ifelse(corValue > 0, 1,-1), lty = 2, col = 'red') +
  geom_text_repel(fontface=3, size=3, min.segment.length = 0, bg.color = 'white', bg.r = 0.1) + 
  labs(subtitle = parse(text = paste0('rho == ', round(corValue,3)))) +
  xlab(parse(text = 'Average~(QL-XII-47~10~mu*M~Samples-MDAMB231)')) +
  ylab('24 h') + 
  theme(panel.border = element_rect(colour = ifelse(corValue >= 0.6, "black", "gray70"), fill=NA, size=2), plot.title = element_text(face = 2))


# png('../Figures/F1A.png', width = 4800, height = 2500, res = 300)
# (A1|(A2/A3))/(A7|(A8/A9))|(A4|(A5/A6))/(A10|(A11/A12))
# dev.off()

# QL-XII-47 Example
T2 <- S1Profiles[,colnames(S1Profiles)[grepl('QL-XII-47_MDA',colnames(S1Profiles))]]
dnT2 <- dimnames(T2)
T2 <- preprocessCore::normalize.quantiles(as.matrix(T2))
dimnames(T2) <- dnT2
T2 <- as.data.frame(T2)
T2$Avg <- rowMeans(T2)
T2 <- T2[order(abs(T2$Avg), decreasing = TRUE),]
T2$G <- rownames(T2)
T2$G[11:nrow(T2)] <- NA
colnames(T2) <- c('X0.08um', 'X0.4um', 'X10um', 'X2um', 'Avg', 'G')

corValue <- cor(T2$Avg, T2$X0.08um, method = 'sp')
B1 <- ggplot(T2, aes(Avg, X0.08um, label = G)) + 
  geom_point(alpha = 0.3, pch = 16) + 
  geom_density2d() + 
  theme_bw() + 
  geom_abline(intercept = 0, slope = ifelse(corValue > 0, 1,-1), lty = 2, col = 'red') +
  geom_text_repel(fontface=3, min.segment.length = 0, bg.color = 'white', bg.r = 0.1) + 
  labs(tag = 'B',subtitle = parse(text = paste0('rho == ', round(corValue,3)))) +
  xlab(parse(text = 'Average~(QL-XII-47~Samples-MDAMB231)')) +
  ylab(expression(atop('Average',('QL-XII-47'~0.08~mu*'M'~'Samples-MDAMB231'))))+
  theme(panel.border = element_rect(colour = ifelse(corValue >= 0.6, "black", "gray70"), fill=NA, size=2), plot.title = element_text(face = 2))

corValue <- cor(T2$Avg, T2$X0.4um, method = 'sp')
B2 <- ggplot(T2, aes(Avg, X0.4um, label = G)) + 
  geom_point(alpha = 0.3, pch = 16) + 
  geom_density2d() + 
  theme_bw() + 
  geom_abline(intercept = 0, slope = ifelse(corValue > 0, 1,-1), lty = 2, col = 'red') +
  geom_text_repel(fontface=3, min.segment.length = 0, bg.color = 'white', bg.r = 0.1) + 
  labs(subtitle = parse(text = paste0('rho == ', round(corValue,3)))) +
  xlab(parse(text = 'Average~(QL-XII-47~Samples-MDAMB231)')) +
  ylab(expression(atop('Average',('QL-XII-47'~0.4~mu*'M'~'Samples-MDAMB231'))))+
  theme(panel.border = element_rect(colour = ifelse(corValue >= 0.6, "black", "gray70"), fill=NA, size=2), plot.title = element_text(face = 2))

corValue <- cor(T2$Avg, T2$X2um, method = 'sp')
B3 <- ggplot(T2, aes(Avg, X2um, label = G)) + 
  geom_point(alpha = 0.3, pch = 16) + 
  geom_density2d() + 
  theme_bw() + 
  geom_abline(intercept = 0, slope = ifelse(corValue > 0, 1,-1), lty = 2, col = 'red') +
  geom_text_repel(fontface=3, min.segment.length = 0, bg.color = 'white', bg.r = 0.1) + 
  labs(subtitle = parse(text = paste0('rho == ', round(corValue,3)))) +
  xlab(parse(text = 'Average~(QL-XII-47~Samples- MDAMB231)')) +
  ylab(expression(atop('Average',('QL-XII-47'~2~mu*'M'~'Samples-MDAMB231'))))+
  theme(panel.border = element_rect(colour = ifelse(corValue >= 0.6, "black", "gray70"), fill=NA, size=2), plot.title = element_text(face = 2))

corValue <- cor(T2$Avg, T2$X10um, method = 'sp')
B4 <- ggplot(T2, aes(Avg, X10um, label = G)) + 
  geom_point(alpha = 0.3, pch = 16) + 
  geom_density2d() + 
  theme_bw() + 
  geom_abline(intercept = 0, slope = ifelse(corValue > 0, 1,-1), lty = 2, col = 'red') +
  geom_text_repel(fontface=3, min.segment.length = 0, bg.color = 'white', bg.r = 0.1) + 
  labs(subtitle = parse(text = paste0('rho == ', round(corValue,3)))) +
  xlab(parse(text = 'Average~(QL-XII-47~Samples-MDAMB231)')) +
  ylab(expression(atop('Average',('QL-XII-47'~10~mu*'M'~'Samples-MDAMB231'))))+
  theme(panel.border = element_rect(colour = ifelse(corValue >= 0.6, "black", "gray70"), fill=NA, size=2), plot.title = element_text(face = 2))

# png('../Figures/F1B.png', width = 4800, height = 1300, res = 300)
# B1 | B2 | B3 | B4
# dev.off()

T3 <- S2Profiles[,S2_names[grepl('QL-XII-47',S2_names)]]
dnT3 <- dimnames(T3)
T3 <- preprocessCore::normalize.quantiles(as.matrix(T3))
T3 <- as.data.frame(T3)
dimnames(T3) <- dnT3
T3$Avg <- rowMeans(T3)
T3 <- T3[order(abs(T3$Avg), decreasing = TRUE),]
T3$G <- rownames(T3)
T3$G[11:nrow(T3)] <- NA
colnames(T3) <- c('BT20', 'HS578T', 'MDAMB231', 'Avg', 'G')

corValue <- cor(T3$BT20, T3$Avg, method = 'sp')
C1 <- ggplot(T3, aes(Avg, BT20, label = G)) + 
  geom_point(alpha = 0.3, pch = 16) + 
  geom_density2d() + 
  theme_bw() + 
  geom_abline(intercept = 0, slope = ifelse(corValue > 0, 1,-1), lty = 2, col = 'red') +
  geom_text_repel(fontface=3, min.segment.length = 0, bg.color = 'white', bg.r = 0.1) + 
  labs(tag = 'C', title = 'QL-XII-47 - BT20', subtitle = parse(text = paste0('rho == ', round(corValue,3)))) +
  xlab(parse(text = 'Average~(QL-XII-47~Samples~Across~Cell~Lines)')) +
  ylab(parse(text = 'Average~(QL-XII-47 - BT20)')) + 
  theme(panel.border = element_rect(colour = ifelse(corValue >= 0.6, "black", "gray70"), fill=NA, size=2), plot.title = element_text(face = 1))

corValue <- cor(T3$HS578T, T3$Avg, method = 'sp')
C2 <- ggplot(T3, aes(Avg, HS578T, label = G)) + 
  geom_point(alpha = 0.3, pch = 16) + 
  geom_density2d() + 
  theme_bw() + 
  geom_abline(intercept = 0, slope = ifelse(corValue > 0, 1,-1), lty = 2, col = 'red') +
  geom_text_repel(fontface=3, min.segment.length = 0, bg.color = 'white', bg.r = 0.1) + 
  labs(title = 'QL-XII-47 - HS578T', subtitle = parse(text = paste0('rho == ', round(corValue,3)))) +
  xlab(parse(text = 'Average~(QL-XII-47~Samples~Across~Cell~Lines)')) +
  ylab(parse(text = 'Average~(QL-XII-47 - HS578T)')) + 
  theme(panel.border = element_rect(colour = ifelse(corValue >= 0.6, "black", "gray70"), fill=NA, size=2), plot.title = element_text(face = 1))

corValue <- cor(T3$MDAMB231, T3$Avg, method = 'sp')
C3 <- ggplot(T3, aes(Avg, MDAMB231, label = G)) + 
  geom_point(alpha = 0.3, pch = 16) + 
  geom_density2d() + 
  theme_bw() + 
  geom_abline(intercept = 0, slope = ifelse(corValue > 0, 1,-1), lty = 2, col = 'red') +
  geom_text_repel(fontface=3, min.segment.length = 0, bg.color = 'white', bg.r = 0.1) + 
  labs(title = 'QL-XII-47 - MDAMB231', subtitle = parse(text = paste0('hat(rho) == ', round(corValue,3)))) +
  xlab(parse(text = 'Average~(QL-XII-47~Samples~Across~Cell~Lines)')) +
  ylab(parse(text = 'Average~(QL-XII-47 - MDAMB231)')) + 
  theme(panel.border = element_rect(colour = ifelse(corValue >= 0.6, "black", "gray70"), fill=NA, size=2), plot.title = element_text(face = 1))

# png('../Figures/F1C.png', width = 4800, height = 1500, res = 300)
# C1 | C2 | C3
# dev.off()

layout <- "
AAABBBDDDEEE
AAACCCDDDFFF
GGGHHHJJJKKK
GGGIIIJJJLLL
MMMNNNOOOPPP
MMMNNNOOOPPP
QQQQRRRRSSSS
QQQQRRRRSSSS"

png('../Figures/F3.png', width = 5600, height = 5000, res = 300)
A1 + A2 + A3 + A4 + A5 + A6 + A7 + A8 + A9 + A10 + A11 + 
  A12 + B1 + B2 + B3 + B4 + C1 + C2 + C3 + patchwork::plot_layout(design = layout, guides = 'collect', tag_level = 'new')
dev.off()
