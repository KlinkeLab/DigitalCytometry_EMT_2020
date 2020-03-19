library(colorspace)
library(MASS)
setwd("~/Documents/Publications/EMTSignature/R")
rm(list = ls())

load(file = "./data/BRCA_CCLE_TPM_HK.rda")

EMTgs.file.name <- "./data/EMT_genelist.csv"
tmp <- read.table(EMTgs.file.name, head=TRUE, sep = ",", stringsAsFactors = FALSE, na.strings = "null")
EMTgs <- tmp$EMT_Genes

BC_RNAseq_EMT <- BC_RNAseq_HK[BC_RNAseq_HK$gene_symbol %in% EMTgs, ]

#Convert from TPM to log2(TPM + 0.001)
SeqRes <- log2(BC_RNAseq_EMT[,c(3:ncol(BC_RNAseq_EMT))] + 0.005)

# Number of Cell lines x Number of genes
# SeqRes gives PCA clustering of genes
# t(SeqRes) gives PCA clustering of cell lines
fs2PCA <- prcomp(SeqRes, retx = TRUE, scale. = FALSE)

# Eigenvals 819.3 172.5 32.4 20.1 12.1 11.0 10.18 8.73 8.05
#Scree plot
EigenVals <- fs2PCA$sdev^2
EigenVals
TotVar <- sum(EigenVals)
TotVar
PrPC <- c(1:10)
for ( i in 1:10 )
{
  PrPC[i] <- EigenVals[i]/TotVar
}

load.Matrix <- data.frame(fs2PCA$rotation)
# n_genes x n_PC          =   n_genes x n_pts      *  n_pts x n_PC
# (square n_genes = n_PC) = (here n_genes < n_pts) * (n_pts > n_PC)
PCA2.data <- as.matrix(SeqRes) %*% fs2PCA$rotation

mp <- barplot(PrPC[c(1:10)]*100, names.arg = c(1:10), ylim = c(0,70), ylab = "Variance (%)", xlab = "Principal Component")
lines(c(0,2,3,4,5,6,7,8,9,12), c(2.75, 2.66, 2.58, 2.54, 2.53, 2.45, 2.41, 2.38, 2.32, 2.30), col="red", lwd = 2, lty = 2)

mp <- barplot(EigenVals[c(1:10)], names.arg = c(1:10), ylim = c(0,900), ylab = "Eigenvalues", xlab = "Principal Component")
lines(c(0,2,3,4,5,6,7,8,9,12), c(37.0, 35.8, 34.7, 34.3, 34.1, 33.1, 32.5, 32.1, 31.3, 31.0), col="red", lwd = 2, lty = 2)

# Keep gene for further analysis if it is expressed above the threshold in more than 
# 5% of the samples
KeepGene <- apply(SeqRes, 1, function(x) length(x[x > 0])/length(x) > 0.05 )

rainColor <- colorRampPalette(c("red", "yellow", "blue"))(100)
meanSeqRes <- rowMeans(SeqRes)
PtCol <- rainColor[cut(meanSeqRes,breaks = 100)]

plot(-PCA2.data[,1], PCA2.data[,2],type = "n", xlim = c(-90, 70), ylim = c(-35, 30))
text(-PCA2.data[,1], PCA2.data[,2], labels = BC_RNAseq_EMT$gene_symbol, cex = 1, col = PtCol)
points(-PCA2.data[BC_RNAseq_EMT$gene_symbol == "CDH1",1], PCA2.data[BC_RNAseq_EMT$gene_symbol == "CDH1",2], pch = 19, col = "blue", cex = 1)
points(-PCA2.data[BC_RNAseq_EMT$gene_symbol == "VIM",1], PCA2.data[BC_RNAseq_EMT$gene_symbol == "VIM",2], pch = 19, col = "red", cex = 1)

plot(PCA2.data[KeepGene,c(3,2)], type = "n", xlim = c(-20, 30), ylim = c(-35, 30))
text(PCA2.data[KeepGene,c(3,2)], labels = BC_RNAseq_EMT$gene_symbol[KeepGene], cex = 0.5, col = "black" ) #PtCol[KeepGene]
points(PCA2.data[BC_RNAseq_EMT$gene_symbol == "CDH1",3], PCA2.data[BC_RNAseq_EMT$gene_symbol == "CDH1",2], pch = 19, col = "blue", cex = 1)
points(PCA2.data[BC_RNAseq_EMT$gene_symbol == "VIM",3], PCA2.data[BC_RNAseq_EMT$gene_symbol == "VIM",2], pch = 19, col = "red", cex = 1)
lines(c(-110, 80), c(12.1, 12.1), col="red", lwd=2, lty=2)
lines(c(-110, 80), c(-12.1, -12.1), col="red", lwd=2, lty=2)
lines(c(-12.1, -12.1), c(-100, 70), col="blue", lwd=2, lty=2)
lines(c(12.1, 12.1), c(-100, 70), col="blue", lwd=2, lty=2)

Esig <- BC_RNAseq_EMT$gene_symbol[PCA2.data[,2] > 12.1 & KeepGene & PCA2.data[,3] < 12.1 & PCA2.data[,3] > -12.1]
Msig <- BC_RNAseq_EMT$gene_symbol[PCA2.data[,2] < -12.1 & KeepGene & PCA2.data[,3] < 12.1 & PCA2.data[,3] > -12.1]

save(Esig, Msig, file = "./data/BRCA_CCLE_Esig_Msig.rda")

##################################################################
#
# EMT signature correlates
#
##################################################################

#Epithelial signature
RNAseq.E <- BC_RNAseq_EMT[match(Esig, BC_RNAseq_EMT$gene_symbol), c(3:ncol(BC_RNAseq_EMT))]

#Mesenchymal signature
RNAseq.M <- BC_RNAseq_EMT[match(Msig, BC_RNAseq_EMT$gene_symbol), c(3:ncol(BC_RNAseq_EMT))]

aRNAseq.E <- log2(RNAseq.E + 0.005)
aRNAseq.M <- log2(RNAseq.M + 0.005)
tmp <- kmeans(t(rbind(aRNAseq.E, aRNAseq.M)), centers = 2)$centers
KDi <- apply(tmp, 2, mean)

KDi.E <- KDi[c(1:length(Esig))]
KDi.M <- KDi[c((length(Esig)+1):length(KDi))]

#Average across genes included in state metric for cell line
EBrCa <- apply(RNAseq.E/(RNAseq.E + 2^KDi.E), 2, mean)
MBrCa <- apply(RNAseq.M/(RNAseq.M + 2^KDi.M), 2, mean)
write.csv(cbind(c(rep("Epi", length(Esig)), rep("Mes", length(Msig))), c(Esig, Msig), KDi), file = "GenePCAresult_TPM2.csv")

##################################################################
#
##################################################################
BrCaTypes.file.name <- "./data/BrCaCellLines-IntrinsicSubtypes.csv"
tmp <- read.csv(file = BrCaTypes.file.name)
BrCa_Types <- data.frame(CellName = tmp$X, Type = tmp$BrCaType)

BrCa_EMT <- BrCa_Types[match(colnames(RNAseq.E), BrCa_Types$CellName), ]

# Color based on sample type
FillColors <- c("red", "yellow", "pink", "skyblue", "black")
OutColors <- c("red", "black", "red", "blue", "black")
ColorFill <- FillColors[as.factor(BrCa_EMT$Type)]
ColorOut <- OutColors[as.factor(BrCa_EMT$Type)]
SlimCL <- lapply(as.character(BrCa_EMT$CellName), function(x) substr(x, 1, nchar(x)-7))

xvals <- 10^seq(-5,0,len = 1000)
yvals <- 1 - xvals

plot(MBrCa, EBrCa, type = "p", pch = 21, col = ColorOut, bg = ColorFill, 
     xlab = "Mesenchymal Signal", ylab = "Epithelial Signal", xlim = c(0, 1), ylim = c(0,1))
lines(xvals, yvals, lty = 2, col = "black")
text(MBrCa, EBrCa+0.05, labels = SlimCL, cex = 0.5)
legend(x="bottomleft", legend = c("Basal", "Claudin Low", "HER2", "Luminal A", "Luminal B"), col = c("red", "black", "red", "blue", "black"), pt.bg = c("red", "yellow", "pink", "skyblue", "black"), pch = 21, cex = 1)


# Plot VIM versus E-cadherin
CDH1 <- as.numeric(BC_RNAseq_EMT[BC_RNAseq_EMT$gene_symbol == "CDH1",c(3:ncol(BC_RNAseq_EMT))])
VIM <- as.numeric(BC_RNAseq_EMT[BC_RNAseq_EMT$gene_symbol == "VIM",c(3:ncol(BC_RNAseq_EMT))])

plot(log2(VIM + 0.001), log2(CDH1+ 0.001), pch = 21, col = ColorOut, bg = ColorFill, xlab = "Vimentin (TPM)", ylab = "E-cadherin (TPM)")#, xlim = c(-4, 0), ylim = c(-4.5,0))
xvals <- 10^seq(-5,0,len = 1000)
yvals <- 1 - xvals
lines(log2(xvals+0.0001), log2(yvals+0.0001), lty = 2, col = "black")
