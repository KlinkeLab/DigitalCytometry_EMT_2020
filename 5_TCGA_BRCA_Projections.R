library(colorspace)
rm(list = ls())

# data: https://github.com/KlinkeLab/DigitalCytometry_EMT_2020/blob/master/Files/TCGA-BRCA.GDC_phenotype.tsv
TCGA_Phenotypes <- read.table("TCGA-BRCA.GDC_phenotype.tsv", head=TRUE, sep = "\t", fill = TRUE, colClasses = c("character"))

#Load variable GEData_HK
# data: generated by 3_TCGA_BRCA_PreProcess.R
load(file = "BRCA_TCGA_TPM_HK.rda")

##################################################################
#
# EMT signature correlates
#
##################################################################

#Epithelial and Mesenchymal state metric genes
# data: generated by 4_TCGA_BRCA_Ridge.R
load(file = "BRCA_TCGA_Esig_Msig.rda")

RNAseq.E <- GEData_HK[match(E_TGenes, row.names(GEData_HK)), ]

E.thresh <- kmeans(t(log2(RNAseq.E + 0.03)), centers = 2)$centers
E.thresh2 <- apply(E.thresh, 2, mean)

EBrCa <- rep(0, dim(RNAseq.E)[2])
for (i in 1:dim(RNAseq.E)[2])
{
  EBrCa[i] <- mean(RNAseq.E[,i]/(RNAseq.E[,i] + 2^E.thresh2), na.rm = TRUE)
}

# Mesenchymal signature
# change C7orf10:SUGCT, LEPRE1:P3H1, LHFP:LHFPL6

RNAseq.M <- GEData_HK[match(M_TGenes, row.names(GEData_HK)), ]

M.thresh <- kmeans(t(log2(RNAseq.M + 0.03)), centers = 2)$centers
M.thresh2 <- apply(M.thresh, 2, mean)

MBrCa <- rep(0, dim(RNAseq.M)[2])
for (i in 1:dim(RNAseq.M)[2])
{
  MBrCa[i] <- mean(RNAseq.M[,i]/(RNAseq.M[,i] + 2^M.thresh2), na.rm = TRUE)
}

#Combine Epithelial and Mesenchymal genes for each sample
BrCa_State <- data.frame(TCGAName = colnames(RNAseq.M), Epithelial = EBrCa, Mesenchymal = MBrCa)#, Intrinsic_Subtype = BCType)
write.csv(BrCa_State, file = "BrCa-TCGA-EMT-State-Mar20.csv")

BrCa_EMT <- read.csv(file = "BrCa-TCGA-EMT-State-Mar20.csv")
##################################################################
#
#
#
##################################################################
# Identify patient subtype
# Primary: 1, Normal: 11, Metastatic: 6
SampleIDs <- substr(colnames(RNAseq.M), 1,16)
TCGA_Pheno_idx <- match(SampleIDs, TCGA_Phenotypes$submitter_id.samples)
ER_status <- TCGA_Phenotypes$breast_carcinoma_estrogen_receptor_status[TCGA_Pheno_idx] 
PR_status <- TCGA_Phenotypes$breast_carcinoma_progesterone_receptor_status[TCGA_Pheno_idx]
HER2_IHC_status <- TCGA_Phenotypes$lab_proc_her2_neu_immunohistochemistry_receptor_status[TCGA_Pheno_idx]
HER2_ISH_status <- TCGA_Phenotypes$lab_procedure_her2_neu_in_situ_hybrid_outcome_type[TCGA_Pheno_idx]
Tissue_type <- TCGA_Phenotypes$sample_type.samples[TCGA_Pheno_idx]
BrCa_Type = rep("TN", length(SampleIDs))
for (i in 1:length(BrCa_Type))
{
  BrCa_Type[i] <- ifelse(ER_status[i] == "Positive" | PR_status[i] == "Positive", "ER/PR+", "TN")
  BrCa_Type[i] <- ifelse(HER2_IHC_status[i] == "Positive" | HER2_ISH_status[i] == "Positive", "HER2+", BrCa_Type[i])
  BrCa_Type[i] <- ifelse(Tissue_type[i] == "Solid Tissue Normal", "Normal", BrCa_Type[i])
}

DotColors <- c("blue", "pink", "green", "red")
ColorDot <- DotColors[as.factor(BrCa_Type)]

xvals <- 10^seq(-5,0,len = 1000)
yvals <- 1 - xvals

plot(BrCa_EMT$Mesenchymal, BrCa_EMT$Epithelial, type = "p", pch = 19, cex = 0.5, col = ColorDot, 
     xlab = "Mesenchymal Signal", ylab = "Epithelial Signal", xlim = c(0, 1), ylim = c(0,1))
lines(xvals, yvals, lty = 2, col = "black")
legend(x="bottomleft", legend = c("Normal", "ER/PR", "HER2", "TN"), col = c("green", "blue", "pink", "red"), pch = 19, cex = 1)

##### Make patient subtype density distribution
Nkeep <- BrCa_Type == "Normal"
TNkeep <- BrCa_Type == "TN"
ERkeep <- BrCa_Type == "ER/PR+"
Hkeep <- BrCa_Type == "HER2+"

Nden = density(BrCa_EMT$Epithelial[Nkeep], adj = 0.4, from = 0, to = 1)
TNden = density(BrCa_EMT$Epithelial[TNkeep], adj = 0.4, from = 0, to = 1)
ERden = density(BrCa_EMT$Epithelial[ERkeep], adj = 0.4, from = 0, to = 1)
Hden = density(BrCa_EMT$Epithelial[Hkeep], adj = 0.4, from = 0, to = 1)

xlimit = c(1, 0)
plot(Nden$x, Nden$y, type = "l", xlim = xlimit, ylim = c(0,11), xlab = "TPM", ylab = "Density", col = "green", lwd = 2)
lines(TNden$x, TNden$y, col = "red", lwd = 2)
lines(ERden$x, ERden$y, col = "blue", lwd = 2)
lines(Hden$x, Hden$y, col = "black", lwd = 2)
lines(Hden$x, Hden$y, col = "pink", lwd = 2, lty = 2)

median(BrCa_EMT$Epithelial[Hkeep])
median(BrCa_EMT$Epithelial[TNkeep])
median(BrCa_EMT$Epithelial[ERkeep])

Nden = density(BrCa_EMT$Mesenchymal[Nkeep], adj = 0.4, from = 0, to = 1)
TNden = density(BrCa_EMT$Mesenchymal[TNkeep], adj = 0.4, from = 0, to = 1)
ERden = density(BrCa_EMT$Mesenchymal[ERkeep], adj = 0.4, from = 0, to = 1)
Hden = density(BrCa_EMT$Mesenchymal[Hkeep], adj = 0.4, from = 0, to = 1)

xlimit = c(0, 1)
plot(Nden$x, Nden$y, type = "l", xlim = xlimit, ylim = c(0,10), xlab = "TPM", ylab = "Density", col = "green", lwd = 2)
lines(TNden$x, TNden$y, col = "red", lwd = 2)
lines(ERden$x, ERden$y, col = "blue", lwd = 2)
lines(Hden$x, Hden$y, col = "black", lwd = 2)
lines(Hden$x, Hden$y, col = "pink", lwd = 2, lty = 2)

median(BrCa_EMT$Mesenchymal[Hkeep])
median(BrCa_EMT$Mesenchymal[TNkeep])
median(BrCa_EMT$Mesenchymal[ERkeep])

##################################################################################################
#
#
# End of script
#
#
##################################################################################################
