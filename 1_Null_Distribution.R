library(colorspace)
setwd("~/Documents/Publications/EMTSignature/R")
rm(list = ls())

RNAseq.file.name <- "./data/CCLE_RNAseq_rsem_genes_tpm_20180929.txt"
RNAseq.dat <- read.table(RNAseq.file.name, head=TRUE, sep = "\t", stringsAsFactors = FALSE, na.strings = "null")
GenesIDs <- sapply(RNAseq.dat$gene_id, function(x) unlist(strsplit(as.character(x), "[.]"))[1])
RNAseq.dat$gene_id <- GenesIDs

# Translate Ensemble gene ids to gene_symbols
#tolk <- read.table("./data/GRCh38p13_geneAnnotation.txt", head=TRUE, sep = "\t", fill = TRUE, as.is = TRUE) 
tolk <- read.table("./data/GRCh37_geneAnnotation.txt", head=TRUE, sep = "\t", fill = TRUE, as.is = TRUE) 
tolk_GS <- tolk$Gene.name[match(RNAseq.dat$gene_id, tolk$Gene.stable.ID)]

BC.RNAseq.dat <- data.frame(gene_id = RNAseq.dat[,1], gene_symbol = tolk_GS, RNAseq.dat[,grep("_BREAST", colnames(RNAseq.dat))], stringsAsFactors = FALSE)

# Housekeeping gene scaling
HKGenes <- read.table("./data/HousekeepingGenes-PMID23810203.txt", strip.white = TRUE, head=FALSE, sep = "\t", colClasses = c("character"))

BC.RNAseq.HK <- RNAseq.dat[tolk_GS %in% HKGenes$V1, grep("_BREAST", colnames(RNAseq.dat))]

SCALE.HK <- apply(BC.RNAseq.HK, 2, median)/34.47 # 34.47 average expression of HK genes in CCLE cell lines

# Just normalize gene expression to median housekeeping gene expression and 
# leave in linear expression space
BC_RNAseq_HK <- BC.RNAseq.dat
for (i in 3:ncol(BC_RNAseq_HK))
{
  BC_RNAseq_HK[ ,i] <- BC.RNAseq.dat[ , i]/SCALE.HK[i-2]
}  

plot(SCALE.HK, type = "p")

save(BC_RNAseq_HK, file = "./data/BRCA_CCLE_TPM_HK.rda")
load(file = "./data/BRCA_CCLE_TPM_HK.rda")

EMTgs.file.name <- "./data/EMT_genelist.csv"
tmp <- read.table(EMTgs.file.name, head=TRUE, sep = ",", stringsAsFactors = FALSE, na.strings = "null")
EMTgs <- tmp$EMT_Genes

BC_RNAseq_EMT <- BC_RNAseq_HK[BC_RNAseq_HK$gene_symbol %in% EMTgs, ]

#Convert from TPM to log2(TPM + 0.001) - minimum value not equal to zero is 0.007
SeqRes <- log2(BC_RNAseq_EMT[,c(3:ncol(BC_RNAseq_EMT))] + 0.001)

#####################################################
#Summary results
Ntot = 1000
NCx = vector('double',783*Ntot)
NCy = vector('double',783*Ntot)

for ( i in 1:Ntot )
{
  #Negative control 
  RS <- sample(simplify2array(SeqRes), 783*57, replace=T)
  A2 <- matrix(c(RS), nrow = 783)
  Rdat <- SeqRes
  Rdat[c(1:783),c(1:57)] <- A2
  
  NC_PCA <- prcomp(Rdat, retx = TRUE, scale. = FALSE)
  PCA2.data <- as.matrix(Rdat) %*% NC_PCA$rotation
  NCx[783*(i-1)+c(1:783)] <- PCA2.data[,1]
  NCy[783*(i-1)+c(1:783)] <- PCA2.data[,2]
}

#Scree plot
EigenVals <- NC_PCA$sdev^2
TotVar <- sum(EigenVals)
PrPC <- c(1:10)
for ( i in 1:10 )
{
  PrPC[i] <- EigenVals[i]/TotVar
}

# Results for Eigenvalues
EigenVals[c(1:10)]
TotVar
# [1] 37.03158 35.80556 34.71875 34.25588 34.08679 33.05218 32.53840 32.05880 31.33336 30.96885
# 1348.094
barplot(EigenVals[c(1:10)], names.arg = c(1:10), ylim = c(0,40), ylab = "Eigenvalues", xlab = "Principal Component")

# Results for PrPC
PrPC[c(1:10)]
#  [1] 0.02746959 0.02656014 0.02575396 0.02541061 0.02528518 0.02451772 0.02413660 0.02378084 0.02324272
# [10] 0.02297232
barplot(PrPC[c(1:10)]*100, names.arg = c(1:10), ylim = c(0,40), ylab = "Variance (%)", xlab = "Principal Component")

# null distribution of PC projections
L1 <- density(c(NCx, NCy), adjust = 0.25, na.rm=TRUE, from = -30, to = 30)
SumY <- sum(L1$y)
CumV <- 0
LoLim <- 0
while (CumV < 0.025){
  LoLim <- LoLim + 1
  CumV <-sum(L1$y[1:LoLim])/SumY
}
CumV <- 0
HiLim <- 0
while (CumV < 0.975){
  HiLim <- HiLim + 1
  CumV <-sum(L1$y[1:HiLim])/SumY
}

ValLoLim <- L1$x[LoLim]
ValLoLim
# -12.15
ValHiLim <- L1$x[HiLim]
ValHiLim
# 12.15

plot(L1, main = "Density distribution of PC projections of null distribution")
lines(c(ValLoLim, ValLoLim), c(0,1), lty = 2, col = "red", lwd = 2)
lines(c(ValHiLim, ValHiLim), c(0,1), lty = 2, col = "red", lwd = 2)
text(x = -22, y = 0.05, labels = sprintf("%.2f%%", ValLoLim), cex=1.0)
text(x = 22, y = 0.05, labels = sprintf("%.2f%%", ValHiLim), cex=1.0)


