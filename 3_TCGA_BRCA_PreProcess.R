library(SummarizedExperiment)
library(TCGAbiolinks)

rm(list = ls())

query.seq <- GDCquery(project = "TCGA-BRCA", 
                       data.category = "Transcriptome Profiling", 
                       data.type = "Gene Expression Quantification",
                       sample.type = c("Primary solid Tumor","Solid Tissue Normal"),
                       workflow.type = "HTSeq - Counts")

GDCdownload(query.seq, files.per.chunk = 4)
seq.brca <- GDCprepare(query = query.seq, save = TRUE, save.filename = "brca-RNAseq-Counts.rda", summarizedExperiment = TRUE)

# data: https://github.com/KlinkeLab/DigitalCytometry_EMT_2020/blob/master/Files/gencode.gene.info.v22.tsv
GeneInfo <- read.table(file = "gencode.gene.info.v22.tsv", sep = "\t", head=TRUE, colClasses = c("character", "character"))

load("brca-RNAseq-Counts.rda")
seq.brca <- data

seqRes <- assay(seq.brca)
rownames(seqRes) <- rowData(seq.brca)$external_gene_name
geneLength <-GeneInfo[match(rowData(seq.brca)$original_ensembl_gene_id, GeneInfo$gene_id),c("gene_id", "gene_name", "exon_length")]

##Calculate the RPK value
RPK <- matrix(0, nrow=dim(seqRes)[1], ncol=dim(seqRes)[2])

for(row in 1:dim(seqRes)[1]){
  for(col in 1:dim(seqRes)[2]){
    RPK[row,col] <- seqRes[row,col]/as.numeric(geneLength$exon_length[row])
  }
}

##Calculate the sums of each column and divide by 1000000
scale_factor <- colSums(RPK)/1000000

##Now divide all values in each column by the scaling factor
TPM <- t(t(RPK)/scale_factor)
colnames(TPM) <- colnames(seqRes)
row.names(TPM) <- row.names(seqRes)

save(TPM, file = "BRCA_TPM.Rda")

load("BRCA_TPM.Rda")
KeepRow <- apply(TPM, 1, function(x) sum(x) != 0)
TCGA.RNAseq.dat <- TPM[KeepRow,]

# Samples to be removed based outlier based on EMT genes (all normal), metastatic samples, and outliers based on HK genes (all tumor)
RejS = c("TCGA-BH-A0B5-11A-23R-A12P-07", "TCGA-E2-A1BC-11A-32R-A12P-07", "TCGA-BH-A0H9-11A-22R-A466-07", 
         "TCGA-BH-A1EO-11A-31R-A137-07", "TCGA-A7-A0DB-11A-33R-A089-07", "TCGA-BH-A0DD-11A-23R-A12P-07", 
         "TCGA-A7-A0D9-11A-53R-A089-07", "TCGA-BH-A0B8-11A-41R-A089-07", "TCGA-E2-A15I-11A-32R-A137-07", 
         "TCGA-E2-A158-11A-22R-A12D-07", "TCGA-BH-A18V-06A-11R-A213-07", "TCGA-BH-A1ES-06A-12R-A24H-07", 
         "TCGA-E2-A15A-06A-11R-A12D-07", "TCGA-E2-A15E-06A-11R-A12D-07", "TCGA-E2-A15K-06A-11R-A12P-07", 
         "TCGA-A7-A13E-01B-06R-A277-07", "TCGA-A7-A0DC-01A-11R-A00Z-07", "TCGA-A7-A0DB-01C-02R-A277-07", 
         "TCGA-A7-A13D-01B-04R-A277-07", "TCGA-A7-A0DB-01A-11R-A277-07", "TCGA-A7-A0DC-01B-04R-A22O-07", 
         "TCGA-A7-A13D-01A-13R-A277-07", "TCGA-A7-A13E-01A-11R-A277-07")

TCGA.RNAseq.dat2 <- TCGA.RNAseq.dat[, !(colnames(TCGA.RNAseq.dat) %in% RejS)]

BCType <- substr(colnames(TCGA.RNAseq.dat2),14,15)
DotColors <- c("blue", "green", "red")
ColorDot <- DotColors[as.factor(BCType)]

# Housekeeping gene scaling
# data: https://github.com/KlinkeLab/DigitalCytometry_EMT_2020/blob/master/Files/HousekeepingGenes-PMID23810203.txt
HKGenes <- read.table("HousekeepingGenes-PMID23810203.txt", strip.white = TRUE, head=FALSE, sep = "\t", colClasses = c("character"))

RNAseq.HK <- TCGA.RNAseq.dat2[rownames(TCGA.RNAseq.dat2) %in% HKGenes$V1, ]
SCALE.HK <- apply(RNAseq.HK, 2, median)/35.33 # 35.33 average expression of HK genes in CCLE cell lines

Median.HK <- rep(0, ncol(RNAseq.HK))
plot(density(log2(RNAseq.HK[,1] + 0.01), adj = 0.1), xlab = "Expression (TPM)", ylab = "Density", ylim = c(0,10), col = ColorDot[1])
Median.HK[1] = median(RNAseq.HK[,1])
for (i in 2:ncol(RNAseq.HK))
{
  tmp <- density(log2(RNAseq.HK[,i] + 0.01), adj = 0.1)
  lines(tmp$x, tmp$y + 0.1*i, col = ColorDot[i])
  Median.HK[i] = median(RNAseq.HK[,i])
}  

plot(Median.HK, type = "p", col = ColorDot)

GEData_HK <- TCGA.RNAseq.dat2
for (i in 1:ncol(GEData_HK))
{
  GEData_HK[,i] <- TCGA.RNAseq.dat2[ , i]/SCALE.HK[i]
}  

save(GEData_HK, file = "BRCA_TCGA_TPM_HK.rda")
