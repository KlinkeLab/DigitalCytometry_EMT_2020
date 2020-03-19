library(tidyverse)
library(caret)
library(glmnet)

setwd("~/Documents/Publications/EMTSignature/R")
rm(list = ls())

#Load List of Genes associatd with Epithelial and Mesenchymal state metrics based on CCLE data
load(file = "./data/BRCA_CCLE_Esig_Msig.rda")

#######################
FIBROgs.file.name <- "./data/Fibroblasts.csv"
Fibro_Genes <- read.table(FIBROgs.file.name, head=TRUE, sep = ",", stringsAsFactors = FALSE, na.strings = "null")

# Mesenchymal signature
# switch LEPRE1:P3H1
Msig[Msig == "LEPRE1"] <- "P3H1"
M_Genes <- Msig[!(Msig %in% Fibro_Genes$Gene_symbol)]

#Epithelial signature
E_Genes <- Esig[!(Esig %in% Fibro_Genes$Gene_symbol)]

#Load variable GEData_HK
load(file = "./data/BRCA_TCGA_TPM_HK.rda")

# Let's do E_Genes first
KeepRow <- match(E_Genes, rownames(GEData_HK))
if (any(is.na(KeepRow))) stop('Couldn\'t find matching gene in list')

ExpMat <- t(GEData_HK[KeepRow, ])
ExpMata <- log2(ExpMat + 0.03)
ExpMatb <- as.data.frame(scale(ExpMata, center = apply(ExpMata, 2, median) ))
ExpMat2 <- data.frame(ExpMatb, Sample_Type = as.factor(substr(rownames(ExpMatb),14,15)))

## Inspect the data
#head(ExpMat2, 4)

# Split the data into training and test set
set.seed(123)

# number of iterations
nit <- 500
nres <- matrix(rep(0, nit*(length(E_Genes) + 5)), nrow = nit)
colnames(nres) <- c("Lambda_Min", "Pred_0", "Pred_Min", "Pred_1", "Intercept", E_Genes)
for (i in 1:nit)
{
  # %>% pipe result to next function 
  training.samples <- ExpMat2[,"Sample_Type"] %>% createDataPartition(p = 0.8, list = FALSE)
  train.data  <- ExpMat2[training.samples, ]
  test.data <- ExpMat2[-training.samples, ]

  # Dummy code categorical predictor variables
  x <- model.matrix(Sample_Type~., train.data)[,-1]
  # Convert the outcome (class) to a numerical variable
  y <- ifelse(train.data$Sample_Type == "01", 1, 0)

  # Find the best lambda using cross-validation, use ridge regression as predictor
  # variables exhibit high colinearity (alpha = 0)

  # Set up a range of lambdas to test
  lambda_seq <- 10^seq(-1, -16, by = -.1)

  # Ridge regression with cross validation
  cv.ridge <- cv.glmnet(x, y, alpha = 0, family = "binomial", lambda = lambda_seq)

  t1 <- cv.ridge$lambda.min

  t5 <- coef(cv.ridge, cv.ridge$lambda.min)@x

  # Final model with lambda = 0
  ridge.model <- glmnet(x, y, alpha = 0, family = "binomial",
                      lambda = 0)
  # Make prediction on test data
  x.test <- model.matrix(Sample_Type~., test.data)[,-1]
  probabilities <- ridge.model %>% predict(newx = x.test)
  predicted.classes <- ifelse(probabilities > 0.5, "01", "11")
  # Model accuracy
  observed.classes <- test.data$Sample_Type
  t2 <- mean(predicted.classes == observed.classes)

  # Final model with lambda.min
  ridge.model <- glmnet(x, y, alpha = 0, family = "binomial", lambda = cv.ridge$lambda.min)
  # Make prediction on test data
  x.test <- model.matrix(Sample_Type~., test.data)[,-1]
  probabilities <- ridge.model %>% predict(newx = x.test)
  predicted.classes <- ifelse(probabilities > 0.5, "01", "11")
  # Model accuracy
  observed.classes <- test.data$Sample_Type
  t3 <- mean(predicted.classes == observed.classes)

  # Final model with lambda = 1
  ridge.model <- glmnet(x, y, alpha = 0, family = "binomial", lambda = 1)
  # Make prediction on test data
  x.test <- model.matrix(Sample_Type~., test.data)[,-1]
  probabilities <- ridge.model %>% predict(newx = x.test)
  predicted.classes <- ifelse(probabilities > 0.5, "01", "11")
  # Model accuracy
  observed.classes <- test.data$Sample_Type
  t4 <- mean(predicted.classes == observed.classes)

 nres[i,] = c(t1, t2, t3, t4, t5)
}

write.csv(nres, file = "BRCA-RidgeRegression-EGenes-Mar20.csv")

RR_E <- read.csv(file = "BRCA-RidgeRegression-EGenes-Mar20.csv")
E_KeepGene <- rep(TRUE, length(E_Genes))

pdf("BRCA-RidgeCoef-EGenes-Mar20.pdf", width = 7, height = 7)
opar <- par(mfrow = c(3,3))
di <- density(RR_E[,6], adj = 0.5, from = -20, to = 20)
PNeg <- sum(di$y[di$x <0])/sum(di$y)
plot(di$x, di$y, type = "l", xlim = c(-20, 20), ylim = c(0,1), main = colnames(RR_E)[6])
text(x = 0, y = 0.8, srt = 0, labels = signif(PNeg,5), pos = 2)
for (i in 7:ncol(RR_E))
{
  di <- density(RR_E[,i], adj = 0.5, from = -20, to = 20)
  plot(di$x, di$y, type = "l", xlim = c(-10, 10), ylim = c(0,1), main = colnames(RR_E)[i])
  PNeg <- sum(di$y[di$x <0])/sum(di$y)
  E_KeepGene[i-6] <- PNeg > 0.95
  text(x = 0, y = 0.8, srt = 0, labels = signif(PNeg,5), pos = 2)
}  
dev.off() #

# Now do M_Genes
KeepRow <- match(M_Genes, rownames(GEData_HK))
if (any(is.na(KeepRow))) stop('Couldn\'t find matching gene in list')
ExpMat <- t(GEData_HK[KeepRow, ])
ExpMata <- log2(ExpMat + 0.03)
ExpMatb <- as.data.frame(scale(ExpMata, center = apply(ExpMata, 2, median) ))
ExpMat2 <- data.frame(ExpMatb, Sample_Type = as.factor(substr(rownames(ExpMatb),14,15)))

# Inspect the data
#head(ExpMat2, 4)

# Split the data into training and test set
set.seed(123)

# number of iterations
nit <- 500
nres <- matrix(rep(0, nit*(length(M_Genes) + 5)), nrow = nit)
colnames(nres) <- c("Lambda_Min", "Pred_0", "Pred_Min", "Pred_1", "Intercept", M_Genes)
for (i in 1:nit)
{
  # %>% pipe result to next function 
  training.samples <- ExpMat2[,"Sample_Type"] %>% createDataPartition(p = 0.8, list = FALSE)
  train.data  <- ExpMat2[training.samples, ]
  test.data <- ExpMat2[-training.samples, ]
  
  # Dummy code categorical predictor variables
  x <- model.matrix(Sample_Type~., train.data)[,-1]
  # Convert the outcome (class) to a numerical variable
  y <- ifelse(train.data$Sample_Type == "01", 1, 0)
  
  # Find the best lambda using cross-validation, use ridge regression as predictor
  # variables exhibit high colinearity (alpha = 0)
  
  # Set up a range of lambdas to test
  lambda_seq <- 10^seq(-1, -16, by = -.1)
  
  # Ridge regression with cross validation
  cv.ridge <- cv.glmnet(x, y, alpha = 0, family = "binomial", lambda = lambda_seq)
  
  t1 <- cv.ridge$lambda.min
  
  t5 <- coef(cv.ridge, cv.ridge$lambda.min)@x
  
  # Final model with lambda = 0
  ridge.model <- glmnet(x, y, alpha = 0, family = "binomial",
                        lambda = 0)
  # Make prediction on test data
  x.test <- model.matrix(Sample_Type~., test.data)[,-1]
  probabilities <- ridge.model %>% predict(newx = x.test)
  predicted.classes <- ifelse(probabilities > 0.5, "01", "11")
  # Model accuracy
  observed.classes <- test.data$Sample_Type
  t2 <- mean(predicted.classes == observed.classes)
  
  # Final model with lambda.min
  ridge.model <- glmnet(x, y, alpha = 0, family = "binomial", lambda = cv.ridge$lambda.min)
  # Make prediction on test data
  x.test <- model.matrix(Sample_Type~., test.data)[,-1]
  probabilities <- ridge.model %>% predict(newx = x.test)
  predicted.classes <- ifelse(probabilities > 0.5, "01", "11")
  # Model accuracy
  observed.classes <- test.data$Sample_Type
  t3 <- mean(predicted.classes == observed.classes)
  
  # Final model with lambda = 1
  ridge.model <- glmnet(x, y, alpha = 0, family = "binomial", lambda = 1)
  # Make prediction on test data
  x.test <- model.matrix(Sample_Type~., test.data)[,-1]
  probabilities <- ridge.model %>% predict(newx = x.test)
  predicted.classes <- ifelse(probabilities > 0.5, "01", "11")
  # Model accuracy
  observed.classes <- test.data$Sample_Type
  t4 <- mean(predicted.classes == observed.classes)
  
  nres[i,] = c(t1, t2, t3, t4, t5)
}

write.csv(nres, file = "BRCA-RidgeRegression-MGenes-Mar20.csv")

RR_M <- read.csv(file = "BRCA-RidgeRegression-MGenes-Mar20.csv")
M_KeepGene <- rep(TRUE, length(M_Genes))

pdf("BRCA-RidgeCoef-MGenes-Mar20.pdf", width = 7, height = 7)
opar <- par(mfrow = c(3,3))
di <- density(RR_M[,6], adj = 0.5, from = -20, to = 20)
PNeg <- sum(di$y[di$x <0])/sum(di$y)
plot(di$x, di$y, type = "l", xlim = c(-20, 20), ylim = c(0,1), main = colnames(RR_M)[6])
text(x = 0, y = 0.8, srt = 0, labels = signif(PNeg,5), pos = 2)
for (i in 7:ncol(RR_M))
{
  di <- density(RR_M[,i], adj = 0.5, from = -20, to = 20)
  plot(di$x, di$y, type = "l", xlim = c(-10, 10), ylim = c(0,1), main = colnames(RR_M)[i])
  PNeg <- sum(di$y[di$x <0])/sum(di$y)
  M_KeepGene[i-6] <- PNeg < 0.05
  text(x = 0, y = 0.8, srt = 0, labels = signif(PNeg,5), pos = 2)
}  
dev.off() #

E_TGenes <- E_Genes[E_KeepGene]
M_TGenes <- M_Genes[M_KeepGene]

save(E_TGenes, M_TGenes, file = "./data/BRCA_TCGA_Esig_Msig.rda")
