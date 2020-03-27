# **An unsupervised feature extraction and selection strategy for identifying epithelial-mesenchymal transition state metrics in breast cancer and melanoma**

This repository supplies the code developed in the study of D.J. Klinke and A. Torang **_"An unsupervised feature extraction and selection strategy for identifying epithelial-mesenchymal transition state metrics in breast cancer and melanoma"_**. The corresponding pre-print can be found on bioRxiv ( doi: https://doi.org/10.1101/865139 ). It can be used to reproduce the results of the study and investigate the methodology to be used for other datasets.

## **Requirements**

* R version 3.6.2.
* R libraries: glmnet, caret, tidyverse, colorspace, MASS, SummarizedExperiment, TCGAbiolinks

## **Data**

All necessary data is provided in the "Files" folder of the repository or linkes to public data have been provided in the scripts.

For more information about data used in this study please refer to [https://portals.broadinstitute.org/ccle] for Files: CCLE_RNAseq_081117.rpkm.gct accessed 12/22/2017 and CCLE_Expression_Entrez_2012-10-18.res accessed 6/15/2018; [https://tcpaportal.org/mclp] for File: MCLP-v1.1-Level4.txt accessed 6/15/2018; [https://www.ncbi.nlm.nih.gov] entries GSE75688, GSE72056, and GSE98394; and [https://www.ebi.ac.uk] entry E-MTAB-6831.

## **Quick start**

To reproduce the results, download the relevant script and data and load the corresponding data into the R workspace. Running the scripts will generate all relevant figures and data tables for the given portion of the study.

# General notes

The code provided in this repository reproduces the main results of the study of D.J. Klinke and A. Torang **_"An unsupervised feature extraction and selection strategy for identifying epithelial-mesenchymal transition state metrics in breast cancer and melanoma"_** but it is not meant as a self-contained module for use in further analysis.

## Citation

Klinke, D.J. & Torang, A. **_"An unsupervised feature extraction and selection strategy for identifying epithelial-mesenchymal transition state metrics in breast cancer and melanoma."_** bioRxiv (2019). [https://doi.org/10.1101/865139]
