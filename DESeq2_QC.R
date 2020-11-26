rm(list = ls())

library(DESeq2)
library(tidyverse)

library(vsn)
library(pheatmap)
library(RColorBrewer)
library(apeglm)
library(svglite)
library(Cairo)

library(hciR)
library(hciRdata)
library(DEGreport)

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

# create working directory
workingDir <- "201119_WorkingDirectory_DESeq2"
dirPath <- file.path(PrimaryDirectory, workingDir)
dir.create(dirPath)
setwd(dirPath)

savePath <- file.path(dirPath, "results")
saveDir <- file.path(dirPath, "results")
dir.create(savePath)

# load counts and coldata
filesPath <- file.path(dirPath, "files")
countFile <- "Counts"

cts <- read_tsv(file.path(filesPath, countFile))
is_tibble(cts)
rowname_cts <- cts
names(rowname_cts)[1] <- "rowname"
rowname_cts <- column_to_rownames(rowname_cts)

# remove D0 samples
grep("D0", colnames(rowname_cts))
rowname_cts <- rowname_cts[, -grep("D0", colnames(rowname_cts))]

# remove DA samples
# grep("DA", colnames(rowname_cts))
# rowname_cts <- rowname_cts[, -grep("DA", colnames(rowname_cts))]
colnames(rowname_cts)

# only integer values
rowname_cts <- round(rowname_cts)

coldataFile <- "coldata.txt"

colData <- read_tsv(file.path(filesPath, coldataFile))
is_tibble(colData)
names(colData)[1] <- "rowname"
col.data <- column_to_rownames(colData)
# remove D0 and DA samples

grep("D0", rownames(col.data))
col.data <- col.data[-grep("D0", rownames(col.data)), ]

# grep("DA", rownames(col.data))
# col.data <- col.data[-grep("DA", rownames(col.data)), ]
rownames(col.data)

# 2 indicates columns, 1 indicates rows
col.data <- as.data.frame(apply(col.data, 2, as.factor))


summary(col.data)
all(colnames(rowname_cts) == rownames(col.data))

# Running DESeq2

# setting data for DESeq2

dds <- DESeqDataSetFromMatrix(countData = rowname_cts, colData = col.data, 
                              design = ~Condition + Region + ID)

# multifactorial design

# group includes BM vs PB and Hyp vs Norm, BM_Norm as a reference
dds$group <- factor(paste0(dds$Region, "_", dds$Condition))
levels(dds$group) <- c("BM_Norm", "BM_Hyp", "PBL_Norm", "PBL_Hyp")
design(dds) <- formula(~ ID + group)

# pre-filtering the dataset
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 10,]
nrow(dds)

# Variance Stabilizing Transformation
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

vsd <- vst(dds, blind = TRUE)
svg(file.path(savePath, "PCAplot_allPts.svg"))
DESeq2::plotPCA(vsd, intgroup = "group")
dev.off()
png(file.path(savePath, "PCAplot_allPts.png"))
DESeq2::plotPCA(vsd, intgroup = "group")
dev.off()

dataPCA <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
(outliers <- dataPCA %>% filter(PC1 < -50))
plot(assay(vsd)[, 1:2], pch = 16, cex = 0.3)
plot(assay(vsd)[, 2:3], pch = 16, cex = 0.3)
plot(assay(vsd)[, c(1,16)], pch = 16, cex = 0.3)

# try and remove DA samples