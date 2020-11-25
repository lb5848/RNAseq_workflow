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
colnames(rowname_cts)

# only integer values
rowname_cts <- round(rowname_cts)

coldataFile <- "coldata.txt"

colData <- read_tsv(file.path(filesPath, coldataFile))
is_tibble(colData)
names(colData)[1] <- "rowname"
col.data <- column_to_rownames(colData)
grep("D0", rownames(col.data))
col.data <- col.data[-grep("D0", rownames(col.data)), ]
rownames(col.data)

# 2 indicates columns, 1 indicates rows
col.data <- as.data.frame(apply(col.data, 2, as.factor))


summary(col.data)
all(colnames(rowname_cts) == rownames(col.data))

# 1- include D0 samples
# Running DESeq2

# setting data for DESeq2
dds <- DESeqDataSetFromMatrix(countData = rowname_cts, colData = col.data, 
                              design = ~Condition + Region + ID)

# multifactorial design

dds$group <- factor(paste0(dds$Region, "_", dds$Condition))
design(dds) <- formula(~group + ID)

# pre-filtering the dataset
dds <- dds[ rowSums(counts(dds)) > 1,]

dds <- DESeq(dds)

# VST transformation - for applications OTHER THAN differential testing
vsd <- vst(dds, blind = TRUE)

plot_pca(vsd, intgroup = "group", pc = c(1,3))
