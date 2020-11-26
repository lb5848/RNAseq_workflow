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


# VST (Variance Stabilizing Transformation) transformation - for applications OTHER THAN differential testing
vsd <- vst(dds, blind = TRUE)
plot(assay(vsd)[, 2:3], pch = 16, cex = 0.3)

plot_dist(vsd, intgroup = "group")

plot_pca(vsd, intgroup = "group", pc = c(1,2))
plot_pca(vsd, intgroup = "Condition", pc = c(1,2))
plot_pca(vsd, intgroup = "Region", pc = c(1,2))
plot_pca(vsd, intgroup = "ID", pc = c(1,2))
dataPCA <- DESeq2::plotPCA(vsd, intgroup = "group", returnData = TRUE)

outliers <- dataPCA %>% filter(dataPCA$PC1 < - 50)

dds <- estimateSizeFactors(dds)

sizeFactors(dds)
# DA_BM samples appear to be outliers -> re-run DESeq after removing them
# ASK GREG
cts <- read_tsv(file.path(filesPath, countFile))
is_tibble(cts)
rowname_cts <- cts
names(rowname_cts)[1] <- "rowname"
rowname_cts <- column_to_rownames(rowname_cts)

# remove D0 samples
grep("D0", colnames(rowname_cts))
rowname_cts <- rowname_cts[, -grep("D0", colnames(rowname_cts))]

# remove DA samples
grep("DA", colnames(rowname_cts))
rowname_cts <- rowname_cts[, -grep("DA", colnames(rowname_cts))]
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

grep("DA", rownames(col.data))
col.data <- col.data[-grep("DA", rownames(col.data)), ]
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


# VST (Variance Stabilizing Transformation) transformation - for applications OTHER THAN differential testing
vsd <- vst(dds, blind = TRUE)
plot_pca(vsd, intgroup = "group", pc = c(1,2))

svg(file.path(savePath,"PCAplot_4pts.svg"))
DESeq2::plotPCA(vsd, intgroup = "group")
dev.off()

# run DESeq2 pipeline
dds <- DESeq(dds, test = "LRT", reduced = ~ ID)
res_LRT <- results(dds)
p.adj.cutoff <- 0.01
sig_res_LRT <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(padj < p.adj.cutoff)
dim(sig_res_LRT)
sigLRT_genes <- sig_res_LRT %>%
  pull(gene)

clustering_sig_genes <- sig_res_LRT %>%
  arrange(padj)

vsd_mat <- assay(vsd)
cluster_vsd <- vsd_mat[clustering_sig_genes$gene,]

clusters <- degPatterns(cluster_vsd, metadata = colData(dds), time = "group")
