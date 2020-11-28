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
library(stats)

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

savePath <- file.path(dirPath, "results_wald")
saveDir <- file.path(dirPath, "results_wald")
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
                              design = ~ID + Region + Condition + Region:Condition)

# multifactorial design

# pre-filtering the dataset
nrow(dds)
dds <- dds[rowSums(counts(dds)) > 10, ]
nrow(dds)
view(counts(dds))

# VST (Variance Stabilizing Transformation) transformation - for applications OTHER THAN differential testing
vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = "Condition")
plotPCA(vsd, intgroup = "Region")
plotPCA(vsd, intgroup = c("Region","Condition"))

svg(file.path(savePath, "PCAplot_4Pts.svg"))
DESeq2::plotPCA(vsd, intgroup = c("Region","Condition"))
dev.off()
png(file.path(savePath, "PCAplot_4Pts.png"))
DESeq2::plotPCA(vsd, intgroup = c("Region","Condition"))
dev.off()

# hierarchical clustering
vsd_mat <- assay(vsd)
vsd_cor <- cor(vsd_mat)
pheatmap(vsd_cor)
plot_dist(vsd, intgroup = c("Region","Condition"))

dds <- DESeq(dds, test = "Wald")
p.adj.cutoff <- 0.05
resBM_PB <- results(dds, contrast = c("Region", "PBL", "BM"))
res_hciR <- resBM_PB %>% data.frame() %>% 
  rownames_to_column(var = "id") %>%
  filter(padj < p.adj.cutoff) %>%
  as_tibble()


x <- top_counts(res_hciR, vsd, top = 1500, sort_fc = TRUE, filter = TRUE)
plot_genes(x, intgroup = "Region", scale = "row", annotation_names_col = FALSE, show_rownames = FALSE)

resHyp_Norm <- results(dds, contrast = c("Condition", "Norm", "Hyp"))
res_hciR <- resHyp_Norm %>% data.frame() %>% 
  rownames_to_column(var = "id") %>%
  filter(padj < p.adj.cutoff) %>%
  as_tibble()


x <- top_counts(res_hciR, vsd, top = 100, sort_fc = TRUE, filter = TRUE)
plot_genes(x, intgroup = "Condition", scale = "row", annotation_names_col = FALSE, show_rownames = FALSE)

res2 <- results(dds, contrast = list( c("Condition_Norm_vs_Hyp", "RegionPBL.ConditionNorm") ))
res_hciR <- res2 %>% data.frame() %>% 
  rownames_to_column(var = "id") %>%
  filter(padj < p.adj.cutoff) %>%
  as_tibble()


x <- top_counts(res_hciR, vsd, top = 100, sort_fc = TRUE, filter = TRUE)
plot_genes(x, intgroup = "Region", scale = "row", annotation_names_col = FALSE, show_rownames = FALSE)


