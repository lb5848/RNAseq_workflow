rm(list = ls())

library(BiocParallel)
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
library(limma)
library(data.table)
library(ggplotify)



# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

# create working directory
workingDir <- "201119_WorkingDirectory_DESeq2_full_workflow"
dirPath <- file.path(PrimaryDirectory, workingDir)
dir.create(dirPath)
setwd(dirPath)

savePath <- file.path(dirPath, "results")
saveDir <- file.path(dirPath, "results")
dir.create(savePath)

register(MulticoreParam(4))

# load files
filesPath <- file.path(dirPath, "files")
countFile <- "Counts"

# load and prep cts file
cts <- read_tsv(file.path(filesPath, countFile))
is_tibble(cts)
rowname_cts <- cts
names(rowname_cts)[1] <- "rowname"
rowname_cts <- column_to_rownames(rowname_cts)
colnames(rowname_cts)
# remove D0 samples
grep("D0", colnames(rowname_cts))
rowname_cts <- rowname_cts[, -grep("D0", colnames(rowname_cts))]

# remove DA samples
grep("DA", colnames(rowname_cts))
rowname_cts <- rowname_cts[, -grep("DA", colnames(rowname_cts))]
colnames(rowname_cts)

# only integer values
rowname_cts <- round(rowname_cts)

# load and prep col.data file
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

################### Running DESeq2 ##########################
# LRT
ddsLRT <- DESeqDataSetFromMatrix(countData = rowname_cts, colData = col.data, 
                                 design = ~ ID + Condition)
# Wald
dds <- DESeqDataSetFromMatrix(countData = rowname_cts, colData = col.data, 
                              design = ~ ID + Condition)
# setup multifactorial design

# create "group" - ?levels "BM_Norm", "PBL_Norm", "BM_Hyp", "PBL_Hyp"

ddsLRT$group <- factor(paste0(ddsLRT$Region, "_", ddsLRT$Condition),
                       levels = c("BM_Norm", "PBL_Norm", "BM_Hyp", "PBL_Hyp"))
design(ddsLRT) <- formula(~ ID + group)

dds$group <- factor(paste0(ddsLRT$Region, "_", ddsLRT$Condition),
                    levels = c("BM_Norm", "PBL_Norm", "BM_Hyp", "PBL_Hyp"))
design(dds) <- formula(~ ID + group)

# Pre-Filtering
keep <- rowSums( counts(ddsLRT) ) >= 15
ddsLRT <- ddsLRT[ keep, ]

# varianceStabilizingTransformation
vsdLRT <- vst(ddsLRT, blind = FALSE)
vsdLRT_mat <- assay(vsdLRT)
plot_pca(vsdLRT, intgroup = "group")


DESeq2::plotPCA(vsdLRT, intgroup = "group")
ggsave(file.path(savePath, "PCA_4pts_activated.svg"), plot = last_plot())
ggsave(file.path(savePath, "PCA_4pts_activated.png"), plot = last_plot())

# run DESeq2, test = LRT, reduced = ~ ID
ddsLRT <- DESeq(ddsLRT, test = "LRT", reduced = ~ ID)

plotDispEsts(ddsLRT)
plot_dist(vsdLRT, intgroup = "group")

resLRT <- results(ddsLRT, name = resultsNames(ddsLRT)[7])
p.adj.cutoff <- 0.01
log2FC.cutoff <- 2

sig_resLRT <- resLRT %>% data.frame() %>% 
  rownames_to_column(var = "id") %>%
  as_tibble() %>%
  filter(padj < p.adj.cutoff) %>%
  filter(abs(log2FoldChange) > log2FC.cutoff)

plot_volcano(sig_resLRT)

sig_resLRT <- sig_resLRT %>% arrange(padj)
cluster_vsd <- vsdLRT_mat[ sig_resLRT$id, ]
meta <- colData(vsdLRT)
meta$Time <- rep("D3", 16)
clusters <- degPatterns(cluster_vsd, metadata = meta, time = "group", minc = 15)

# rldLRT <- rlog(ddsLRT)
sig_resLRThciR <- sig_resLRT %>% data.frame() %>% column_to_rownames( var = "id" ) %>%
  rownames_to_column(var = "id") %>%
  as_tibble()

x <- top_counts(sig_resLRThciR, vsdLRT, top = 1000, filter = TRUE, sort_fc = TRUE)

plot_genes(x, intgroup = "group", scale = "row", show_rownames = FALSE, annotation_names_col = FALSE, 
           show_colnames = TRUE)
out <- plot_genes(x, intgroup = "group", scale = "row", show_rownames = FALSE, 
                  annotation_names_col = FALSE, show_colnames = FALSE, output = "pheatmap")
plot <- plot_genes(x, intgroup = "group", scale = "row", show_rownames = FALSE, 
                   annotation_names_col = FALSE, show_colnames = FALSE, output = "pheatmap")
plot <- as.ggplot(plot, scale = 1, hjust = 0, vjust = 0)
ggsave(file.path(savePath, "heatmap1000genes.svg"), plot = plot)
ggsave(file.path(savePath, "heatmap1000genes.png"), plot = plot)

out.clust <- cbind(x, cluster = sort(cutree(out$tree_row, k = 6)))
out.clust <- as_tibble(out.clust)
out.clust.2 <- out.clust %>% filter(out.clust$cluster == 6) %>% select(-cluster)
plot.clust <- x %>% filter(id %in% out.clust.2$id)
plot_genes(plot.clust, intgroup = "group", scale = "row", show_rownames = TRUE, 
           annotation_names_col = FALSE, show_colnames = FALSE, output = "pheatmap")


# DESeq2 w/ Wald test

# Wald
dds <- DESeqDataSetFromMatrix(countData = rowname_cts, colData = col.data, 
                              design = ~ ID + Condition)
# setup multifactorial design

# create "group" - ?levels "BM_Norm", "PBL_Norm", "BM_Hyp", "PBL_Hyp"

dds$group <- factor(paste0(ddsLRT$Region, "_", ddsLRT$Condition),
                    levels = c("BM_Norm", "PBL_Norm", "BM_Hyp", "PBL_Hyp"))
design(dds) <- formula(~ ID + group)

# Pre-Filtering
keep <- rowSums( counts(dds) ) >= 15
dds <- dds[ keep, ]

# varianceStabilizingTransformation
vsd <- vst(dds, blind = TRUE)
DESeq2::plotPCA(vsd, intgroup = "group")

dds <- DESeq(dds, test = "Wald")
resultsNames(dds)

res <- results_all(dds, vs = "all", trt = "group")
plot_volcano(res[[1]])
res_sig <- res %>% lapply(function(x) {
  x <- x %>% as_tibble() %>% filter(padj < p.adj.cutoff) %>% filter(abs(log2FoldChange) > log2FC.cutoff)})
plot_volcano(res_sig[[6]])
res_sig
# save csv files - significant genes only - all contrasts
lapply(1:length(res_sig), function(i){
  res_sig[[i]] %>% as_tibble() %>% arrange(padj) %>% fwrite(file.path(savePath, paste0(names(res_sig[i]), ".csv")))
})

names(res_sig)
res_sighciR <- res_sig %>% lapply(function(x){
  x <- x %>% data.frame() %>% column_to_rownames(var = "id") %>% rownames_to_column(var = "id") %>% as_tibble()
})

names(res_sighciR)
# BM_Norm vs PBL_Norm
x <- top_counts(res_sighciR[[1]], vsd, top = 1000, filter = TRUE, sort_fc = TRUE)

plot_genes(x, intgroup = "group", scale = "row", show_rownames = FALSE, annotation_names_col = FALSE, 
           show_colnames = TRUE)

#BM_Norm vs BM_Hyp
x <- top_counts(res_sighciR[[2]], vsd, top = 150, filter = TRUE, sort_fc = TRUE)
noMS <- x %>% select(-grep("MS", colnames(x)))
plot_genes(x, intgroup = "group", scale = "row", show_rownames = FALSE, annotation_names_col = FALSE, 
           show_colnames = FALSE)
plot_genes(noMS, intgroup = "group", scale = "row", show_rownames = FALSE, annotation_names_col = FALSE, 
           show_colnames = FALSE)

# PBL_Norm vs BM_Hyp
x <- top_counts(res_sighciR[[4]], vsd, top = 2800, filter = TRUE, sort_fc = FALSE)

plot_genes(x, intgroup = "group", scale = "row", show_rownames = FALSE, annotation_names_col = FALSE, 
           show_colnames = TRUE)

# BM_Hyp vs PBL_Hyp
x <- top_counts(res_sighciR[[6]], vsd, top = 40, filter = TRUE, sort_fc = TRUE)

plot_genes(x, intgroup = "group", scale = "row", show_rownames = TRUE, annotation_names_col = FALSE, 
           show_colnames = TRUE)
