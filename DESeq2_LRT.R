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
                              design = ~ID + Condition + Region)

# multifactorial design

# group includes BM vs PB and Hyp vs Norm, BM_Norm as a reference
dds$group <- factor(paste0(dds$Region, "_", dds$Condition))
levels(dds$group) <- c("BM_Norm", "BM_Hyp", "PBL_Norm", "PBL_Hyp")
design(dds) <- formula(~ ID + group)

# pre-filtering the dataset
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 10,]
nrow(dds)
view(counts(dds))

# VST (Variance Stabilizing Transformation) transformation - for applications OTHER THAN differential testing
vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = "group")

svg(file.path(savePath, "PCAplot_4Pts.svg"))
DESeq2::plotPCA(vsd, intgroup = "group")
dev.off()
png(file.path(savePath, "PCAplot_4Pts.png"))
DESeq2::plotPCA(vsd, intgroup = "group")
dev.off()

# hierarchical clustering
vsd_mat <- assay(vsd)
vsd_cor <- cor(vsd_mat)
pheatmap(vsd_cor)
plot_dist(vsd, intgroup = "group")

dds <- DESeq(dds, test = "LRT", reduced = ~ID)
res <- results(dds)
summary(res)

boxplot(log10(assays(dds)[["cooks"]]), range = 0, las = 2)

# W <- res$stat
# maxCooks <- apply(assays(dds)[["cooks"]],1,max)
# idx <- !is.na(W)
# plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic", 
#      ylab="maximum Cook's distance per gene",
#      ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
# m <- ncol(dds)
# p <- 3
# abline(h=qf(.99, p, m - p))

p.adj.cutoff <- 0.05
log2FC.cutoff <- 1

# filter only significant genes
rld <- rlog(dds)

res_hciR <- res %>% data.frame() %>% 
  rownames_to_column(var = "id") %>%
  filter(padj < p.adj.cutoff) %>%
  as_tibble()


x <- top_counts(res_hciR, rld, top = 1500, sort_fc = FALSE, filter = TRUE)
plot_genes(x, intgroup = "group", scale = "row", annotation_names_col = FALSE, show_rownames = FALSE)
ggsave(file.path(savePath, "heatmap1500.svg"))


x <- top_counts(res_hciR, rld, top = 500, sort_fc = TRUE, filter = TRUE)
plot_genes(x, intgroup = "group", scale = "row", annotation_names_col = FALSE, show_rownames = FALSE)
ggsave(file.path(savePath, "heatmap500.svg"))

out <- plot_genes(x, intgroup = "group", scale = "row", annotation_names_col = FALSE, show_rownames = FALSE)

out.clust <- cbind(x, cluster = sort(cutree(out$tree_row, k = 6)))
out.clust <- as_tibble(out.clust)
summary(out.clust)

out.rld <- rld[out.clust$id, ]
out.rld <- assay(out.rld) %>% data.frame() %>% rownames_to_column(var = "id") %>% as_tibble()
# try to plot different clusters!

sig_res <- out.res %>% 
  data.frame() %>% 
  rownames_to_column(var = "gene_name") %>%
  as_tibble() %>%
  filter(padj < p.adj.cutoff)

clustering_sig_genes <- sig_res %>%
  arrange(padj)


cluster_vsd <- vsd_mat[clustering_sig_genes$gene_name,]
clusters <- degPatterns(cluster_vsd, metadata = colData(vsd), time = "group", col = NULL, consensusCluster = TRUE)
class(clusters)
clusters$df

