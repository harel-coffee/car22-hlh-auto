#!/usr/bin/env Rscript

options(warn=1)
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("RColorBrewer"))

parser <- ArgumentParser()
parser$add_argument(
    "--eset-file", type="character", help="input eset file rds", required=TRUE
)
args <- parser$parse_args()

prior_count <- 1
k <- 2
fig_dim <- 5
fig_res <- 300
ylim <- c(-0.6, 0.6)

cat("Loading:", args$eset_file, "\n")
eset <- readRDS(args$eset_file)

dataset_name <- gsub(
    "_counts", "", tools::file_path_sans_ext(basename(args$eset_file))
)

if (!dir.exists("results")) dir.create("results")

colors <- brewer.pal(5, "Set2")
legend <- c("Batch1", "Batch2", "Batch4", "Batch5", "Batch6")
batch <- as.integer(eset$Batch)

# counts
cat("Analyzing counts\n")
results_file <- paste0(
    "results/", paste(dataset_name, "counts_rle.png", sep="_")
)
cat("Creating:", results_file, "\n")
png(file=results_file, width=fig_dim, height=fig_dim, units="in", res=fig_res)
plotRLE(exprs(eset), col=colors[batch], outline=FALSE, ylim=c(-3, 3))
title(paste(dataset_name, "Counts RLE"), cex.main=0.8)
legend("topright", legend=legend, col=colors, pch=15, cex=0.5)
invisible(dev.off())

results_file <- paste0(
    "results/", paste(dataset_name, "counts_pca.png", sep="_")
)
cat("Creating:", results_file, "\n")
png(file=results_file, width=fig_dim, height=fig_dim, units="in", res=fig_res)
plotPCA(exprs(eset), col=colors[batch])
title(paste(dataset_name, "Counts PCA"), cex.main=0.8)
legend("topright", legend=legend, col=colors, pch=15, cex=0.5)
invisible(dev.off())

results_file <- paste0(
    "results/", paste(dataset_name, "counts_mds.png", sep="_")
)
cat("Creating:", results_file, "\n")
png(
    file=results_file, width=fig_dim, height=fig_dim, units="in", res=fig_res
)
plotMDS(exprs(eset), col=colors[batch])
title(paste(dataset_name, "Counts MDS"), cex.main=0.8)
legend("topright", legend=legend, col=colors, pch=15, cex=0.5)
invisible(dev.off())

# DESeq2
cat("Analyzing DESeq2 RLE-VST\n")
dds <- DESeqDataSetFromMatrix(exprs(eset), pData(eset), ~Batch + Class)
dds <- estimateSizeFactors(dds, quiet=TRUE)
dds <- estimateDispersions(dds, fitType="mean", quiet=TRUE)
vst <- varianceStabilizingTransformation(dds, blind=FALSE)
vsd <- assay(vst)

results_file <- paste0(
    "results/", paste(dataset_name, "deseq2_vsd.rds", sep="_")
)
cat("Saving:", results_file, "\n")
saveRDS(vsd, file=results_file)

results_file <- paste0(
    "results/", paste(dataset_name, "deseq2_rle.png", sep="_")
)
cat("Creating:", results_file, "\n")
png(file=results_file, width=fig_dim, height=fig_dim, units="in", res=fig_res)
plotRLE(vsd, col=colors[batch], outline=FALSE, ylim=ylim)
title(paste(dataset_name, "DESeq2 RLE-VST RLE"), cex.main=0.8)
legend("topright", legend=legend, col=colors, pch=15, cex=0.5)
invisible(dev.off())

results_file <- paste0(
    "results/", paste(dataset_name, "deseq2_pca.png", sep="_")
)
cat("Creating:", results_file, "\n")
png(file=results_file, width=fig_dim, height=fig_dim, units="in", res=fig_res)
plotPCA(vsd, col=colors[batch])
title(paste(dataset_name, "DESeq2 RLE-VST PCA"), cex.main=0.8)
legend("topright", legend=legend, col=colors, pch=15, cex=0.5)
invisible(dev.off())

results_file <- paste0(
    "results/", paste(dataset_name, "deseq2_mds.png", sep="_")
)
cat("Creating:", results_file, "\n")
png(
    file=results_file, width=fig_dim, height=fig_dim, units="in", res=fig_res
)
plotMDS(vsd, col=colors[batch])
title(paste(dataset_name, "DESeq2 RLE-VST MDS"), cex.main=0.8)
legend("topright", legend=legend, col=colors, pch=15, cex=0.5)
invisible(dev.off())
