#!/usr/bin/env Rscript

options(warn=1)
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("EnhancedVolcano"))
suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("RUVSeq"))

parser <- ArgumentParser()
parser$add_argument(
    "--eset-file", type="character", help="input eset file rds", required=TRUE
)
args <- parser$parse_args()

prior_count <- 1
k <- 2
fc <- 1.0
lfc <- log2(fc)
fdr <- 0.2
fig_dim <- 8
fig_res <- 300
xlim <- c(-3, 2.5)
ylim <- c(0, 4)

cat("Loading:", args$eset_file, "\n")
eset <- readRDS(args$eset_file)

dataset_name <- gsub(
    "_counts", "", tools::file_path_sans_ext(basename(args$eset_file))
)

if (!dir.exists("results")) dir.create("results")

hk_gene_idxs <- which(
    fData(eset)$Code.Class %in% c("Control", "Housekeeping", "housekeeping")
)

# DESeq2
cat("Running DESeq2\n")
dds <- DESeqDataSetFromMatrix(exprs(eset), pData(eset), ~Batch + Class)
rowData(dds) <- fData(eset)
dds <- DESeq(dds, fitType="mean", quiet=TRUE)
results <- as.data.frame(results(
    dds, name=resultsNames(dds)[length(resultsNames(dds))], alpha=fdr,
    lfcThreshold=lfc, altHypothesis="greaterAbs", pAdjustMethod="BH"
))
results$Code.Class <- rowData(dds)$Code.Class
results$Accession <- rowData(dds)$Accession
results <- results[, c(
    "Code.Class", "Accession",
    colnames(results)[!(colnames(results) %in% c("Code.Class", "Accession"))]
)]
if (lfc > 0 && ("svalue" %in% colnames(results))) {
    results <- results[order(results$svalue), ]
} else {
    results <- results[order(results$pvalue), ]
}
results_file <- paste0(
    "results/", paste(dataset_name, "deseq2_results.tsv", sep="_")
)
cat("Writing:", results_file, "\n")
write.table(
    data.frame("Symbol"=row.names(results), results), file=results_file,
    sep="\t", quote=FALSE, row.names=FALSE
)
results_file <- paste0(
    "results/", paste(dataset_name, "deseq2_volcano.png", sep="_")
)
cat("Creating:", results_file, "\n")
png(file=results_file, width=fig_dim, height=fig_dim, units="in", res=fig_res)
EnhancedVolcano(
    results,
    lab=row.names(results),
    x="log2FoldChange",
    y=ifelse(lfc > 0 && ("svalue" %in% colnames(results)), "svalue", "padj"),
    ylab=bquote(~-Log[10]~adjusted~italic(P)),
    xlim=xlim,
    ylim=ylim,
    pCutoff=fdr,
    FCcutoff=lfc,
    pointSize=3.0,
    labSize=3.0,
    drawConnectors=FALSE,
    colConnectors="grey50",
    title=paste(dataset_name, "differential expression"),
    subtitle=paste(
        "DESeq2: RLE normalization, Wald test",
        ifelse(lfc > 0, ", lfcThreshold", "")
    ),
    caption=paste("FC cutoff = ", fc, "; FDR cutoff = ", fdr, sep="")
)
invisible(dev.off())
