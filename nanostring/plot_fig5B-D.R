#!/usr/bin/env Rscript

options(warn=1)
debug <- TRUE

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("NanoStringNorm"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(ComplexHeatmap))

parser <- ArgumentParser()
parser$add_argument(
  "--vsd-file", type="character", help="VSD file.", required=TRUE
)
parser$add_argument(
  "--deseq2-res", type="character", help="DESeq2 output results"
)
parser$add_argument(
  "--meta_file", type="character", help="RCC meta data file.", required=TRUE
)
parser$add_argument(
  "--eset-file", type="character", help="input eset file rds", required=TRUE
)
parser$add_argument(
  "--dir-results", type="character", help="Results directory", required=FALSE
)
args <- parser$parse_args()

if(debug){print(args)}

# Read in vsd.
vsd <- readRDS(args$vsd_file)
if(debug){print(str(vsd))}

# Plot MDS (Figure 5B)
mat_mds <- plotMDS(vsd)

mds_df <- tibble(
  patient_name = names(mat_mds$x),
  MDS1 = mat_mds$x, 
  MDS2 = mat_mds$y) %>%
  left_join(read_tsv(args$meta_file), 
            by="patient_name") %>%
  dplyr::mutate(HLH = as.factor(HLH)) %>%
  dplyr::mutate(TCS = as.factor(TCS)) %>%
  dplyr::mutate(batch = as.factor(batch))

if(debug){print(mds_df)}

# Plot MDS
g <- ggplot(data=mds_df, aes(x=MDS1, y=MDS2))
g <- g + geom_point(aes(x=MDS1, y=MDS2, fill=TCS, color=TCS,  pch=HLH), 
                    size=2, stroke=1)
g <- g + scale_shape_manual(values=c(1, 3))
g <- g + scale_fill_manual(values=c("#3498db", "#e74c3c", "black"))
g <- g + scale_color_manual(values=c("#3498db", "#e74c3c", "black"))
#g <- g + geom_text_repel(data=mds_df, aes(x=MDS1, y=MDS2, label=patient_name), size=3)
g <- g + theme_bw()
g <- g + theme(
  aspect.ratio=1,
  axis.text = element_text(size = 12))

pdf(width=5, height=5, 
    file=paste(args$dir_results, "Figure5B_deseq2_MDS_mode=symbol_colorbyTCS.pdf", sep=""),
    useDingbats = FALSE)
plot(g)
dev.off()

# Vulcano plot (Figure 5C)
# Select significant genes.

# settings
# mode (TCS/HLH)
mode <- "TCS"
# pvalue cutoff
cutoffpvaladj <- 0.2
# log fold change cutoff
log2fc <- 0

deseq2 <- read_tsv(args$deseq2_res) %>%
  dplyr::mutate(significant = ifelse(padj < cutoffpvaladj & (log2FoldChange > log2fc | log2FoldChange < log2fc), "yes", "no")) %>%
  dplyr::mutate(significant=replace(significant, (padj < cutoffpvaladj & log2FoldChange > log2fc), "up")) %>%
  dplyr::mutate(significant=replace(significant, (padj < cutoffpvaladj & log2FoldChange < log2fc), "down"))

# Plot vulcano
g <- ggplot(data=deseq2, aes(x=log2FoldChange, y=-log10(padj)))
g <- g + geom_vline(xintercept = 0, linetype="dashed", color="grey", size=0.7)
#g <- g + geom_vline(xintercept = 0.26, linetype="dashed", color="grey", size=0.7)
#g <- g + geom_vline(xintercept = -0.26, linetype="dashed", color="grey", size=0.7)
g <- g + geom_hline(yintercept = -log10(cutoffpvaladj), linetype="dashed", color="grey", size=0.7)
g <- g + geom_point(aes(x=log2FoldChange, y=-log10(padj), fill=significant, color=significant),
                    pch=21, size=2, alpha=0.5)
g <- g + scale_fill_manual(values=c("blue", "grey", "#e41a1c"))
g <- g + scale_color_manual(values=c( "blue", "grey", "#e41a1c"))
g <- g + geom_text_repel(
  data=dplyr::filter(deseq2, significant!="no"), 
  aes(x=log2FoldChange, y=-log10(padj), label=Symbol),
  arrow = arrow(length = unit(0.02, "npc"), type = "closed", ends = "first"),
  segment.color = "grey50",
  size=3)
g <- g + xlim(c(-2.5, 2.5))
g <- g + ylim(c(0, 2.5))
g <- g + xlab("Log2(Fold change)")
g <- g + ylab("-Log10(Adjusted P value)")
g <- g + theme_bw()
g <- g + theme(
  aspect.ratio=1,
  axis.text = element_text(size = 12))

pdf(width=5, height=5, 
    file=paste(args$dir_results, "Figure5C_deseq2_replot_sig_genes=0p2_log2=0_mode=", mode, ".pdf", sep=""), 
    useDingbats = FALSE)
plot(g)
dev.off()

# Heatmap (Figure 5D)
# Filter
sig <- dplyr::filter(deseq2, significant!="no") %>% dplyr::arrange(log2FoldChange)
vsd_sig <- vsd[sig$Symbol, ]
vsd_sig_scaled <- t(apply(vsd_sig, 1, scale))
colnames(vsd_sig_scaled) <- colnames(vsd_sig)

# Read in meta file
eset <- readRDS(args$eset_file)
if(debug){print(eset)}

# Add annotation
TCS <- eset$TCS
names(TCS) <- paste("P", eset$randomized_id, sep="")
HLH <- eset$HLH
names(HLH) <- paste("P", eset$randomized_id, sep="")
CRS <- eset$CRS
names(CRS) <- paste("P", eset$randomized_id, sep="")
outcome <- eset$outcome
names(outcome) <- paste("P", eset$randomized_id, sep="")
batch <- eset$batch2
names(batch) <- paste("P", eset$randomized_id, sep="")

# for publication
mat_anno2 <- HeatmapAnnotation(
  TCS=as.factor(TCS),
  HLH=as.factor(HLH),
  col = list(TCS = c("0"="#3498db", "1"="#e74c3c"), 
             HLH = c("0"="#3498db", "1"="#e74c3c")),
  gp = gpar(col = "white"),
  na_col = "black")

# Plot heatmap
pdf(file=paste(args$dir_results, "Figure5D_deseq2_heatmap_clustering=kmeans_signficant_genes=0p2_log2_TCS_full_anno.pdf", sep=""), 
    height=6, width=13, useDingbats = FALSE)
Heatmap(
  vsd_sig_scaled,
  name="z-score", 
  cluster_rows=TRUE,
  top_annotation = mat_anno2,
  rect_gp = gpar(col = "white", lwd = 1),
  width = unit(13, "cm"), 
  height = unit(6, "cm"))
dev.off()
