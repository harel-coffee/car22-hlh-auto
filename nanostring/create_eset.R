#!/usr/bin/env Rscript

options(warn=1)
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("NanoStringNorm"))

parser <- ArgumentParser()
parser$add_argument(
    "--rcc-dir", type="character", help="rcc markup file dir"
)
parser$add_argument(
    "--counts-file", type="character", help="counts file tsv"
)
parser$add_argument(
    "--meta-file", type="character", help="metadata file tsv", required=TRUE
)
parser$add_argument(
    "--out-file", type="character", help="output eset file rds", required=TRUE
)
parser$add_argument(
    "--class-type", type="character", help="class label type",
    choices=c("hlh", "tcs", "response"), required=TRUE
)
args <- parser$parse_args()

if (is.null(args$rcc_dir)) {
    if (is.null(args$counts_file)) stop("Specify --rcc-dir or --counts-file")
    counts <- read.delim(
        args$counts_file, row.names=NULL, stringsAsFactors=FALSE
    )
    counts[, 1:3] <- counts[, c("code_class", "gene_name", "accession")]
    colnames(counts)[1:3] <- c("Code.Class", "Name", "Accession")
} else {
    rcc <- read.markup.RCC(args$rcc_dir)
    cat("\n")
    counts <- rcc$x
    sample_cols <- vapply(
        strsplit(colnames(counts)[4:length(colnames(counts))], "_", fixed=TRUE),
        "[", "", 1
    )
    colnames(counts) <- c(
        "Code.Class", "Name", "Accession",
        gsub("^X", "P", sample_cols, perl=TRUE)
    )
}

row.names(counts) <- counts[, 2]
fdata <- counts[, c(1, 3), drop=FALSE]
counts[, 1:3] <- NULL
counts <- counts[, order(colnames(counts)), drop=FALSE]
if (!identical(row.names(fdata), row.names(counts)))
    stop("counts and fdata row names not identical")

pdata <- read.delim(args$meta_file, row.names=1)
pdata <- pdata[order(row.names(pdata)), , drop=FALSE]
if (!identical(row.names(pdata), colnames(counts)))
    stop("pdata row names and counts col names not identical")

if (args$class_type == "hlh") {
    pdata$Class <- pdata$HLH
} else if (args$class_type == "tcs") {
    pdata$Class <- pdata$TCS
} else if (args$class_type == "response") {
    pdata$Class <- ifelse(
        pdata$outcome == "OR", 1, ifelse(pdata$outcome == "NR", 0, NA)
    )
}

pdata$Class <- as.factor(pdata$Class)
pdata$Batch <- as.factor(pdata$batch)
pdata$Group <- NULL

eset <- ExpressionSet(
    assayData=as.matrix(counts),
    phenoData=AnnotatedDataFrame(pdata),
    featureData=AnnotatedDataFrame(fdata)
)

eset <- eset[, !is.na(pdata$Class)]

# drop outliers
eset <- eset[, !(colnames(eset) %in% c("P3333", "P4661", "P9691"))]
if (args$class_type == "hlh") {
    # drop cases with no CRS
    eset <- eset[, !(colnames(eset) %in% c("P992", "P1056", "P1639"))]
}
pData(eset) <- gdata::drop.levels(pData(eset))
saveRDS(eset, file=args$out_file)
