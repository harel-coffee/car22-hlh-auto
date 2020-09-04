#!/usr/bin/env Rscript

options(warn=1)
debug <- TRUE

# Load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(clusterProfiler)) 
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(ReactomePA))
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(pathview))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument(
  "--deseq2-res", type="character", help="DESeq2 output results"
)
parser$add_argument(
  "--var", type="character", help="Variable of interest.", choices=c("HLH", "TCS"), required=TRUE
)
parser$add_argument(
  "--sortby", type="character", help="Sort gene list by pvalue or wald statistic.", choices=c("PVAL", "STAT"), required=TRUE
)
parser$add_argument(
  "--pval-threshold", type="double", help="P-value threshold (default: 0.2)", required=FALSE
)
parser$add_argument(
  "--fc-threshold", type="double", help="Fold change threshold (default: 0.2)", required=FALSE
)
parser$add_argument(
  "--dir-results", type="character", help="Results directory", required=FALSE
)
args <- parser$parse_args()

if(debug){print(args)}

if (is.null(args$pval_threshold)) {
  args$pval_threshold <- 0.2
}

if (is.null(args$fc_threshold)) {
  args$fc_threshold <- 0.2
}

# Read in deseq2 output file
deseq2 <- read_tsv(args$deseq2_res) %>%
  dplyr::mutate(significant = ifelse(
    padj < args$pval_threshold & (log2FoldChange > args$fc_threshold | log2FoldChange < -(args$fc_threshold)), "yes", "no"))

if(debug){print(deseq2)}

# Load Mart
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

# Extract ids and gene description
ids <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol', 
                 'entrezgene_id', 'chromosome_name', 'description'), 
  filter = 'hgnc_symbol',
  values = deseq2$Symbol, 
  mart = ensembl)

# Filter (only keep ids with numeric chromosome location)
ids_filtered <- as_tibble(ids) %>%
  dplyr::filter(!grepl("CHR", chromosome_name)) %>%
  dplyr::distinct(hgnc_symbol, .keep_all=TRUE)

# Try mapping IDs with org.Hs.eg.db
ids2 <- mapIds(org.Hs.eg.db, deseq2$Symbol, 'ENTREZID', 'SYMBOL')
idsHsdb <- as.data.frame(ids2)
idsHsdb$symbol <- names(ids2)
idsHsdb$ids2 <- as.character(idsHsdb$ids2)

# Manually fill in missing annotation
manual_anno <- tibble(
  symbol = c("GARS", "TRBC1/2", "GLUD1/2", "KIR3DL1/2", "CCL4/L1",
             "CCL3/L1", "TRAV36", "XCL1/2", "TRAV38-2", "TRAV29",
             "FCGR3A/B", "TPSAB1/B2", "TRAV23", "CD45RA", "CD45R0",
             "MINOS1", "KLRC1/2", "TRAV14", "TRGV3/5"),
  symbol_fix = c("GARS1", "TRBC2", "GLUD1", "KIR3DL1", "CCL4", 
                 "CCL3", "TRAV36DV7", "XCL1", "TRAV38.2", "TRAV29DV5",
                 "FCGR3A", "TPSAB1", "TRAV23DV6", "PTPRC", "PTPRC",
                 "MICOS10", "KLRC1", "TRAV14DV4", "TRGV3"),
  entrez_ids_manual = c("2617", "28638", "2746", "3811", "6351", 
                        "6348", "28646", "6375", "100652412", "28653",
                        "2214", "7177", "28660", "5788", "5788", 
                        "440574", "3821", "28669", "6976"))

# Merge with deseq2
deseq2match <- deseq2 %>%
  dplyr::left_join(ids_filtered, by = c("Symbol"="hgnc_symbol")) %>%
  dplyr::left_join(idsHsdb, by = c("Symbol"="symbol")) %>%
  dplyr::left_join(manual_anno, by = c("Symbol"="symbol")) %>%
  dplyr::filter(Code.Class=="Endogenous") %>%
  dplyr::mutate(entrez_ids = ids2) %>%
  dplyr::mutate(entrez_ids = ifelse(is.na(ids2), entrez_ids_manual, ids2)) %>%
  dplyr::rename(gene_symbol=Symbol, alt_gene_symbol=symbol_fix,
                code_class=Code.Class, accession=Accession) %>%
  dplyr::select(gene_symbol, alt_gene_symbol, entrez_ids, accession, 
                description, code_class, baseMean, log2FoldChange, lfcSE,
                stat, pvalue, padj, significant)

write_excel_csv(
  as.data.frame(deseq2match), 
  path=paste(args$dir_results, "DESEQ2_NANOSTRING_GENE_MAPPING_",
             args$var, ".csv", sep=""))

# Only select differentially expressed genes (downregulated)
deseq2_sig_down <- dplyr::filter(deseq2match, significant=="yes" & log2FoldChange < 0)
genes <- deseq2_sig_down$log2FoldChange
names(genes) <- deseq2_sig_down$entrez_ids

# Original gene list & rank gene list by log2fc/logfcSE
universe <- deseq2match$stat
names(universe) <- deseq2match$entrez_ids
universe_sorted <- sort(universe, decreasing = TRUE)

#################
# GO ENRICHMENT #
#################

# set max set size no larger than 10% of panel size
maxGSSize <- 770 * 0.1

# Overrepresentation Analysis (ORA)
go_results1 <- enrichGO(
  gene = deseq2_sig_down$entrez_ids,
  universe = names(universe_sorted),
  OrgDb = "org.Hs.eg.db",
  ont='ALL',
  pvalueCutoff = 1,
  minGSSize    = 15,
  maxGSSize    = maxGSSize,
  qvalueCutoff = 1, 
  readable=TRUE)

write_excel_csv(
  as.data.frame(go_results1), 
  path=paste(args$dir_results, 
             "GO_ORA_UNIVERSE=NANOSTRING_DOWNREG_TRESH=LOG2_0_ONTOLOGY=ALL_MIN=10_MAX=", 
             maxGSSize, "_SORT=", args$sortby, "_COMP=", args$var, ".csv", sep=""))

# Gene Set Enrichment Analysis (GSEA)
gsea_go <- gseGO(
  geneList     = universe_sorted,
  OrgDb        = org.Hs.eg.db,
  ont          = "ALL",
  nPerm        = 1000,
  minGSSize    = 15,
  maxGSSize    = maxGSSize,
  pvalueCutoff = 1,
  verbose      = FALSE)

gsea_go_y <- setReadable(
  gsea_go, 
  OrgDb = org.Hs.eg.db, 
  keyType="ENTREZID")

# Plot enrichment plots
# for (i in 1:30){
#  p1 <- enrichplot::gseaplot2(
#    gsea_go_y, 
#    geneSetID = i,
#    title = gsea_go_y$Description[i],
#    pvalue_table = TRUE)
#  
#  pdf(width=9, height=5, 
#      file=paste(args$dir_results, 
#                 "enrichment_plots/", "GO_", 
#                 gsea_go_y$Description[i],"_SORT=", 
#                 args$sortby, "_COMP=", args$var, ".pdf", sep=""))
#  plot(p1)
#  dev.off()
#}

# Write results to file
write_excel_csv(
  as.data.frame(gsea_go_y), 
  path=paste(args$dir_results, 
             "GO_GSEA_UNIVERSE=NANOSTRING_ONTOLOGY=ALL_MIN=15_MAX=", maxGSSize, "_SORT=", 
             args$sortby, "_COMP=", args$var,".csv", sep=""))

# Enrichment bar plot (ORA)
go_enrichment <- as_tibble(gsea_go_y) %>%
  dplyr::mutate(direction=ifelse(NES<0, "down", "up")) %>%
  dplyr::filter(qvalues < 0.1) %>%
  dplyr::arrange(NES) 

go_top10 <- head(go_enrichment, 10)
if(debug){print(go_top10)}

# Plot barplot
#pdf(file=paste(args$dir_results, 
#               "barplot_go_gsea.pdf", sep=""), 
#    height=2.5, width=8, useDingbats = FALSE)
#ggbarplot(go_top10,
#               x = "Description", y = "NES",
#               color = "black",            # Set bar border colors to white
#               fill = "direction",
#               palette = "jco",            # jco journal color palett. see ?ggpar
#               sort.val = "desc",          # Sort the value in descending order
#               sort.by.groups = FALSE,     # Don't sort inside each group
#               x.text.angle = 90,          # Rotate vertically x axis texts
#               ylab = "NES (Enrichment Score)",
#               width=0.85, 
#               legend.title = "Enrichment",
#               rotate = TRUE,
#               ggtheme = theme_minimal())
#dev.off()

###################
# KEGG ENRICHMENT #
###################

# overrepresentation analysis ORA
kk <- enrichKEGG(
  gene = deseq2_sig_down$entrez_ids,
  universe = names(universe_sorted), 
  organism  = 'hsa',
  minGSSize = 15,
  maxGSSize = maxGSSize,
  pvalueCutoff = 1,
  qvalueCutoff = 1)

kk_y <- setReadable(
  kk, 
  OrgDb = org.Hs.eg.db, 
  keyType="ENTREZID")

# plot barplot
pdf(file=paste(args$dir_results, 
               "barplot_kegg_ORA.pdf", sep=""), 
    height=1.5, width=7, useDingbats = FALSE)
barplot(kk_y, showCategory=5)
dev.off()

write_excel_csv(
  as.data.frame(kk_y), 
  path=paste(args$dir_results, 
             "KEGG_ORA_UNIVERSE=NANOSTRING_DOWNREG_TRESH=LOG2_0_MIN=10_MAX=", maxGSSize, "_SORT=",
             args$sortby, "_COMP=", args$var, ".csv",  sep=""))

# GSEA
kk2 <- gseKEGG(
  geneList     = universe_sorted,
  organism     = 'hsa',
  nPerm        = 1000,
  minGSSize = 15,
  maxGSSize = maxGSSize,
  pvalueCutoff = 1,
  verbose      = FALSE,
  use_internal_data = FALSE)

kk2_y <- setReadable(
  kk2, 
  OrgDb = org.Hs.eg.db, 
  keyType="ENTREZID")

for (i in 1:10){
  p1 <- enrichplot::gseaplot2(
    kk2_y, 
    geneSetID = i,
    title = kk2_y$Description[i],
    pvalue_table = TRUE)
  
  pdf(width=9, height=5, 
      file=paste(args$dir_results, "KEGG_", 
                 gsub(" ", "_", kk2_y$Description[i]),"_SORT=", args$sortby, 
                 "_COMP=", args$var, ".pdf", sep=""))
  plot(p1)
  dev.off()
}

write_excel_csv(
  as.data.frame(kk2_y), 
  path=paste(args$dir_results, 
             "KEGG_GSEA_MIN=10_MAX=",maxGSSize, "_SORT=",
             args$sortby,"_COMP=",args$var, ".csv", sep=""))

# Enrichment plot (GSEA)
kegg_enrichment <- as_tibble(kk2_y) %>%
  dplyr::mutate(direction=ifelse(NES<0, "down", "up")) %>%
  dplyr::filter(qvalues < 0.2) %>%
  dplyr::arrange(NES) 

kegg_top10 <- head(kegg_enrichment, 10)


pdf(file=paste(args$dir_results, 
               "barplot_kegg_gsea.pdf", sep=""), 
    height=1.8, width=6, useDingbats = FALSE)
ggbarplot(kegg_top10,  
          x = "Description", y = "NES",
          color = "white",            # Set bar border colors to white
          fill = "direction",
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "desc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "NES (Enrichment Score)",
          width=1, 
          legend.title = "Enrichment",
          rotate = TRUE,
          ggtheme = theme_minimal())
dev.off()

#######################
# REACTOME ENRICHMENT #
#######################
reac <- enrichPathway(
  gene = deseq2_sig_down$entrez_ids,
  universe = names(universe_sorted), 
  organism = "human",
  minGSSize = 15,
  maxGSSize = maxGSSize,
  pvalueCutoff = 1,
  qvalueCutoff = 1)

reac_y <- setReadable(
  reac, 
  OrgDb = org.Hs.eg.db, 
  keyType="ENTREZID")

write_excel_csv(
  as.data.frame(reac_y), 
  path=paste(args$dir_results, "REACTOME_ORA_DOWNREG_THRES_LOG2_0_MIN=10_MAX=", maxGSSize, "_SORT=", 
             args$sortby, "_COMP=",args$var,".csv", sep=""))

gsea_reac <- gsePathway(
  gene = universe_sorted,
  organism = "human",
  nPerm = 1000,
  minGSSize = 10,
  maxGSSize = maxGSSize,
  pvalueCutoff = 1)

gsea_reac_y <- setReadable(
  gsea_reac, 
  OrgDb = org.Hs.eg.db, 
  keyType="ENTREZID")

#for (i in 1:10){
#  p1 <- enrichplot::gseaplot2(
#    gsea_reac_y, 
#    geneSetID = i,
#    title = gsea_reac_y$Description[i],
#    pvalue_table = TRUE)
#  
#  pdf(width=9, height=5, 
#      file=paste(dir_results, "enrichment_plots/", "REACTOME_", 
#                 gsea_reac_y$Description[i], "_SORT=", diff_sort, 
#                 "_COMP=",comparison,".pdf", sep=""))
#  plot(p1)
#  dev.off()
#}

write_excel_csv(
  as.data.frame(gsea_reac_y), 
  path=paste(args$dir_results, "REACTOME_GSEA_MIN=10_MAX=", maxGSSize, "_SORT=", args$sortby, ".csv", sep=""))
