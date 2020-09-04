## NanoString Analysis

Unpack NanoString raw data files:

```bash
unzip data/De-Identified_CD22_Frozen_Whole_Cohort.zip -d data
```

Create ExpressionSets:

```bash
./create_eset.R --class-type hlh \
--rcc-dir data/De-Identified_CD22_Frozen_Whole_Cohort \
--meta-file data/car22_tox_nano_meta.tsv \
--out-file data/car22_tox_hlh_nano_counts.rds

./create_eset.R --class-type tcs \
--rcc-dir data/De-Identified_CD22_Frozen_Whole_Cohort \
--meta-file data/car22_tox_nano_meta.tsv \
--out-file data/car22_tox_tcs_nano_counts.rds
```

Analyze batch effects and differential expression and produce plots:

```bash
./analyze_batch_effects.R --eset-file data/car22_tox_hlh_nano_counts.rds
./analyze_diff_expr.R --eset-file data/car22_tox_hlh_nano_counts.rds
```

```bash
./analyze_batch_effects.R --eset-file data/car22_tox_tcs_nano_counts.rds
./analyze_diff_expr.R --eset-file data/car22_tox_tcs_nano_counts.rds
```

Reproduce Figure 5B, 5C and 5D:
```bash
./plot_fig5B-D.R \
--vsd-file data/car22_tox_tcs_no_outliers_nano_deseq2_ruvg_vsd.rds \
--dir-results ..results/ \
--deseq2-res data/car22_tox_tcs_no_outliers_nano_deseq2_results.tsv \
--meta_file data/car22_tox_nano_meta.tsv \
--eset-file data/car22_tox_tcs_no_outliers_nano_counts.rds
```

Perform gene set enrichment analysis (GO, KEGG, REACTOME) (Figure 5E):
```bash
# if directory enrichment_analysis does not exist, create it
$ mkdir ../results/enrichment_analysis

# for TCS
./enrichment_analysis.R \
--deseq2-res data/car22_tox_tcs_no_outliers_nano_deseq2_results.tsv \
--var TCS \
--dir-results ../results/enrichment_analysis/ \
--sortby STAT

# for HLH
./enrichment_analysis.R \
--deseq2-res data/car22_tox_hlh_no_outliers_nano_deseq2_results.tsv \
--var HLH \
--dir-results ../results/enrichment_analysis/ \
--sortby STAT
```
