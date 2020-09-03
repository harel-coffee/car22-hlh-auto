# Script to reproducte parts of supplementary table 3

# Load libraries
library(tidyverse)

# Start workflow
crs <- read_csv(snakemake@input[[1]])
crs$log10_cytokine_levels <- log10(crs$cytokine_levels)

crs_stats <- crs %>%
  dplyr::group_by(days_rel_to_CRS, cytokines, HLH) %>%
  dplyr::summarise(
    mean = mean(cytokine_levels, na.rm=TRUE), 
    sd = sd(cytokine_levels, na.rm=TRUE),
    median = median(cytokine_levels, na.rm=TRUE),
    IQR_stat = IQR(cytokine_levels, na.rm = TRUE),
    log_mean = mean(log10(cytokine_levels), na.rm=TRUE),
    log_sd = sd(log10(cytokine_levels), na.rm=TRUE),
    log_median = median(log10(cytokine_levels), na.rm=TRUE),
    log_IQR_stat = IQR(log10(cytokine_levels), na.rm = TRUE),
    n = n())

#write_excel_csv(
#  crs_stats, 
#  path=paste(dir_data, "rel_to_CRS_stats.csv", sep=""))

# Significant timepoints
timepoints <- sort(unique(crs$days_rel_to_CRS))
cytokines <- unique(crs$cytokines)

# test
#timepoint <- -10
#cytokine <- "GM-CSF"

# initialize
summary_table <- data.frame(
  timepoint_rel_CRS = integer(),
  cytokine = character(),
  wilcox_pvalue = numeric(),
  log10_wilcox_pvalue = numeric())
summary_table$cytokine <- as.character(summary_table$cytokine)
i <- 1

for (timepoint in timepoints){
  #print(timepoint)
  for (cytokine in cytokines){
    #print(cytokine)
    
    df <- crs %>% dplyr::filter(days_rel_to_CRS==timepoint & cytokines==cytokine)
    wilcox_res <- wilcox.test(
      dplyr::filter(df, HLH=="carHLH-" )$cytokine_levels, 
      dplyr::filter(df, HLH=="carHLH+" )$cytokine_levels,
      alternative = c("two.sided"))
    
    log10_wilcox_res <- wilcox.test(
      dplyr::filter(df, HLH=="carHLH-" )$log10_cytokine_levels, 
      dplyr::filter(df, HLH=="carHLH+" )$log10_cytokine_levels,
      alternative = c("two.sided"))
    
    summary_table[i,]$timepoint_rel_CRS <- timepoint
    summary_table[i,]$cytokine <- cytokine
    summary_table[i,]$wilcox_pvalue <- wilcox_res$p.value
    summary_table[i,]$log10_wilcox_pvalue <- log10_wilcox_res$p.value
    
    i <- i + 1
    
    #print(df)
  }
}

# Merge tables
crs_stats_hlh_pos <- crs_stats %>% dplyr::filter(HLH == "carHLH+") %>%
  dplyr::rename(
    carHLHpos_days_rel_to_CRS = days_rel_to_CRS,
    carHLHpos_cytokines = cytokines,
    carHLHpos_mean = mean, 
    carHLHpos_sd = sd, 
    carHLHpos_median = median,
    carHLHpos_IQR_stat = IQR_stat,
    carHLHpos_log_mean = log_mean, 
    carHLHpos_log_sd = log_sd, 
    carHLHpos_log_median = log_median,
    carHLHpos_log_IQR_stat = log_IQR_stat,
    carHLHpos_n = n) %>%
  dplyr::select(-HLH)


crs_stats_hlh_neg <- crs_stats %>% dplyr::filter(HLH == "carHLH-") %>%
  dplyr::rename(
    carHLHneg_days_rel_to_CRS = days_rel_to_CRS,
    carHLHneg_cytokines = cytokines,
    carHLHneg_mean = mean, 
    carHLHneg_sd = sd, 
    carHLHneg_median = median,
    carHLHneg_IQR_stat = IQR_stat,
    carHLHneg_log_mean = log_mean, 
    carHLHneg_log_sd = log_sd, 
    carHLHneg_log_median = log_median,
    carHLHneg_log_IQR_stat = log_IQR_stat,
    carHLHneg_n = n) %>%
  dplyr::select(-HLH)

summary_table <- as_tibble(summary_table) %>%
  dplyr::mutate(padj_wilcox = p.adjust(wilcox_pvalue, method="BH")) %>%
  dplyr::mutate(padj_log10wilcox = p.adjust(log10_wilcox_pvalue, method="BH"))

master_table <- bind_cols(crs_stats_hlh_pos, crs_stats_hlh_neg, summary_table) %>%
  dplyr::select(carHLHpos_cytokines,
                carHLHpos_days_rel_to_CRS, 
                carHLHpos_mean, carHLHneg_mean,
                carHLHpos_sd, carHLHneg_sd,
                carHLHpos_median, carHLHneg_median,
                carHLHpos_IQR_stat, carHLHneg_IQR_stat,
                wilcox_pvalue, padj_wilcox,
                carHLHpos_log_mean, carHLHneg_log_mean,
                carHLHpos_log_sd, carHLHneg_log_sd,
                carHLHpos_log_median, carHLHneg_log_median,
                carHLHpos_log_IQR_stat, carHLHneg_log_IQR_stat,
                log10_wilcox_pvalue, padj_log10wilcox,
                carHLHpos_n, carHLHneg_n) %>%
  dplyr::rename(days_relative_to_CRS = carHLHpos_days_rel_to_CRS, 
                cytokines = carHLHpos_cytokines)

write_excel_csv(
  master_table, 
  path=snakemake@output[[1]])