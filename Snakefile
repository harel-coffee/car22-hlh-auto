from os.path import join

####################################
# SETTINGS, FILES, AND DIRECTORIES #
####################################
# Directories
DATA_DIR = 'data'
RESULTS_DIR = 'results'
FIGURES_DIR = join(RESULTS_DIR, 'figs')

# Data files
CLINICAL_DATA = join(DATA_DIR, 'clinical_data_05-12-20_v1.xlsx')

# Output files
MODEL_A_SUMMARY = join(RESULTS_DIR, 'predictive-modelA-summary.txt')
MODEL_B_SUMMARY = join(RESULTS_DIR, 'predictive-modelB-summary.txt')
REL_TO_CRS = join(RESULTS_DIR, 'relative_to_CRS.csv')
REL_TO_CRS_STATS = join(RESULTS_DIR, 'relative_to_CRS_stats.csv')

# Figures
FIGS_PREFIX = join(FIGURES_DIR, 'fig')
FIG2 = '%s_main_Figure2.pdf' % (FIGS_PREFIX)
FIG2_4CYTO = '%s_main_Figure2_4cytokines.pdf' % (FIGS_PREFIX)
FIG4C = '%s_main_modelA_Fig4C.pdf' % (FIGS_PREFIX)
FIG4D = '%s_main_modelB_Fig4D.pdf' % (FIGS_PREFIX)
SUPPL_FIG2 = '%s_suppl_Fig2.pdf' % (FIGS_PREFIX)

rule all:
    input:
        FIG2,
        FIG2_4CYTO,
        SUPPL_FIG2,
        REL_TO_CRS,
        FIG4C,
        FIG4D,
        REL_TO_CRS_STATS

rule relative_to_CRS_main:
    input:
        CLINICAL_DATA
    output:
        FIG2,
        FIG2_4CYTO
    script:
        "relative_to_CRS/main_figures.py"

rule relative_to_CRS_suppl_fig:
    input:
        CLINICAL_DATA
    output:
        SUPPL_FIG2
    script:
        "relative_to_CRS/supplementary_figures.py"

rule relative_to_CRS_suppl_table:
    input:
        CLINICAL_DATA
    output:
        REL_TO_CRS
    script:
        "relative_to_CRS/supplementary_tables.py"

rule relative_to_CRS_suppl_table3:
    input:
        REL_TO_CRS
    output:
        REL_TO_CRS_STATS
    script:
        "relative_to_CRS/supplementary_tables.R"

rule modelA:
    input:
        CLINICAL_DATA
    output:
        MODEL_A_SUMMARY,
        FIG4C
    script:
        "predictive_models/model_A.py"

rule modelB:
    input:
        CLINICAL_DATA
    output:
        MODEL_B_SUMMARY,
        FIG4D
    script:
        "predictive_models/model_B.py"
