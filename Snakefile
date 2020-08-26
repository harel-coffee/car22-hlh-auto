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
MODEL_SUMMARY = join(RESULTS_DIR, 'predictive-models-summary.tsv')

# Figures
FIGS_PREFIX = join(FIGURES_DIR, 'fig')
FIG1 = '%s_main_1.pdf' % (FIGS_PREFIX)
FIG2 = '%s_main_2.pdf' % (FIGS_PREFIX)
FIG3 = '%s_main_modelA_3.pdf' % (FIGS_PREFIX)
FIG4 = '%s_main_modelB_4.pdf' % (FIGS_PREFIX)
SUPPL_FIG1 = '%s_suppl_3.pdf' % (FIGS_PREFIX)

rule all:
    input:
        FIG1,
        FIG2,
        FIG3,
        FIG4,
        SUPPL_FIG1

rule relative_to_CRS_main:
    input:
        CLINICAL_DATA
    output:
        FIG1,
        FIG2
    script:
        "relative_to_CRS/main_figures.py"

rule relative_to_CRS_suppl:
    input:
        CLINICAL_DATA
    output:
        SUPPL_FIG1
    script:
        "relative_to_CRS/supplementary_figures.py"

rule modelA:
    input:
        CLINICAL_DATA
    output:
        FIG3
    script:
        "predictive_models/model_A.py"

rule modelB:
    input:
        CLINICAL_DATA
    output:
        FIG4
    script:
        "predictive_models/model_B.py"
