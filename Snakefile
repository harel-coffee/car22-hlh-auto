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
FIG1 = '%s1.pdf' % (FIGS_PREFIX)
FIG2 = '%s2.pdf' % (FIGS_PREFIX)
#FIG3 = '%s3.pdf' % (FIGS_PREFIX)

rule all:
    input:
        FIG1,
        FIG2

rule relative_to_CRS:
    input:
        CLINICAL_DATA
    output:
        FIG1,
        FIG2
    script:
        "relative_to_CRS/main_figures.py"
