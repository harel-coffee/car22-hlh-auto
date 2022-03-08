# car22-hlh
Code to reproduce analysis for 

### Characterization of HLH-like manifestations as a CRS variant in patients receiving CD22 CAR T cells
Daniel A. Lichtenstein, Fiorella Schischlik, Lipei Shao, Seth M. Steinberg, Bonnie Yates, Hao-Wei Wang, Yanyu Wang, Jon Inglefield, Alina Dulau-Florea, Francesco Ceppi, Leandro C. Hermida, Kate Stringaris, Kim Dunham, Philip Homan, Parthav Jailwala, Justin Mirazee, Welles Robinson, Karen M. Chisholm, Constance Yuan, Maryalice Stetler-Stevenson, Amanda K. Ombrello, Jianjian Jin, Terry J. Fry, Naomi Taylor, Steven L. Highfill, Ping Jin, Rebecca A. Gardner, Haneen Shalabi, Eytan Ruppin, David F. Stroncek, Nirali N. Shah. Blood (2021) 138 (24): 2469–2484. https://doi.org/10.1182/blood.2021011898

### Set up environment
1. Install and set up [Miniconda3](https://docs.conda.io/en/latest/miniconda.html)
```
# Create environment
$ conda env create -f envs/car22-hlh.yaml
$ conda activate car22-hlh

# Install pyhumboldt as follows
# Using pip with the following command (go to pyhumboldt directory):
$ cd pyhumboldt
$ pip install -e. --user

# (optional) uninstall pyhumboldt
$ pip uninstall
```

###  Run analysis
```
# The following command will reproduce:
# Figure 2 (Figure2_cytokines_over_time_rel_CRS.pdf)
# Figure 4C (Figure4C_modelA.pdf)
# Figure 4D (Figure4D_modelB.pdf)
# Supplementary Figure 2 (Suppl_Figure2_cytokines_over_time.pdf)
# Supplementary Table 3 (relative_to_CRS_stats.csv)
# Summary of model A (predictive-modelA-summary.txt)
# Summary of model B (predictive-modelB-summary.txt)
$ cd ..
$ snakemake --cores 1

# For nanostring analysis, go to nanostring directory and follow README instructions:
$ cd nanostring

# For changes related to revisons to to subdirectory /revisions
# File structure
# revisions/figures:
main_figures.py # code for reproducing figure 3 of the final manuscripts

# revisions/predictive_models:
logreg_relative_to_CRS.py # code for reproducing the predictive models in figure 4 of the final manuscript

# revisions/stats: # code used for reproducing stats/supplementary tables of the final manuscript
day_after_infusion_heatmap.rmd
day_after_infusion_stats.rmd
relative_to_CRS_stats.rmd
```

