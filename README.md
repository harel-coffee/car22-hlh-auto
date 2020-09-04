# car22-hlh
Code to reproduce analysis for 

### Biomarkers Associated with Hemophagocytic Lymphohistiocytosis-Like Toxicities Provide Novel Insights into CD22 CAR T-Cell Inflammatory Responses

Daniel A. Lichtenstein, Fiorella Schischlik, Lipei Shao, Seth M. Steinberg, Bonnie Yates, Hao-Wei Wang, Francesco Ceppi, Leandro C. Hermida, Yanyu Wang, Jon Inglefield, Shakuntala Rampertaap, Welles Robinson, Karen M. Chisholm, Constance Yuan, Maryalice Stetler-Stevenson, Amanda K. Ombrello, Sergio Rosenzweig, Jianjian Jin, Terry J. Fry, Steven L. Highfill, Ping Jin, Haneen Shalabi, Rebecca Gardner, Eytan Ruppin, David F. Stroncek, Nirali N. Shah (Under review, 2020)

## Set up environment
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
# For reproducing figure 2 (relative to CRS) and figure 4 (predictive models) run:
cd ..
snakemake --cores 1
# For nanostring analysis, go to nanostring directory and follow README instructions:
cd nanostring
```

