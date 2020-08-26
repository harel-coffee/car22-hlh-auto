# car22-hlh
Code to reproduce analysis for Lichtenstein et al., 2020

## Set up

```
# Create environment
$ conda env create -f envs/car22-hlh.yaml
$ conda activate car22-hlh

# Install pyhumboldt as follows
# Using pip with the following command (go to pyhumboldt directory):
$ cd pyhumboldt
$ pip install -e. --user

# uninstall pyhumboldt
$ pip uninstall
```

###  Run analysis
```
# For reproducing figure 2 (relative to CRS)  and figure 4 (predictive models) run:
cd ..
snakemake --cores 1
# For nanostring analysis, go to nanostring directory and follow README instructions:
cd nanostring
```

