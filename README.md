# car22-hlh
Code to reproduce analysis for Lichtenstein et al., 2020

## Set up

```
# Create environment
$ conda env create -f envs/environment.yaml
$ conda activate car22-hlh

# Install pyhumboldt as follows
# Using pip with the following command (go to pyhumboldt directory):
$ cd pyhumboldt
$ pip install -e. --user

# uninstall pyhumboldt
$ pip uninstall

# Run analysis
snakemake --cores 1

```

