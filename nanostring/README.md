# NanoString Analysis

## Installation

Install and set up [Miniconda3](https://docs.conda.io/en/latest/miniconda.html)

Install the conda environment:

```bash
conda env create -f envs/car22-tox-nano.yml
```

Activate the environment:

```bash
conda activate car22-tox-nano
```

Unpack NanoString raw data files:

```bash
unzip data/De-Identified_CD22_Frozen_Whole_Cohort.zip -d data
```

## Analysis

Create ExpressionSets:

```bash
./create_eset.R --class-type hlh \
--rcc-dir data/De-Identified_CD22_Frozen_Whole_Cohort \
--meta-file data/car22_tox_nano_meta.tsv
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
