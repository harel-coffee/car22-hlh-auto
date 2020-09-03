#!/usr/bin/env phython
import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pyhumboldt.car22 as pyh22
from pyhumboldt.car22 import Cytokines
import pyhumboldt.ploty as pyhty
import pyhumboldt.unsupervised as pyhed
import missingno as msno


if __name__ == '__main__':

    # Set variables and path ---------------------------------------------------
    save_figs = True
    filename = snakemake.input[0]

    # New data version 1 -------------------------------------------------------
    data = pyh22.load_car22(filename, access_var="all", drop=True, version=1)
    cytokines = pyh22.Cytokines(data.cytokines)
    cytokines.df = data.cytokines_days_num

    # Plot missing data
    #msno.matrix(cytokines.df)

    # Get days relative to CRS
    days_to_index = data.cytokines_days_num.unstack().unstack(level=1)
    days_to_index.index.rename(names='date', level=0, inplace=True)

    # Merge with clinical data
    days_outcome = pd.merge(
        data.secondary_outcome.reset_index(),
        days_to_index.reset_index(),
        on='patient_id')

    # Calculate days relative to CRS
    days_outcome['days_rel_to_CRS'] = days_outcome.date - days_outcome.date_CRS

    # Select specific days
    cytokines_p2 = days_outcome.loc[days_outcome['days_rel_to_CRS'].isin(range(-10, 11))].reset_index()

    # Add HLH info
    Xclin = pd.merge(
        cytokines_p2,
        data.primary_outcome.reset_index(),
        on='patient_id')

    # Related to supplementary table 3------------------------------------------
    # Calculate a separate table: with median, mean, IQR and SD with stats for
    # each timepoint to determine difference, and the number evaluable at each
    # timepoint (stratified by HLH and not-HLH)

    # Melt for plotting
    days_melt = pd.melt(
        Xclin,
        id_vars=['index', 'patient_id', 'date_CRS', 'max_grade_CRS', 'CRS',
                'HLH', 'days_HLH', 'date', 'days_rel_to_CRS'],
        var_name='cytokines',
        value_name='cytokine_levels')

    # Recode HLH column
    days_melt = days_melt.replace(
        {'HLH': {0:'carHLH-', 1:'carHLH+'}})

    # Save dataframe
    days_melt.to_csv(snakemake.output[0], index=False)
