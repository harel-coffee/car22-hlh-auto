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
    #figs_path = "/Users/schischlikf2/datasets/CAR-T/figs_main"

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

    # Plot A -------------------------------------------------------------------

    # Select specific cytokines
    days_melt_sel = days_melt.loc[days_melt['cytokines'].isin(['IFN-gamma',
        'IL-6', 'TNF-alpha', 'IL-1B'])]
    days_melt_sel['cytokine_levels'] = days_melt_sel['cytokine_levels'].astype(float)
    days_melt_sel['log_cytokine_levels'] = np.log(days_melt_sel['cytokine_levels'])

    # Plot patients over time
    axes = pyhty.features_over_time(
        rows=1,
        columns=4,
        features=days_melt_sel,
        hue='HLH',
        feature_column="cytokines",
        y_axis='log_cytokine_levels',
        x_axis='days_rel_to_CRS',
        estimator=np.mean,
        err_style='band',
        ci=68,
        units=None,
        figsize_height=3,
        figsize_width=13,
        fig_title=None,
        fig_title_pos=1.02,
        xlabel="Day relative to CRS onset",
        ylabel="log$_{10}$(Cytokine levels) (pg/mL)",
        palette=["#3498db", "#e74c3c"],
        draw_line=True)

    plt.savefig(snakemake.output[0])

    # Plot B -------------------------------------------------------------------
    # Main figure revised
    # Main figure (cytokine levels relative to CRS onset)
    # Re-order the cytokines
    # Top row: (R to L): IFNy, IL-1Beta, IL-6, IL-8, IL-10
    # Middle row: (R to L): IL12p70, IL13, MIP1a, TNFa, IL4
    # Bottom row: (R to L). IL2, GMCSF, IL-15, IL-18
    # Leave with Error band
    # Remove legends from each individual file. # post modification
    # y-axis revise to: log10 (cytokine level) (pg/mL) [PLEASE NOTE I ADDED the “10”]
    # x-axis revise to: Day relative to CRS onset [Please note that I changed “Days” to “Day”]

    # Plot all cytokines
    days_melt['cytokine_levels'] = days_melt['cytokine_levels'].astype(float)
    days_melt['log_cytokine_levels'] = np.log(days_melt['cytokine_levels'])
    #days_melt['log_cytokine_levels'] = np.log(days_melt['cytokine_levels'] + 1)

    # Sort by specific cytokine order
    days_melt['cytokines'] = days_melt['cytokines'].astype('category')
    days_melt['cytokines'].cat.reorder_categories(
        ['IFN-gamma', 'IL-1B', 'IL-6', 'IL-8','IL-10',
        'IL-12p70', 'IL-13', 'MIP-alpha', 'TNF-alpha', 'IL-4',
        'IL-2', 'GM-CSF', 'IL-15', 'IL-18'], inplace=True)
    days_melt = days_melt.sort_values(by=['cytokines'])

    # sort HLH category
    days_melt['HLH'] = days_melt['HLH'].astype('category')
    days_melt['HLH'].cat.reorder_categories(['carHLH-', 'carHLH+'], inplace=True)

    # Plot patients over time
    axes = pyhty.features_over_time(
        rows=3,
        columns=5,
        features=days_melt,
        hue='HLH',
        feature_column="cytokines",
        y_axis='log_cytokine_levels',
        x_axis='days_rel_to_CRS',
        estimator=np.mean,
        err_style='band',
        ci=68,
        units=None,
        figsize_height=10,
        figsize_width=20,
        fig_title=None,
        fig_title_pos=1.02,
        xlabel="Day relative to CRS onset",
        ylabel="log$_{10}$(Cytokine levels) (pg/mL)",
        palette=["#3498db", "#e74c3c"],
        draw_line=True)

    plt.savefig(snakemake.output[1])

    # TABLE A ------------------------------------------------------------------
    # Calculate a separate table: with median, mean, IQR and SD with stats for
    # each timepoint to determine difference, and the number evaluable at each
    # timepoint (stratified by HLH and not-HLH)
    # Are you able to do this table in BOTH: Actual value and log10 of cytokines

    # Melt for plotting
    #days_melt = pd.melt(
    #    Xclin,
    #    id_vars=['index', 'patient_id', 'date_CRS', 'max_grade_CRS', 'CRS',
    #            'HLH', 'days_HLH', 'date', 'days_rel_to_CRS'],
    #    var_name='cytokines',
    #    value_name='cytokine_levels')

    # Recode HLH column
    #days_melt = days_melt.replace(
    #    {'HLH': {0:'carHLH-', 1:'carHLH+'}})

    # save dataframe
    #days_melt.to_csv(
    #os.path.join(figs_path,"rel_to_CRS_export.csv"),
    #index=False)
