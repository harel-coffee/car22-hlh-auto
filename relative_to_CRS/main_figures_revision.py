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
from matplotlib import gridspec


if __name__ == '__main__':

    # Set variables and path ---------------------------------------------------
    save_figs = True
    filename = '~/src/2020-CAR22-toxicities/data/REAL Full Cytokine (EVEN and odd post CAR)_De-identified CD22 Data for Bioinformatics_5-12-20_v1.xlsx'
    #figs_path = "/Users/schischlikf2/datasets/CAR-T/figs/"
    figs_path = '/Users/schischlikf2/src/2020-CAR22-toxicities/figures/figs_revision/'

    # New data version 1 -------------------------------------------------------
    data = pyh22.load_car22(filename, access_var="all", drop=True, version=1)
    cytokines = pyh22.Cytokines(data.cytokines)
    cytokines.df = data.cytokines_days_num

    # Add FDR pvalues to plot---------------------------------------------------
    pvalues = pd.read_csv(os.path.join(figs_path, "rel_to_CRS_stats_carHLH_with_pvalues.csv"))

    # Plot missing data
    msno.matrix(cytokines.df)

    # Get days relative to CRS
    # crs_transform does not work on new dataset
    days_to_index = data.cytokines_days_num.unstack().unstack(level=1)
    days_to_index.index.rename(names='date', level=0, inplace=True)

    # Merge with clinical data
    days_outcome = pd.merge(
        data.secondary_outcome.reset_index(),
        days_to_index.reset_index(),
        on='patient_id')

    # Calc days relative to CRS
    days_outcome['days_rel_to_CRS'] = days_outcome.date - days_outcome.date_CRS

    # select specific days
    cytokines_p2 = days_outcome.loc[days_outcome['days_rel_to_CRS'].isin(range(-10, 11))].reset_index()
    #cytokines_p2 = days_outcome.loc[days_outcome['days_rel_to_CRS'].isin([2])].reset_index()
    #cyto_p2_drop = cytokines_p2.drop(columns=['IL-18', 'max_grade_CRS', 'CRS']).dropna()

    # add HLH info
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

    # TABLE A ------------------------------------------------------------------
    # Calculate a separate table: with median, mean, IQR and SD with stats for
    # each timepoint to determine difference, and the number evaluable at each
    # timepoint (stratified by HLH and not-HLH)
    # Are you able to do this table in BOTH: Actual value and log10 of cytokines

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

    # save dataframe
    days_melt.to_csv(os.path.join(figs_path,"rel_to_CRS_export.csv"), index=False)

    # Final version for revision -----------------------------------------------
    # This version has 3 cytokines in a row.

    days_melt['cytokine_levels'] = days_melt['cytokine_levels'].astype(float)
    days_melt['log10_cytokine_levels'] = np.log10(days_melt['cytokine_levels'])

    # Change type of category
    days_melt['cytokines'] = days_melt['cytokines'].astype('category')

    # sort HLH category
    days_melt['HLH'] = days_melt['HLH'].astype('category')
    days_melt['HLH'].cat.reorder_categories(['carHLH-', 'carHLH+'], inplace=True)
    days_melt = days_melt.sort_values(by=['HLH'])

    #features = ['IFN-gamma', 'IL-1B', 'IL-6', 'IL-8','IL-10',
    #    'IL-12p70', 'IL-13', 'MIP-alpha', 'TNF-alpha', 'IL-4',
    #    'IL-2', 'GM-CSF', 'IL-15', 'IL-18']
    features = ['IFN-gamma', 'IL-12p70', 'IL-2',
                'IL-1B', 'IL-13', 'GM-CSF',
                'IL-6', 'MIP-alpha', 'IL-15',
                'IL-8', 'TNF-alpha', 'IL-18',
                'IL-10', 'IL-4']
    palette = ["#3498db", "#e74c3c"]

    # Select only signficant datapoints
    fdr = pvalues[pvalues['padj_FDR'] < 0.05]

    fig, axes = plt.subplots(
        nrows=5,
        ncols=3,
        figsize=(11, 17), # width, height
        sharey=False)

    plt.subplots_adjust(
        wspace=0.3,
        hspace=0.5)

    # set font
    sns.set_style(
        {'font.family': 'sans-serif',
        'font.sans-serif': ['Arial']})

    for ax, feature in zip(axes.flatten(), features):

        # Get max value for cytokine
        # TODO: Better: get max of 68 confidence interval
        cyto_fdr = fdr[fdr['cytokine']==feature]
        cyto_pvalues = pvalues[pvalues['cytokine']==feature]
        if np.logical_not(cyto_fdr.empty):
            cyto_max = max(cyto_pvalues['carHLH_pos_log10_mean'])

        sns.lineplot(
            data=days_melt[days_melt['cytokines']==feature],
            x='days_rel_to_CRS',
            y='log10_cytokine_levels',
            hue='HLH',
            marker='o',
            estimator=np.mean,
            err_style='band',
            units=None,
            ci=68,
            palette=palette,
            ax=ax,
            sort=True)

        # Plot significant points
        if np.logical_not(cyto_fdr.empty):
            for i in cyto_fdr['days_rel_to_CRS']:
                ax.annotate('*', xy=(i, cyto_max), annotation_clip=False,
                    fontsize='large')

        ax.legend([],[], frameon=False)
        ax.set(title=feature)
        ax.tick_params(axis='both', top=False, right=False, reset =True)
        # set y axis label
        ax.set_ylabel("log$_{10}$(Cytokine levels) (pg/mL)")
        ax.set_xlabel("Day relative to CRS onset")
        ax.axvline(0, color="grey", linestyle="--")

    plt.savefig(os.path.join(figs_path,"log10_day_rel_CRS_pvalues_main.pdf"))
