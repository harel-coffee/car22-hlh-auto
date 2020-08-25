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


if __name__ == '__main__':

    # Set variables and path ---------------------------------------------------
    save_figs = True
    filename = snakemake.input[0]

    # New data version 1 -------------------------------------------------------
    data = pyh22.load_car22(filename, access_var="all", drop=True, version=1)

    cytokine_plot = pyhty.transform_cytokines(
        data.cytokines_days_num,
        primary_outcome = data.primary_outcome)

    # Recode HLH column
    cytokine_plot = cytokine_plot.replace(
        {'HLH': {0:'carHLH-', 1:'carHLH+'}})

    cytokine_plot['cytokine_levels'] = cytokine_plot['cytokine_levels'].astype(float)
    cytokine_plot['log_cytokine_levels'] = np.log(cytokine_plot['cytokine_levels'])

    # Sort by specific cytokine order
    cytokine_plot['cytokines'] = cytokine_plot['cytokines'].astype('category')
    cytokine_plot['cytokines'].cat.reorder_categories(
        ['IFN-gamma', 'IL-1B', 'IL-6', 'IL-8','IL-10',
        'IL-12p70', 'IL-13', 'MIP-alpha', 'TNF-alpha', 'IL-4',
        'IL-2', 'GM-CSF', 'IL-15', 'IL-18'], inplace=True)
    cytokine_plot = cytokine_plot.sort_values(by=['cytokines'])

    # sort HLH category
    cytokine_plot['HLH'] = cytokine_plot['HLH'].astype('category')
    cytokine_plot['HLH'].cat.reorder_categories(['carHLH-', 'carHLH+'], inplace=True)

    # Select specific days
    cytokine_plot_days = cytokine_plot.loc[cytokine_plot['days'].isin(
        [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,99,100])]

    cytokine_plot_days = cytokine_plot_days.replace(
            {'days': {99:21, 100:28}})

    cytokine_plot_days['days'].unique()

    axes = pyhty.features_over_time(
        rows=3,
        columns=5,
        features=cytokine_plot_days,
        hue='HLH',
        feature_column="cytokines",
        y_axis='log_cytokine_levels',
        x_axis='days',
        estimator=np.mean,
        err_style='band',
        ci=68,
        units=None,
        figsize_height=10,
        figsize_width=20,
        fig_title=None,
        fig_title_pos=1.02,
        xlabel="Day after infusion",
        ylabel="log$_{10}$(Cytokine levels) (pg/mL)",
        palette=["#3498db", "#e74c3c"],
        draw_line=False)

    plt.savefig(snakemake.output[0])
