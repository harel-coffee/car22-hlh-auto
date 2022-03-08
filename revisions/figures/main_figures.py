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
    filename = 'revisions/data/clinical_data_05-12-20_v1.xlsx'
    figs_path = 'figs_revision/' # choose figure path

    # Load new CRS grading (ASBMT)
    file_asbmt_crs = 'revisions/data/ID Key_15C0029 CD22 CAR CRS Data_4-30-21_export.csv'
    asbmt_crs = pd.read_csv(file_asbmt_crs)
    asbmt_crs['Max Grade CRS_ASBMT'].value_counts()

    # New data version 1 -------------------------------------------------------
    data = pyh22.load_car22(filename, access_var="all", drop=True, version=1)
    cytokines = pyh22.Cytokines(data.cytokines)
    cytokines.df = data.cytokines_days_num

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

    # Plot A -------------------------------------------------------------------

    # Select specific cytokines
    days_melt_sel = days_melt.loc[days_melt['cytokines'].isin(['IFN-gamma',
        'IL-6', 'TNF-alpha', 'IL-1B'])]
    days_melt_sel['cytokine_levels'] = days_melt_sel['cytokine_levels'].astype(float)
    days_melt_sel['log10_cytokine_levels'] = np.log10(days_melt_sel['cytokine_levels'])

    # Plot patients over time
    axes = pyhty.features_over_time(
        rows=1,
        columns=4,
        features=days_melt_sel,
        hue='HLH',
        feature_column="cytokines",
        y_axis='log10_cytokine_levels',
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

    plt.savefig(os.path.join(figs_path,"4_log10cytokines_rel_to_CRS_day10to10.pdf"))

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
    days_melt['log10_cytokine_levels'] = np.log10(days_melt['cytokine_levels'])
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
    days_melt = days_melt.sort_values(by=['HLH'])

    # Plot patients over time
    axes = pyhty.features_over_time(
        rows=3,
        columns=5,
        features=days_melt,
        hue='HLH',
        feature_column="cytokines",
        y_axis='log10_cytokine_levels',
        x_axis='days_rel_to_CRS',
        estimator=np.mean,
        err_style='band',
        sharey=False,
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

    plt.savefig(os.path.join(figs_path,"log10_cytokines_rel_to_CRS_day10to10_ordered.pdf"))

    # For revision -------------------------------------------------------------
    # Create same table as HLH non HLH, but instead CRS

    # Recode CRS
    days_melt['CRS grade'] = days_melt['max_grade_CRS']
    days_melt = days_melt.replace(
        {'CRS grade': {1:'CRS 1-2', 2:'CRS 1-2', 3:'CRS 3-4', 4:'CRS 3-4'}})

    axes = pyhty.features_over_time(
        rows=3,
        columns=5,
        features=days_melt,
        hue='CRS grade',
        feature_column="cytokines",
        y_axis='log10_cytokine_levels',
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
        palette=["orange", "green"],
        draw_line=True)

    plt.savefig(os.path.join(figs_path,"log10_cytokines_rel_to_CRS_day10to10_CRS_grade_ordered.pdf"))

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
    #days_melt.to_csv(
    #os.path.join(figs_path,"rel_to_CRS_export.csv"),
    #index=False)

    # Loop over strata ---------------------------------------------------------
    # Plot all cytokines
    days_melt['cytokine_levels'] = days_melt['cytokine_levels'].astype(float)
    days_melt['log10_cytokine_levels'] = np.log10(days_melt['cytokine_levels'])

    # sort HLH category
    days_melt['HLH'] = days_melt['HLH'].astype('category')
    days_melt['HLH'].cat.reorder_categories(['carHLH+', 'carHLH-'], inplace=True)
    days_melt = days_melt.sort_values(by=['HLH'])

    # Recode CRS
    days_melt['CRS grade'] = days_melt['max_grade_CRS']
    days_melt = days_melt.replace(
        {'CRS grade': {1:'CRS 1-2', 2:'CRS 1-2', 3:'CRS 3-4', 4:'CRS 3-4'}})

    # sort CRS grade category
    days_melt['CRS grade'] = days_melt['CRS grade'].astype('category')
    days_melt['CRS grade'].cat.reorder_categories(['CRS 1-2', 'CRS 3-4'], inplace=True)
    days_melt = days_melt.sort_values(by=['CRS grade'])

    stratas = ['HLH', 'CRS grade']
    color_palettes = [["#3498db", "#e74c3c"], ["orange", "green"]]
    features = ['IFN-gamma', 'IL-1B', 'IL-6', 'IL-8','IL-10',
    'IL-12p70', 'IL-13', 'MIP-alpha', 'TNF-alpha', 'IL-4',
    'IL-2', 'GM-CSF', 'IL-15', 'IL-18']

    s = stratas * len(features)
    c = color_palettes * len(features)
    f = list(np.repeat(features, len(stratas)))

    fig = plt.figure(figsize=(7, 45))
    gs = fig.add_gridspec(14, 2, hspace=0.3)
    axes = gs.subplots(sharex='col', sharey='row')

    # set font
    sns.set_style(
        {'font.family': 'sans-serif',
        'font.sans-serif': ['Arial']})

    for ax, strata, palette, feature in zip(axes.flatten(), s, c, f):

        sns.lineplot(
            data=days_melt[days_melt['cytokines']==feature],
            x='days_rel_to_CRS',
            y='log10_cytokine_levels',
            hue=strata,
            marker='o',
            estimator=np.mean,
            err_style='band',
            units=None,
            ci=68,
            palette=palette,
            ax=ax,
            sort=True)

        ax.set(title=feature)
        ax.tick_params(axis='both', top=False, right=False, reset =True)
        # set y axis label
        ax.set_ylabel("log$_{10}$(Cytokine levels) (pg/mL)")
        ax.set_xlabel("Day relative to CRS onset")
        ax.axvline(0, color="grey", linestyle="--")

    plt.savefig(os.path.join(figs_path,"log10_cytokines_rel_to_CRS_day10to10_fixed_y_axis.pdf"))

    days_melt.to_csv(os.path.join(figs_path,"day_relative_to_CRS_export.csv"))

    # Loop over strata ---------------------------------------------------------
    # use ASBMT CRS grading
    # colors are sometimes swapped, try running sort HLH and CRS category
    # separately

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

    days_melt = asbmt_crs.merge(
        days_melt, left_on='Randomized ID', right_on='patient_id')

    # sort HLH category
    days_melt['HLH'] = days_melt['HLH_y']
    days_melt['HLH'] = days_melt['HLH'].astype('category')
    days_melt['HLH'].cat.reorder_categories(['carHLH+', 'carHLH-'], inplace=True)
    days_melt = days_melt.sort_values(by=['HLH'])

    # Recode CRS
    days_melt['CRS grade ASBMT'] = days_melt['Max Grade CRS_ASBMT']

    days_melt['CRS grade ASBMT'].value_counts()

    # Remove patients with max_grade_CRS = 0
    days_melt = days_melt[days_melt['CRS grade ASBMT'] != 0 ]
    days_melt = days_melt.replace(
        {'CRS grade ASBMT': {1:'CRS 1-2', 2:'CRS 1-2', 3:'CRS 3-4', 4:'CRS 3-4'}})

    # sort CRS grade category
    days_melt['CRS grade ASBMT'] = days_melt['CRS grade ASBMT'].astype('category')
    days_melt['CRS grade ASBMT'].cat.reorder_categories(['CRS 1-2', 'CRS 3-4'], inplace=True)
    days_melt = days_melt.sort_values(by=['CRS grade ASBMT'])

    # Plot all cytokines
    days_melt['cytokine_levels'] = days_melt['cytokine_levels'].astype(float)
    days_melt['log10_cytokine_levels'] = np.log10(days_melt['cytokine_levels'])

    stratas = ['HLH', 'CRS grade ASBMT']
    color_palettes = [[ "#e74c3c", "#3498db"], ["orange", "green"]]
    features = ['IFN-gamma', 'IL-1B', 'IL-6', 'IL-8','IL-10',
    'IL-12p70', 'IL-13', 'MIP-alpha', 'TNF-alpha', 'IL-4',
    'IL-2', 'GM-CSF', 'IL-15', 'IL-18']

    s = stratas * len(features)
    c = color_palettes * len(features)
    f = list(np.repeat(features, len(stratas)))

    fig = plt.figure(figsize=(7, 45))
    gs = fig.add_gridspec(14, 2, hspace=0.3)
    axes = gs.subplots(sharex='col', sharey='row')

    # set font
    sns.set_style(
        {'font.family': 'sans-serif',
        'font.sans-serif': ['Arial']})

    for ax, strata, palette, feature in zip(axes.flatten(), s, c, f):

        sns.lineplot(
            data=days_melt[days_melt['cytokines']==feature],
            x='days_rel_to_CRS',
            y='log10_cytokine_levels',
            hue=strata,
            marker='o',
            estimator=np.mean,
            err_style='band',
            units=None,
            ci=68,
            palette=palette,
            ax=ax,
            sort=True)

        ax.set(title=feature)
        ax.tick_params(axis='both', top=False, right=False, reset =True)
        # set y axis label
        ax.set_ylabel("log$_{10}$(Cytokine levels) (pg/mL)")
        ax.set_xlabel("Day relative to CRS onset")
        ax.axvline(0, color="grey", linestyle="--")

    plt.savefig(os.path.join(figs_path,"log10_cytokines_rel_to_CRS_day10to10_fixed_y_axis_ASBMT_CRS.pdf"))


    # Final version for revision -----------------------------------------------
    # This version has 3 cytokines in a row.
    # Add pvalues to plot
    pvalues = pd.read_csv("revisions/data/rel_to_CRS_stats_carHLH_with_pvalues.csv"))
    fdr = pvalues[pvalues['padj_FDR'] < 0.05]

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


    # IL18 BP ------------------------------------------------------------------

    no_HLH = pd.read_excel('revisions/data/IL-18BP_No_HLH.xlsx')
    HLH = pd.read_excel('revisions/data/IL-18BP_HLH.xlsx')

    no_HLH.head()
    no_HLH['Subject_ID']

    split = no_HLH["CD22 ID"].str.split(" ", n = 1, expand = True)
    split = split.rename(columns={0:'Protocol', 1:'Subject_ID'})
    no_HLH = split.join(no_HLH)
    no_HLH = no_HLH.drop(columns=["CD22 ID"])
    no_HLH['HLH'] = 'carHLH-'

    no_HLH.head()
    no_HLH.shape

    HLH = HLH.iloc[:, 0:11]
    HLH.head()

    split = HLH["CD22 ID Number"].str.split(" ", n = 1, expand = True)
    split = split.rename(columns={0:'Protocol', 1:'Subject_ID'})
    HLH = split.join(HLH)
    HLH = HLH.drop(columns=["CD22 ID Number"])
    HLH['HLH'] = 'carHLH+'


    il18bp = pd.concat([HLH, no_HLH])
    il18bp.dtypes
    il18bp[il18bp['IL-18BP pg/nl']== '>85000']
    il18bp = il18bp.replace(to_replace = ['>85000', '>85,000'], value=85001)
    il18bp['IL-18BP pg/nl'] = pd.to_numeric(il18bp['IL-18BP pg/nl'])

    il18bp[il18bp['IL-18BP pg/nl']== 85001]

    # log10
    il18bp['log10_IL-18BP'] = np.log10(il18bp['IL-18BP pg/nl'])

    # Plot
    fig = plt.figure(figsize=(4.5,3))

    # Set style
    sns.set_style(
        {'font.family': 'sans-serif',
        'font.sans-serif': ['Arial']})

    ax = sns.lineplot(
        data=il18bp,
        x='Day Post Infusion',
        y='log10_IL-18BP',
        hue='HLH',
        marker='o',
        estimator=np.mean,
        err_style='band',
        units=None,
        ci=68,
        palette=["#e74c3c", "#3498db"],
        sort=True)

    # Set title
    ax.set_title('Serial IL-18BP Stratified by carHLH')

    # set y axis label
    ax.set_ylabel("log$_{10}$(IL-18BP levels) (pg/nl)")
    ax.set_xlabel("Day after infusion")

    fig.savefig(os.path.join(figs_path,'IL-18bp_HLH_63CI_standard_error.pdf'))

    # Plot all patients
    fig = plt.figure(figsize=(6,4))

    # Set style
    sns.set_style(
        {'font.family': 'sans-serif',
        'font.sans-serif': ['Arial']})

    ax = sns.lineplot(
        data=il18bp,
        x='Day Post Infusion',
        y='log_IL-18BP',
        hue='HLH',
        estimator=None,
        err_style=None,
        units='Subject_ID',
        palette=["#e74c3c", "#3498db"],
        sort=True)

    # Set title
    ax.set_title('Serial IL-18BP Stratified by carHLH')

    # set y axis label
    ax.set_ylabel("log$_{10}$(IL-18BP levels) (pg/nl)")
    ax.set_xlabel("Day after infusion")

    fig.savefig(os.path.join(figs_path,'IL-18bp_HLH_units.pdf',
        bbox_inches='tight'))
