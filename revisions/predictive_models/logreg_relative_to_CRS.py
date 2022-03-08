#!/usr/bin/env phython
import sys
import os
import pandas as pd
import numpy as np
import itertools as itertools
import pyhumboldt.car22 as pyh22
import pyhumboldt.ploty as pyhty
import pyhumboldt.unsupervised as pyhed
import pyhumboldt.supervised as pyhsu
from pyhumboldt.car22 import Cytokines
import missingno as msno
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import cross_val_predict
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import FunctionTransformer
from sklearn.preprocessing import OneHotEncoder
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_validate
from sklearn.metrics import recall_score
from sklearn import metrics
from sklearn.preprocessing import scale
# for nested cross validation
from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import permutation_test_score
from sklearn.preprocessing import StandardScaler
# draw AUC curves
from sklearn.metrics import auc
from sklearn.metrics import plot_roc_curve
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve
from sklearn.model_selection import StratifiedKFold
from scipy import interp
# statsmodels
import statsmodels.api as sm
# others
from sklearn.datasets import load_breast_cancer
from sklearn.datasets import load_boston
from sklearn.preprocessing import MinMaxScaler, PolynomialFeatures
from sklearn.model_selection import train_test_split
# feature selection
from sklearn.feature_selection import RFE


if __name__ == '__main__':

    # Set variables and path
    save_figs = True
    #filename_old = 'revisions/data/Full Cytokine_De-identified CD22 Data for Bioinformatics_3-30-20_v1.xlsx'
    filename = 'revisions/data/clinical_data_05-12-20_v1.xlsx'
    figs_path = 'datasets/CAR-T/new_data/' # choose path

    # Print available access variables
    print(pyh22.load_car22.__doc__)

    # Read in data -------------------------------------------------------------
    data = pyh22.load_car22(filename, access_var="all", drop=True, version=1)
    cytokines = Cytokines(data.cytokines).days_to_int()


    # Define features and target variable for baseline -------------------------
    X = pd.concat(
        [data.clinical,
        cytokines.days_to_int().df[0],
        data.pb_tbnk['0'],
        data.bm_tbnk['0'],
        data.inflammatory['0']],
        axis=1,
        sort=False)
    y = data.outcome['HLH']

    # Try building predictor for CRS+2 -----------------------------------------
    data = pyh22.load_car22(filename_old, access_var="all", drop=False, version=0)
    cytokines = pyh22.Cytokines(data.cytokines).days_to_int().crs_transform(data.outcome)

    # only keep day CRS +2
    cytokines_crs_plus2 = cytokines.crs.loc[cytokines.crs['day-post-CRS']==2].reset_index()
    cytokine_clinical = pd.merge(
        cytokines_crs_plus2,
        data.outcome.reset_index(),
        on='patient_id')

    # add BM TBNK
    cyM = pd.merge(cytokine_clinical, data.bm_tbnk['0'], on='patient_id')
    Xm = cyM.loc[:,['IFN-gamma', 'IL-1B', 'IL-6', 'TNF-alpha']].dropna()
    ym = cytokine_clinical['HLH'][Xm.index]

    X_transformed = StandardScaler().fit_transform(np.log1p(Xm))
    logreg = LogisticRegression(C=1).fit(X_transformed, ym)
    logreg.coef_
    logreg.intercept_

    #Make predictions for training set
    y_pred = logreg.predict(X_transformed)

    # Get metrics:
    pyhsu.get_metrics(ym, y_pred)

    # Plot AUC and PR
    roc_data = pyhsu._prepare_roc_curve(ym, X_transformed, logreg, annotation='AUC',
        drop_intermediate=False)

    pyhsu.ROC_curve(roc_data)
    plt.savefig(os.path.join(figs_path,"ROC_validation_seaborn_CRS2_4VAR.pdf"))

    pr_data = pyhsu._prepare_pr_curve(ym,Xm, logreg, annotation='AUC')
    pyhsu.PR_curve(pr_data)
    plt.savefig(os.path.join(figs_path,"PR_validation_seaborn_CRS2_4VAR.pdf"))

    # Try LOOCV ----------------------------------------------------------------
    data = pyh22.load_car22(filename, access_var="all", drop=False)
    cytokines = pyh22.Cytokines(data.cytokines).days_to_int().crs_transform(data.outcome)

    # Define variables to loop over
    days = [2]
    vars_1 = [
        ['BM_ratio_T_NK'],
        ['max_grade_CRS'],
        ['IFN-gamma', 'IL-1B', 'IL-6', 'TNF-alpha'],
        ['IFN-gamma', 'IL-1B', 'IL-6', 'TNF-alpha', 'BM_ratio_T_NK'],
        ['IFN-gamma', 'IL-1B', 'IL-6', 'TNF-alpha', 'max_grade_CRS'],
        ['IFN-gamma', 'IL-1B', 'IL-6', 'TNF-alpha', 'max_grade_CRS','BM_ratio_T_NK'],
        ['TCS', 'max_grade_CRS', 'BM_ratio_T_NK'],
        ['max_grade_CRS', 'BM_ratio_T_NK']]
    vars_2 = [
        ['BM_ratio_T_NK'],
        ['IFN-gamma'],
        ['IL-1B'],
        ['IL-6'],
        ['TNF-alpha'],
        ['IFN-gamma', 'IL-1B'],
        ['IFN-gamma', 'IL-6'],
        ['IFN-gamma', 'TNF-alpha'],
        ['IFN-gamma', 'BM_ratio_T_NK'],
        ['IL-1B', 'BM_ratio_T_NK'],
        ['TNF-alpha', 'BM_ratio_T_NK'],
        ['IFN-gamma', 'max_grade_CRS'],
        ['IFN-gamma', 'max_grade_CRS', 'BM_ratio_T_NK'],
        ['max_grade_CRS', 'BM_ratio_T_NK'],
        ['IFN-gamma', 'IL-1B', 'IL-6', 'TNF-alpha', 'BM_ratio_T_NK']]
    vars_3 = [
            ['BM_ratio_T_NK'],
            ['IFN-gamma', 'IL-1B', 'IL-6', 'TNF-alpha'],
            ['IFN-gamma', 'IL-1B', 'IL-6', 'TNF-alpha', 'BM_ratio_T_NK']]
    # These models were finally used
    vars = [
        ['IFN-gamma', 'BM_ratio_T_NK'],
        ['BM_ratio_T_NK'],
        ['IFN-gamma'],
        ['IL-1B'],
        ['IL-6'],
        ['TNF-alpha']]

    legend_vars_1 = ['T/NK','mgCRS','4Cyto','4Cyto+T/NK', '4Cyto+mgCRS',
        '4Cyto+T/NK+mgCRS','TCS+mgCRS+T/NK', 'mgCRS+T/NK']
    legend_vars_2 = ['T/NK', 'IFN-g', 'IL-1B', 'IL-6', 'TNF-alpha',
        'IFN-g+IL-1B', 'IFN-g+IL-6', 'IFN-g+TNF-alpha', 'IFN-g+T/NK',
        'IL-1B+T/NK', 'TNF-alpha+T/NK', 'IFN-g+mgCRS', 'IFN-g+mgCRS+T/NK',
        'mgCRS+T/NK','4Cyto+T/NK']
    legend_vars_3 = ['T/NK','4Cyto','4Cyto+T/NK']
    legend_vars = ['IFN-g+T/NK', 'T/NK', 'IFN-g', 'IL-1B', 'IL-6', 'TNF-alpha']

    legend_days = ['CRS2']

    roc_data = pd.DataFrame()
    pr_data = pd.DataFrame()

    # Define classifier and CV
    clf = LogisticRegression(C=1)
    cv = LeaveOneOut()

    # for testing purposes
    #day = 1
    #var = ['IFN-gamma', 'IL-1B', 'IL-6', 'TNF-alpha',
    #    'BM_ratio_T_NK', 'HLH']
    var = ['IFN-gamma', 'BM_ratio_T_NK']
    legend_var = ['IFN-g+T/NK']
    legend_day = ['CRS2']

    for day, legend_day in zip(days, legend_days):

        print(day)
        print(legend_day)

        cytokines_crs_plus2 = cytokines.crs.loc[
            cytokines.crs['day-post-CRS'] == day].reset_index()
        cyto_clin = pd.merge(
            cytokines_crs_plus2,
            data.clinical.reset_index(),
            on='patient_id')

        # add BM TBNK
        cyM = pd.merge(cyto_clin, data.bm_tbnk['0'], on='patient_id')

        # log transform continous variables
        cyM['BM_ratio_T_NK'] = np.log(cyM['BM_ratio_T_NK'])
        cyM['IFN-gamma'] = np.log(cyM['IFN-gamma'])
        cyM['IL-1B'] = np.log(cyM['IL-1B'])
        cyM['IL-6'] = np.log(cyM['IL-6'])
        cyM['TNF-alpha'] = np.log(cyM['TNF-alpha'])

        for var, legend_var in zip(vars, legend_vars):
            print(var)
            print(legend_var)

            Xm = cyM.loc[:,var].dropna()

            ym = cyM['HLH'][Xm.index]

            print("Shape of Xm, ym:")
            print(Xm.shape)
            print(ym.shape)

            #Xm_transformed = StandardScaler().fit_transform(Xm)
            Xm_transformed = Xm.copy()

            print("Shape of Xm_transformed:")
            print(Xm_transformed.shape)

            X_scale = pd.DataFrame(
                data=Xm_transformed,
                columns=Xm.columns)

            print("Shape of X_scale, ym:")
            print(X_scale.shape)
            print(ym.shape)

            # Leave one out cross validation
            all_y = []
            all_probs = []
            for train_index, test_index in cv.split(X_scale, ym):
                #print(train_index)
                #print(test_index)
                Xm_train = X_scale.iloc[train_index]
                ym_train = ym.iloc[train_index]
                Xm_test = X_scale.iloc[test_index]
                ym_test = ym.iloc[test_index]
                all_y.append(ym_test.to_list())
                all_probs.append(list(clf.fit(Xm_train, ym_train).decision_function(Xm_test)))

            all_y = [num for elem in all_y for num in elem]
            all_scores = [num for elem in all_probs for num in elem]

            len(ym)
            len(all_scores)
            len(fpr)
            fpr, tpr, thresholds = roc_curve(
                ym, all_scores, drop_intermediate=False)

            annotation = legend_day + '_' + legend_var
            print(annotation)
            #annotation = 'INF-g+T/NK'

            pd_roc = pd.concat(
                [pd.Series(fpr, name='fpr'),
                pd.Series(tpr, name='tpr'),
                pd.Series(thresholds, name='thresholds'),
                pd.Series(list(itertools.repeat(annotation, len(fpr))),
                    name='Legend')],
                axis=1)

            roc_data = roc_data.append(pd_roc)

            precision, recall, thresholds = precision_recall_curve(
                ym, all_scores)
            # Merge results to dataframe for plotting
            pd_pr = pd.concat(
                [pd.Series(recall, name='recall'),
                pd.Series(precision, name='precision'),
                pd.Series(thresholds, name='thresholds'),
                pd.Series(list(itertools.repeat(annotation, len(recall))),
                name='Legend')],
                axis=1)
            pr_data = pr_data.append(pd_pr)

    pyhsu.ROC_curve(
        roc_data,
        alpha=0.8,
        figsize=(5.5, 5.5),
        add_defaults=False,
        style='Legend')
    plt.savefig(os.path.join(figs_path,"4Cyto_ROC_cross-validation_seaborn_model_comp_day2.pdf"))

    pyhsu.PR_curve(
        pr_data,
        add_defaults=False,
        figsize=(5.5, 5.5))
    plt.savefig(os.path.join(figs_path,"4Cyto_PR_cross-validation_seaborn_model_comp_day2.pdf"))

    # Scatterplot --------------------------------------------------------------
    import seaborn as sns
    X_plot = Xm.replace({'HLH': {0:'carHLH-', 1:'carHLH+'}})
    sns.pairplot(X_plot, hue='HLH')
    plt.savefig(os.path.join(figs_path,"pairplot_HLH.pdf"))

    sns.pairplot(X_plot, kind='reg')
    plt.savefig(os.path.join(figs_path,"pairplot_linreg.pdf"))

    # Get values for 2 variable model ------------------------------------------
    # These models were finally used
    var = ['IFN-gamma', 'BM_ratio_T_NK']
    var2 = ['IFN-gamma', 'BM_ratio_T_NK', 'HLH', 'patient_id']

    cytokines_crs_plus2 = cytokines.crs.loc[
        cytokines.crs['day-post-CRS'] == 2].reset_index()
    cyto_clin = pd.merge(
        cytokines_crs_plus2,
        data.clinical.reset_index(),
        on='patient_id')

    # add BM TBNK
    cyM = pd.merge(cyto_clin, data.bm_tbnk['0'], on='patient_id')

    # log transform continous variables
    cyM['BM_ratio_T_NK'] = np.log(cyM['BM_ratio_T_NK'])
    cyM['IFN-gamma'] = np.log(cyM['IFN-gamma'])
    cyM['IL-1B'] = np.log(cyM['IL-1B'])
    cyM['IL-6'] = np.log(cyM['IL-6'])
    cyM['TNF-alpha'] = np.log(cyM['TNF-alpha'])

    Xm = cyM.loc[:,var].dropna()
    Xm_HLH = cyM.loc[:,var2].dropna()
    ym = cyM['HLH'][Xm.index]

    clf = LogisticRegression(C=1).fit(Xm, ym)
    scores = clf.decision_function(Xm)
    clf.coef_
    clf.intercept_

    pd_roc['youden_index'] = pd_roc['tpr'] + (1-pd_roc['fpr'])
    max(pd_roc['youden_index'])

    # Waterfall plot for 2 variable model --------------------------------------
    pd_roc['HLH'] =
    Xm_HLH['scores'] = scores
    Xm_HLH['rank'] = Xm_HLH['scores'].rank(ascending=False)
    Xm_HLH = Xm_HLH.replace({'HLH': {0:'carHLH-', 1:'carHLH+'}})

    fig, ax = plt.subplots(figsize=(5.5,3.5))
    sns.set_context("notebook")
    sns.set_style(
        {'font.family': 'sans-serif',
        'font.sans-serif': ['Arial']})
    sns.barplot(
        data=Xm_HLH,
        x='rank',
        y='scores',
        hue='HLH',
        palette=["#e74c3c", "#3498db"], dodge=False)
    ax.xaxis.set_ticklabels([])
    ax.set_xticks([])
    plt.ylim([-6,6])
    plt.ylabel('Classifier value')
    plt.xlabel('Patients')
    plt.savefig(os.path.join(figs_path,"waterfall_plot_2var_model.pdf"))

    # LOOCV waterfall plot -----------------------------------------------------
    pd_roc = pd.concat(
        [pd.Series(all_y, name='all_y'),
        pd.Series(all_scores, name='all_scores'),
        pd.Series(list(itertools.repeat(annotation, len(all_y))),
            name='Legend')],
        axis=1)

    pd_roc['rank'] = pd_roc['all_scores'].rank(ascending=False)
    pd_roc = pd_roc.replace({'all_y': {0:'carHLH-', 1:'carHLH+'}})

    pd_roc = pd_roc.rename(columns={'all_y':'HLH'})

    fig, ax = plt.subplots(figsize=(5.5,3.5))
    sns.set_context("notebook")
    sns.set_style(
        {'font.family': 'sans-serif',
        'font.sans-serif': ['Arial']})
    sns.barplot(
        data=pd_roc,
        x='rank',
        y='all_scores',
        hue='HLH',
        palette=["#e74c3c", "#3498db"], dodge=False)
    ax.xaxis.set_ticklabels([])
    ax.set_xticks([])
    plt.ylim([-6,6])
    plt.ylabel('Classifier value')
    plt.xlabel('Patients')

    plt.savefig(os.path.join(figs_path,"LOOCV_waterfall_plot_2var_model.pdf"))


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
    cytokines_p2 = days_outcome.loc[days_outcome['days_rel_to_CRS'].isin(range(-10, 10))].reset_index()
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

    # Select specific cytokines
    days_melt = days_melt.loc[days_melt['cytokines'].isin(['IFN-gamma',
        'IL-6', 'TNF-alpha', 'IL-1B'])]
    days_melt['cytokine_levels'] = days_melt['cytokine_levels'].astype(float)
    days_melt['log_cytokine_levels'] = np.log(days_melt['cytokine_levels'])

    # Recode HLH column
    days_melt = days_melt.replace(
        {'HLH': {0:'carHLH-', 1:'carHLH+'}})

    # Plot patients over time
    axes = pyhty.features_over_time(
        rows=1,
        columns=4,
        features=days_melt,
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
        xlabel="Days relative to CRS onset",
        ylabel="log(Cytokine levels) (pg/mL)",
        palette=["#3498db", "#e74c3c"],
        draw_line=True)

    plt.savefig(os.path.join(figs_path,"4_logcytokines_rel_to_CRS_day10to10.pdf"))

    # Predictive modeling ------------------------------------------------------
    # model selection
    # Single cytokine model
    models2 = days_melt['cytokines'].unique()
    days = [0, 1, 2, 3]
    folds = 3
    print_plot=True

    df_model = []
    df_day = []
    df_numpat = []
    df_meanROC = []
    df_sdROC = []

    # set specific model and cytokine
    #day = 2
    #model = ['IFN-gamma', 'IL-1B', 'IL-6', 'TNF-alpha', 'BM_ratio_T_NK']

    for model in models2:

        for day in days:

            df_model.append(model)
            df_day.append(day)

            # select specific days
            cytokines_p2 = days_outcome.loc[
                days_outcome['days_rel_to_CRS'].isin([day])].reset_index()

            # add HLH info
            Xclin = pd.merge(
                cytokines_p2,
                data.primary_outcome.reset_index(),
                on='patient_id')

            # add BM TBNK
            cyM = pd.merge(
                Xclin,
                data.bm_tbnk['0'].reset_index(),
                on='patient_id')

            # Select variables
            Xm = cyM.loc[:,model].dropna()

            # subset
            #Xm = Xm.loc[:,['IFN-gamma', 'BM_ratio_T_NK']]

            # Number of patients
            num_patients = Xm.shape[0]
            df_numpat.append(Xm.shape[0])
            print("Number of patients: ", num_patients)

            # data transformation
            Xm = np.log(Xm)
            #Xm['BM_ratio_T_NK'] = np.log(Xm['BM_ratio_T_NK'])
            #Xm['IFN-gamma'] = np.log(Xm['IFN-gamma'].astype(float))
            #Xm['IL-1B'] = np.log(Xm['IL-1B'].astype(float))
            #Xm['IL-6'] = np.log(Xm['IL-6'].astype(float))
            #Xm['TNF-alpha'] = np.log(Xm['TNF-alpha'].astype(float))

            # select ym
            ym = cyM['HLH'][Xm.index]

            # With ROC curve!
            cv = StratifiedKFold(
                n_splits=folds,
                shuffle=True,
                random_state=0)

            # Collect data for ROC curve
            roc_data = pd.DataFrame()
            #target_var = 'HLH'

            if(print_plot):
                for i, (train_index, test_index) in enumerate(cv.split(Xm, ym)):

                    Xm_train = Xm.iloc[train_index]
                    ym_train = ym.iloc[train_index]
                    Xm_test = Xm.iloc[test_index]
                    ym_test = ym.iloc[test_index]

                    #print(Xm_train)
                    #print(Xm_test)
                    if(single):
                        Xtr = np.array(Xm_train).reshape(-1,1)
                        Xte = np.array(Xm_test).reshape(-1,1)
                    else:
                        Xtr = Xm_train.copy()
                        Xte = Xm_test.copy()

                    clf = LogisticRegression(
                        random_state=0, C=1).fit(
                            Xtr,
                            ym_train)

                    roc_data = roc_data.append(
                            pyhsu._prepare_roc_curve(
                                ym_test,
                                Xte,
                                clf,
                                annotation='ROC fold {}'.format(i),
                                drop_intermediate=False))

                rocs = list(roc_data['Legend'].unique())
                roc_auc = []
                for roc in rocs:
                    roc_sel = roc_data.loc[roc_data['Legend']==roc]
                    roc_auc.append(metrics.auc(roc_sel['fpr'], roc_sel['tpr']))

                df_meanROC.append(np.mean(roc_auc))
                df_sdROC.append(np.std(roc_auc))

                #pyhsu.ROC_curve(
                #    roc_data,
                #    alpha=0.5,
                #    style='Legend',
                #    lw=2,
                #    mean_roc=True,
                #    sd=True)

                #plt.savefig(os.path.join(
                #            figs_path,
                #            'fair_IFN-gamma+BM_ratio_T_NK_3fold_shuffle=True.pdf'))


                #plt.savefig(os.path.join(
                #    figs_path,
                #    '%s_day%s_cv%s_fold_shuffle=True.pdf' % (model, str(day), str(folds))))

    df_t = pd.DataFrame(
        list(zip(df_model, df_day, df_numpat, df_meanROC, df_sdROC)),
        columns = ['model', 'day', 'num_patients', 'meanROC', 'sdROC'])

    # Plot AUCs for cytokines
    fig, ax = plt.subplots(figsize=(8,3.5))
    sns.set_style(
        {'font.family': 'sans-serif',
        'font.sans-serif': ['Arial']})
    sns.barplot(
        data=df_t,
        x='day',
        y='meanROC',
        hue='model',
        palette = sns.color_palette("muted", n_colors=14))
    plt.ylim([0, 1])
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.ylabel('mean AUC (3-fold CV)')
    plt.xlabel('Days relative to CRS onset')

    plt.savefig(os.path.join(
            figs_path,
            'barplot_all_cytokines_3fold_shuffle=True.pdf'))

    # Run permutation test -----------------------------------------------------
    Xm.shape
    n_classes = np.unique(ym).size

    clf = LogisticRegression(random_state=0, C=1)
    cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=0)

    score, permutation_scores, pvalue = permutation_test_score(
        clf, Xm, ym,
        scoring="roc_auc",
        cv=cv,
        n_permutations=100,
        n_jobs=1)

    # View histogram of permutation scores
    fig, ax = plt.subplots(figsize=(4,4))
    sns.set_style(
        {'font.family': 'sans-serif',
        'font.sans-serif': ['Arial']})
    sns.distplot(
        permutation_scores,
        hist_kws=dict(edgecolor="white", linewidth=2),
        label='Permutation score')
    ylim = plt.ylim()
    ax.plot(2 * [score], ylim, '--g', linewidth=1,
         label='Mean AUC score (IFN-g+BM_T/NK)'
         ' (p-value %0.4f)' % pvalue)
    ax.plot(2 * [1. / n_classes], ylim, '--k', linewidth=1, label='Random')
    ax.legend(loc='upper left')
    plt.xlim([0, 1])
    plt.xlabel('AUC scores')

    plt.savefig(os.path.join(figs_path,"2varmodel_3fold_permutation.pdf"))


    # Get coefficients and intercept -------------------------------------------

    logreg = LogisticRegression(C=1).fit(Xm, ym)
    logreg.coef_
    logreg.intercept_

    #Make predictions for training set
    y_pred = logreg.predict(Xm)

    # Get metrics:
    pyhsu.get_metrics(ym, y_pred)

    # Plot AUC and PR
    roc_data = pyhsu._prepare_roc_curve(ym, Xm, logreg, annotation='AUC',
        drop_intermediate=False)

    pyhsu.ROC_curve(roc_data)

    # For manuscript -----------------------------------------------------------
    model1 = ['IFN-gamma', 'IL-1B', 'IL-6', 'TNF-alpha']
    model2 = ['IFN-gamma', 'BM_ratio_T_NK']
    model3 = ['BM_ratio_T_NK']
    model4 = ['IFN-gamma']
    day = 2

    # select specific days
    cytokines_p2 = days_outcome.loc[
        days_outcome['days_rel_to_CRS'].isin([day])].reset_index()

    # add HLH info
    Xclin = pd.merge(
        cytokines_p2,
        data.primary_outcome.reset_index(),
        on='patient_id')

    # add BM TBNK
    cyM = pd.merge(
        Xclin,
        data.bm_tbnk['0'].reset_index(),
        on='patient_id')

    # Select variables
    Xm = cyM.loc[:,model2].dropna()
    #Xm = Xm.loc[:,model3]
    Xm = np.log(Xm)
    ym = cyM['HLH'][Xm.index]

    clf = LogisticRegression(C = 10e9).fit(Xm, ym)
    clf.coef_
    clf.intercept_
    # Print metrics
    pyhsu.get_metrics(ym, clf.predict(Xm))

    # cross-validation
    roc_data = pd.DataFrame()
    sensitivity = []
    specificity = []

    # set cross validation
    cv = StratifiedKFold(
        n_splits=3,
        shuffle=True,
        random_state=0)

    for i, (train_index, test_index) in enumerate(cv.split(Xm, ym)):

        Xm_train = Xm.iloc[train_index]
        ym_train = ym.iloc[train_index]
        Xm_test = Xm.iloc[test_index]
        ym_test = ym.iloc[test_index]

        clf = LogisticRegression(C=10e9).fit(Xm_train, ym_train)
        print('ROC fold:{}'.format(i))
        cm = metrics.confusion_matrix(ym_test, clf.predict(Xm_test))
        sensitivity.append(cm[1,1]/(cm[1,1] + cm[1,0]))
        specificity.append(cm[0,0]/(cm[0,0] + cm[0,1]))
        print(sensitivity)
        print(specificity)

        roc_data = roc_data.append(
            pyhsu._prepare_roc_curve(
                ym_test, Xm_test, clf,
                annotation='ROC fold {}'.format(i),
                drop_intermediate=False))

    # get mean sensitivity & specificity
    np.mean(sensitivity)
    np.std(sensitivity)
    np.mean(specificity)
    np.std(specificity)

    pyhsu.ROC_curve(roc_data, alpha=0.5, style='Legend', lw=2, mean_roc=True, sd=True)
    plt.savefig(os.path.join(figs_path,"BM_ratio_T_NK=39_3fold.pdf"))

    # Same with sm.Logit
    Xm['intercept'] = 1
    logit = sm.Logit(ym, Xm[model2 + ['intercept']])
    results = logit.fit()
    print(results.summary())

    # Get coefficients and intercept
    clf = LogisticRegression(C = 10e9).fit(Xm, ym)
    clf.coef_
    clf.intercept_

    fpr, tpr, threshold = roc_curve(
            ym, clf.decision_function(Xm), drop_intermediate=False)

     pd_roc = pd.concat(
        [pd.Series(fpr, name='fpr'),
        pd.Series(tpr, name='tpr'),
        pd.Series(threshold, name='thresholds')],
        axis=1)

    pd_roc['youden_index'] = pd_roc['tpr'] + (1-pd_roc['fpr'])
    max(pd_roc['youden_index'])
