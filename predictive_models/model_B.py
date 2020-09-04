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
    filename = snakemake.input[0]

    # Print available access variables
    # print(pyh22.load_car22.__doc__)

    # Read in data -------------------------------------------------------------
    data = pyh22.load_car22(filename, access_var="all", drop=True, version=1)
    cytokines = Cytokines(data.cytokines).days_to_int()

    # Plot missing data
    #msno.matrix(cytokines.df)

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

    # 2 parameter model
    model = ['IFN-gamma', 'BM_ratio_T_NK']
    day = 2

    # select specific days
    cytokines_p2 = days_outcome.loc[days_outcome['days_rel_to_CRS'].isin([day])].reset_index()

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
    Xm = np.log(Xm)
    ym = cyM['HLH'][Xm.index]

    # Cross-validation
    roc_data = pd.DataFrame()
    sensitivity = []
    specificity = []

    # Initiate cross validation
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
        print('Sensitivity:', sensitivity)
        print('Specificity:', specificity)

        roc_data = roc_data.append(
            pyhsu._prepare_roc_curve(
                ym_test, Xm_test, clf,
                annotation='ROC fold {}'.format(i),
                drop_intermediate=False))

    # Get mean sensitivity & specificity
    print('Mean sensitivity:', np.mean(sensitivity))
    print('Std sensitivity:', np.std(sensitivity))
    print('Mean specificity:', np.mean(specificity))
    print('Std specificity', np.std(specificity))

    # Print ROC curve
    print('AUC for each fold:')
    pyhsu.ROC_curve(roc_data, alpha=0.5, style='Legend', lw=2, mean_roc=True, sd=True)
    plt.savefig(snakemake.output[1])

    # Same with sm.Logit
    Xm['intercept'] = 1
    logit = sm.Logit(ym, Xm[model + ['intercept']])
    results = logit.fit()

    # Write results to file
    f = open(snakemake.output[0], 'w')
    f.write(str(results.summary()))
    f.close()

    # Same with sklearn logistic regression
    # Get coefficients and intercept
    clf = LogisticRegression(C = 10e9).fit(Xm, ym)
    #clf.coef_
    #clf.intercept_

    fpr, tpr, threshold = roc_curve(
            ym, clf.decision_function(Xm), drop_intermediate=False)

    pd_roc = pd.concat(
        [pd.Series(fpr, name='fpr'),
        pd.Series(tpr, name='tpr'),
        pd.Series(threshold, name='thresholds')],
        axis=1)

    pd_roc['youden_index'] = pd_roc['tpr'] + (1-pd_roc['fpr'])
    print('Youden index:', max(pd_roc['youden_index']))
