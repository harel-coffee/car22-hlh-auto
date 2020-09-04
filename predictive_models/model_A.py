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

    # Visualize features in histogram ------------------------------------------
    #pyhsu.feature_hist(X, y_target='HLH')

    # Three variable model -----------------------------------------------------
    # For the three variable model, no log tranformation was performed!
    # To generate Model A, following an initial screening by univariate methods,
    # for those parameters for which p<0.05, multivariable logistic regression
    # analysis using both backward and stepwise selection was used to identify a
    # set of factors which could jointly impact development of carHLH. Because
    # all but one of the carHLH occurrences were after day 7, prediction models
    # concentrated on factors which were known at or before day 7, with
    # exception for the maximum grade of CRS. Patients with missing values were
    # excluded from analysis.
    # Backward and stepwise selection was performed using SAS.
    model_features = ['TCS', 'max_grade_CRS', 'BM_ratio_T_NK']
    Xm = X.loc[:, model_features].dropna()
    ym = y[Xm.index]

    # Fit model using sklearn
    # set C value very high --> no regularization
    clf = LogisticRegression(C = 10e9).fit(Xm, ym)
    # clf.coef_
    # clf.intercept_
    # Print metrics
    # pyhsu.get_metrics(ym, clf.predict(Xm))

    # Fit model using sm.Logit
    model_features = ['TCS', 'max_grade_CRS', 'BM_ratio_T_NK']
    Xm = X.loc[:, model_features].dropna()
    ym = y[Xm.index]

    # Add intercept
    Xm['intercept'] = 1
    logit = sm.Logit(ym, Xm[model_features + ['intercept']])
    results = logit.fit()

    # Write results to file
    f = open(snakemake.output[0], 'w')
    f.write(str(results.summary()))
    f.close()

    # Get metrics
    auc = roc_auc_score(
        ym,
        results.predict(Xm[model_features + ['intercept']]))
    print("AUC:", auc)

    # Get values for ROC curve
    fpr, tpr, threshold = roc_curve(
        ym, results.predict(Xm[model_features + ['intercept']]))

    # Cross validation
    cv = StratifiedKFold(
        n_splits=4,
        shuffle=True,
        random_state=0)

    crossval_scores = cross_val_score(
        LogisticRegression(C=10e9),
        Xm[model_features],
        ym,
        scoring='roc_auc',
        cv=cv)

    # print(crossval_scores, crossval_scores.mean())

    # Collect data for ROC curve
    roc_data = pd.DataFrame()
    sensitivity = []
    specificity = []
    target_var = 'HLH'

    # Initiate cross validation
    cv = StratifiedKFold(
        n_splits=4,
        shuffle=True,
        random_state=0)

    # For each fold...
    for i, (train_index, test_index) in enumerate(cv.split(Xm, ym)):

        Xm_train = Xm.reset_index().loc[train_index].drop(columns='patient_id')
        ym_train = ym.reset_index().loc[train_index].drop(columns='patient_id')
        Xm_test = Xm.reset_index().loc[test_index].drop(columns='patient_id')
        ym_test = ym.reset_index().loc[test_index].drop(columns='patient_id')

        clf = LogisticRegression(C=10e9).fit(Xm_train, ym_train[target_var])
        print('ROC fold:{}'.format(i))
        cm = metrics.confusion_matrix(ym_test[target_var], clf.predict(Xm_test))
        sensitivity.append(cm[1,1]/(cm[1,1] + cm[1,0]))
        specificity.append(cm[0,0]/(cm[0,0] + cm[0,1]))
        print('Sensitivity:', sensitivity)
        print('Specificity:', specificity)

        roc_data = roc_data.append(
            pyhsu._prepare_roc_curve(
                ym_test[target_var], Xm_test, clf,
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
