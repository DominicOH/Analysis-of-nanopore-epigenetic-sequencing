import pandas as pd
from math import sqrt
import pyranges as pr
import numpy as np
from scipy import stats
from sklearn import metrics 
from read_modbed import readBismarkZeroCov, readModkit

def compareStats(x, y):
    """
    Compares two Series and outputs a Series of summary statistics. 
    """
    r = stats.pearsonr(x, y)
    p = stats.spearmanr(x, y)
    ks = stats.ks_2samp(x, y)
    mw = stats.mannwhitneyu(x, y)
    cvm = stats.cramervonmises_2samp(x, y)
    ad = stats.anderson_ksamp([x, y])
    t_test = stats.ttest_ind(x, y)

    d = {
        "Pearson_r" : round(r.statistic, 2),
        "Spearman_p" : round(p.statistic, 2),
        "Correlation_p_values" : [round(r.pvalue, 2), round(p.pvalue, 2)],
        "Kendall" : round(stats.kendalltau(x, y).statistic, 2), 
        "RMSE" : metrics.mean_squared_error(x, y, squared=False),
        "Mean Absolute Error" : metrics.mean_absolute_error(x, y),
        "Median Absolute Error" : metrics.median_absolute_error(x, y),
        "KS" : [ks.statistic, ks.pvalue],
        "MW" : [mw.statistic, mw.pvalue],
        "CVM" : [cvm.statistic, cvm.pvalue],
        "AD" : [ad.statistic, ad.pvalue],
        "T-Test" : [t_test.statistic, t_test.pvalue] 
    }
    return pd.Series(d)

