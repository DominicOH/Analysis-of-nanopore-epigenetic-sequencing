import pandas as pd
from scipy import stats
from sklearn import metrics 
import numpy as np
from concurrent import futures
from AnalysisTools import read_modbed

def onlyAutosomal(df):
    df = df.loc[df.loc[:, "Chromosome"].str.match("^(chr)\d+$")]
    return df

def fetch_data(dry_run: bool, **kwargs):
    if not dry_run:
        cbm2_path = "data/modbases/nanopore/cbm2/"
        cbm3_path = "data/modbases/nanopore/cbm3/"

        oxbs_path = "data/modbases/public/CRR008808_oxBS/masked/"
        tab_path = "data/modbases/public/CRR008807_TAB/masked/"
    else:
        cbm2_path = "data/dryruns/cbm2/"
        cbm3_path = "data/dryruns/cbm3/" 

        oxbs_path = "data/dryruns/oxbs/"
        tab_path = "data/dryruns/tab/"     

    cbm2_auto = onlyAutosomal(read_modbed.openReps(cbm2_path, **kwargs))
    cbm3_auto = onlyAutosomal(read_modbed.openReps(cbm3_path, **kwargs))

    oxbs_auto = onlyAutosomal(read_modbed.openReps(oxbs_path, modbase="5mC", **kwargs))
    tab_auto = onlyAutosomal(read_modbed.openReps(tab_path, modbase="5hmC", **kwargs))

    return cbm2_auto, cbm3_auto, tab_auto, oxbs_auto

def fetch_data_Parallel(dry_run: bool, split_biorep=False, **kwargs):
    """
    Fetches modkit/bismark data. Order of return is: cbm2, cbm3, tab, oxbs. 

    ::bool dry_run:: Whether to return data for internal testing.  
    """

    if not dry_run:
        cbm2_path = "data/modbases/nanopore/cbm2/"
        cbm3_path = "data/modbases/nanopore/cbm3/"

        oxbs_path = "data/modbases/public/CRR008808_oxBS/masked/"
        tab_path = "data/modbases/public/CRR008807_TAB/masked/"
    else:
        cbm2_path = "data/dryruns/cbm2/"
        cbm3_path = "data/dryruns/cbm3/" 

        oxbs_path = "data/dryruns/oxbs/"
        tab_path = "data/dryruns/tab/"

    if split_biorep:
        bio_reps = ["Nanopore 1", "Nanopore 2"]
    else:
        bio_reps = ["Nanopore", "Nanopore"]

    with futures.ThreadPoolExecutor(4) as ppe:
        all_futures = [ppe.submit(read_modbed.openReps_Parallel, path, insert_cols={"Technique" : bio_rep}, **kwargs) for path, bio_rep in zip([cbm2_path, cbm3_path], bio_reps)]
        all_futures.append(ppe.submit(read_modbed.openReps_Parallel, tab_path, insert_cols={"Technique" : "TAB"}, modbase="5hmC", **kwargs))
        all_futures.append(ppe.submit(read_modbed.openReps_Parallel, oxbs_path, insert_cols={"Technique" : "oxBS"},  modbase="5mC", **kwargs))
        future_dfs = [onlyAutosomal(future.result()) for future in all_futures]

    return future_dfs

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
