import pandas as pd
from math import sqrt
import pyranges as pr
import numpy as np
from scipy import stats
from sklearn import metrics 

def filterDepth(df, 
                min_depth: int = 10, 
                apply_max_depth: bool = False):
    """
    Filters the dataframe to only rows within the minimum and maximum coverage depth. 
    """
    filtered_df = df.copy()
    average = filtered_df["readCount"].mean()
    filtered_df = filtered_df.loc[filtered_df.loc[:, "readCount"] >= min_depth]

    if apply_max_depth:
        filtered_df = filtered_df.loc[filtered_df.loc[:, "readCount"] < (average + 3*sqrt(average))]

    return filtered_df

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

def readModkit(
        path: str, 
        min_depth: int = 10,
        apply_max_depth: bool = True,
        incl_raw_counts: bool = False
):
    """
    Reads the bedmMethyl output of Modkit pileup into a pd.DataFrame. 
    NOTE: It's important that modkit pileup is run with the --only-tabs flag. Otherwise important data columns are separated only by spaces and missed by tab-parsing. 

    :param str path: Filepath of the bedMethyl file. 
    :param int min_depth: The minimum readcount of CpG sites.
    :param bool apply_max_depth: Whether to filter out modbases with a depth greater than d + 3*sqrt(d); where d is the mean depth.
    :param bool incl_raw_counts: Whether the raw count of modified basecalls should be kept in the resulting dataframe.

    """
    colnames = ["Chromosome", "Start", "End", "modBase", "modScore", "Strand", "rem1", "rem2", "rem3", "readCount", "percentMeth", "N_mod", "N_canonical", "N_other", "N_delete", "N_fail", "N_diff", "N_nocall"]
    df_init = pd.read_csv(path, 
                          sep="\t", 
                          names=colnames)
    
    df_filtered = filterDepth(df_init, min_depth, apply_max_depth)

    if incl_raw_counts:
        pivot_cols =  ["readCount", "percentMeth", "N_mod", "N_canonical", "N_other"]
    else:
        pivot_cols = ["readCount", "percentMeth"]

    df_pivot = df_filtered.pivot(index=["Chromosome", "Start", "End", "Strand"], columns="modBase", values=pivot_cols)
    df_pivot.columns = df_pivot.columns.to_flat_index()
    df_pivot = df_pivot.reset_index()
    
    df_pivot = df_pivot.rename(columns={
        ('readCount', 'h') : "readCount",
        "h" : "percentMeth_5hmC",
        "m" : "percentMeth_5mC",
        ("percentMeth", "h") : "percentMeth_5hmC",
        ("percentMeth", "m") : "percentMeth_5mC",
        ("N_mod", "h") : "N_hmC",
        ("N_mod", "m") : "N_mC",
        ("N_canonical", "h") : "N_C"
    }, errors="ignore")

    df_pivot = df_pivot.drop(columns=[
        ("readCount", "m"),
        ("N_canonical", "m"),
        ("N_other", "h"),
        ("N_other", "m")
    ], errors="ignore")

    return df_pivot 

def readBismarkZeroCov(
        path: str, 
        mod: str, 
        min_depth: int = 10,
        apply_max_depth: bool = False,
        incl_raw_counts: bool = False):
    """
    Reads the output file of Bismark methylation extractor. Requires the bed format output produced with the options: -p --bedGraph --zero_based --comprehensive

    :param str path: Path to the bed file of modified bases.
    :param str mod: Type of target modification ["5mC", "5hmC"] 
    :param bool filter: Whether to filter the dataframe by depth [default=True]
    :return: Pandas DataFrame object 

    """
    if mod == "5mC":
        meth_col = "percentMeth_5mC"
    elif mod == "5hmC":
        meth_col = "percentMeth_5hmC"
    else:
        raise ValueError("Please enter a mod type: '5mC' or '5hmC'")

    df = pd.read_csv(path, sep="\t", names=[
            "Chromosome", "Start", "End", meth_col, "N_mod", "N_unmod"]
            ).assign(readCount = lambda row: row.N_mod + row.N_unmod)
        
    if min_depth:
        df = filterDepth(df, min_depth, apply_max_depth)

    if incl_raw_counts:
        return df
    else: 
        return df.drop(columns=["N_mod", "N_unmod"])
    
def pieData(annotated_peak_data, 
            count_column: str,
            min_depth: int = 0, 
            max_depth: int = np.inf
            ): 
    """
    Small function that outputs labels and values for a matplotlib pie chart. Requires the feature_type column. 
    Optional min/max values filter out peaks of below minimum/above maximum depth. 

    :returns:
        - pie_values - Values for the pie chart. 
        - pie_labels - Labels for each value in the pie chart. 
    """
    if min_depth > 0 or max_depth < np.inf:
        try: 
            annotated_peak_data = annotated_peak_data.query(f"total_peakDepth >= {min_depth} & total_peakDepth <= {max_depth}")
        except:
            print("Warning: Not filtering peak depth; not selected or not possible.")

    pie_labels = annotated_peak_data[f"{count_column}"].value_counts().index
    pie_values = annotated_peak_data[f"{count_column}"].value_counts().values

    return pie_values, pie_labels     

def readModbam2bed(path: str,
                   min_depth: int = 10, 
                   apply_max_depth: bool = False,
                   incl_raw_counts: bool = False):
    """
    Opens Modbam2bed bedMethyl files in an appropriate format for this analysis. 
    
    Note: Requires modbam2bed extended output with options: -e --cpg -m 5mC (other mod is assumed to be 5hmC)
    """
    modbed = pd.read_csv(path, sep="\t", 
                names=["Chromosome", "Start", "End", "mod_type", "score", "Strand", "i1", "i2", "i3", "readCount", "percentMeth_mC", "N_C", "N_mC", "N_filt", "N_NA", "N_hmC"])
    modbed["readCount_T"] = modbed.loc[:, ("N_C", "N_mC", "N_hmC")].sum(axis="columns")

    modbed.drop(columns="readCount", inplace=True)
    modbed.rename(columns={"readCount_T" : "readCount"}, inplace=True)

    if min_depth:
        modbed = filterDepth(modbed, min_depth, apply_max_depth)

    modbed["percentMeth_C"] = modbed.loc[:, "N_C"].divide(modbed.loc[:, "readCount"]).multiply(100)
    modbed["percentMeth_5hmC"] = modbed.loc[:, "N_hmC"].divide(modbed.loc[:, "readCount"]).multiply(100)
    modbed["percentMeth_5mC"] = modbed.loc[:, "N_mC"].divide(modbed.loc[:, "readCount"]).multiply(100)

    if incl_raw_counts:
        return modbed.loc[:, ("Chromosome", "Start", "End", "Strand", "readCount", "N_C", "N_mC", "N_hmC", "percentMeth_C", "percentMeth_5mC", "percentMeth_5hmC")]
    else: 
        return modbed.loc[:, ("Chromosome", "Start", "End", "Strand", "readCount", "percentMeth_C", "percentMeth_5mC", "percentMeth_5hmC")]

def asPyRanges(df):
    """
    Function to change pandas DataFrame colnames for PyRanges compatibility. 
    """
    print("Changing colnames to be PyRanges compatible...")
    try:
        df = df.rename(columns={
            "chromosome" : "Chromosome",
            "chromStart" : "Start",
            "chromEnd" : "End"
        }, errors="ignore")
        print("Done")
        return pr.PyRanges(df)
    except:
        return print("Failed")

def asPyRangesDecorator(func):
    """
    Decorator function to change pandas DataFrame colnames for PyRanges compatibility. Same as the above but in decorator form!
    """
    def wrapper(*args, **kwargs):
        df = func(*args, **kwargs)
        return asPyRanges(df)
    return wrapper

@asPyRangesDecorator
def Modkit2Pr(path, min_depth=10, max_depth=True, keep_raw=False):
    return readModkit(path, min_depth, max_depth, keep_raw)

@asPyRangesDecorator
def Modbam2Pr(path, min_depth=10, max_depth=False, keep_raw=False):
    return readModbam2bed(path, min_depth, max_depth, keep_raw)

@asPyRangesDecorator
def Bismark2Pr(path, mod, min_depth=10, max_depth=False, keep_raw=False):
    return readBismarkZeroCov(path, mod, min_depth, max_depth, keep_raw)

def loadChromSize():
    path = "./feature_references/mm39.chrom.sizes"

    df = pd.read_csv(path, sep="\t", names=["Chromosome", "End"])
    df["Start"] = 0

    return df[["Chromosome", "Start", "End"]]