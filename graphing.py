import seaborn as sns
import pandas as pd
import numpy as np

def linePlot(df, ax, reverse=False):
    if df.columns.__contains__("percentMeth_subtraction_5hmC"):
        df = df.loc[df["percentMeth_subtraction_5hmC"] >= 0]
        if reverse==False:
            df["Bins"] = pd.cut(df["percentMeth_TAB_5hmC"], 101, labels=np.arange(0, 101, 1))
            return sns.lineplot(df, x="Bins", y="percentMeth_subtraction_5hmC", errorbar=("pi", 50), estimator="median", label="vs. TAB", ls=":", legend="brief", ax=ax)
        else:
            df["Bins"] = pd.cut(df["percentMeth_Nanopore_5hmC"], 101, labels=np.arange(0, 101, 1))
            return sns.lineplot(df, x="Bins", y="percentMeth_subtraction_5hmC", errorbar=("pi", 50), estimator="median", label="vs. Nanopore", ls="-", legend="brief", ax=ax)
            
    elif df.columns.__contains__("percentMeth_subtraction_5mC"):
        df = df.loc[df["percentMeth_subtraction_5mC"] >= 0]
        if reverse==False:
            df["Bins"] = pd.cut(df["percentMeth_oxBS_5mC"], 101, labels=np.arange(0, 101, 1))
            return sns.lineplot(df, x="Bins", y="percentMeth_subtraction_5mC", errorbar=("pi", 50), estimator="median", label="vs. oxBS", ls=":", legend="brief", ax=ax)
        else:
            df["Bins"] = pd.cut(df["percentMeth_Nanopore_5mC"], 101, labels=np.arange(0, 101, 1))
            return sns.lineplot(df, x="Bins", y="percentMeth_subtraction_5mC", errorbar=("pi", 50), estimator="median", label="vs. Nanopore", ls="-", legend="brief", ax=ax)
            
    elif df.columns.__contains__("percentMeth_WGBS"):
        df["Bins"] = pd.cut(df["percentMeth_WGBS"], 101, labels=np.arange(0, 101, 1))
        return sns.lineplot(df, x="Bins", y="percentMeth_Nanopore", errorbar=("pi", 50), estimator="median", label="vs. WGBS", ls="-", legend="brief", ax=ax)
        
    elif df.columns.__contains__("percentMeth_oxBS_5mC"):
        df["Bins"] = pd.cut(df["percentMeth_oxBS_5mC"], 101, labels=np.arange(0, 101, 1))
        return sns.lineplot(df, x="Bins", y="percentMeth_Nanopore_5mC", errorbar=("pi", 50), estimator="median", label="vs. oxBS", ls="--", legend="brief", ax=ax)
        
    elif df.columns.__contains__("percentMeth_TAB_5hmC"):
        if reverse==False:
            df["Bins"] = pd.cut(df["percentMeth_TAB_5hmC"], 101, labels=np.arange(0, 101, 1))
            return sns.lineplot(df, x="Bins", y="percentMeth_Nanopore_5hmC", errorbar=("pi", 50), estimator="median", label="vs. TAB", ls=":", legend="brief", ax=ax)
        else:
            df["Bins"] = pd.cut(df["percentMeth_Nanopore_5hmC"], 101, labels=np.arange(0, 101, 1))
            return sns.lineplot(df, x="Bins", y="percentMeth_TAB_5hmC", errorbar=("pi", 50), estimator="median", label="vs. Nanopore", ls="-", legend="brief", ax=ax)