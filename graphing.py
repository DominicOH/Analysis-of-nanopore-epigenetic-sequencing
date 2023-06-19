import seaborn as sns
import pandas as pd
import numpy as np
from sklearn import preprocessing
from sklearn import metrics

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

class ROCPlot:
    """
    Main object for building a ROC plot from the merged wideform dataframes. 
    """

    def __init__(self, dataframe, threshold=None):
        self.dataframe = dataframe
        self.threshold = threshold

    def ROCbinariser(self):
        """
        Binarises the bisulphite dataset (mCpG/readCount < 50% = unmodified, mCpG/readCount > 50% = modified).
        """
        if not self.threshold:
            binariser = preprocessing.Binarizer(threshold=50)
        else:
            binariser = preprocessing.Binarizer(threshold=self.threshold)

        self.dataframe["binarised"] = binariser.fit_transform(np.reshape(self.dataframe.iloc[:, 5].to_numpy(), (-1, 1)))
        return self.dataframe

    def ROC(self):
        """
        Retrieves the false positive and true positive rates of Nanopore modcalls relative to the binarised bisulphite data.
        """
        binarised_df = self.ROCbinariser()
        fpr, tpr, threshold = metrics.roc_curve(binarised_df["binarised"], binarised_df.iloc[:, 7])
        return fpr, tpr
    
    def plotROC(self, ax, label):
        """
        Uses the methods above to plot a ROC curve on a given ax. 
        """
        fpr, tpr = self.ROC()
        return ax.plot(fpr, tpr, "-", label=label, lw=2)
    
    def calculateAUC(self):
        """
        For a given ROC curve, calculates the area under that curve (AUC).
        """
        binarised_df = self.ROCbinariser()
        return metrics.roc_auc_score(binarised_df["binarised"], binarised_df.iloc[:, 7])