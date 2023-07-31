import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pyranges as pr
from FeatureReferences import *
from CpGIntersects import CpGIntersects
from scipy import stats
   
class GroupedDF:
    """
    Dataframe-type objects where CpG positions are grouped. Child classes contain additional functionality. Contains interface for relevant Seaborn plotting functions.  
    """
    def __init__(self, df, cpg_threshold=None):
        self.df = df
        if "percentMeth_TAB_5hmC" in self.df.columns:
            self.df.rename(columns={"percentMeth_TAB_5hmC" : "percentMeth_TAB"}, inplace=True)
            
        self.cpg_threshold = cpg_threshold
      
    def dfWithLogCols(self, include_zeros=False):
        """
        Adds columns indicating the log2scale enrichment/depletion of grouped CpG sites relative to the mean of each method. 

        1 is added to the ratio of the group mean to the genomic mean to avoid ratios of 0 and resulting NaNs. 
        """
        filtered_df = self.df.copy()
        filtered_df = filtered_df.loc[filtered_df["CpG_count"].ge(self.cpg_threshold)] # filters out groups with fewer than 'cpg_threshold' CpGs present in either dataset. 

        if not include_zeros:
            add = 0
        else: 
            add = 1

        with np.errstate(divide="ignore"):
            filtered_df["Log2FromMean_TAB"] = np.log2(
                np.add(add, np.divide(filtered_df["percentMeth_TAB"], filtered_df["percentMeth_TAB"].mean()))
                )

            filtered_df["Log2FromMean_Nanopore"] = np.log2(
                np.add(add, np.divide(filtered_df["percentMeth_Nanopore"], filtered_df["percentMeth_Nanopore"].mean()))
                )
        
        if not include_zeros:
            filtered_df = filtered_df.replace(-np.inf, np.nan).dropna()

        return filtered_df
    
    def methodComparison(self):
        """
        Adds "Average" and "Difference" to the dataframe, displaying the average level of enrichment and difference between method enrichment levels respectively.
        """
        df = self.dfWithLogCols(False)

        df["Average"] = df[["Log2FromMean_TAB", "Log2FromMean_Nanopore"]].mean(axis=1)
        df["Difference"] = np.subtract(df["Log2FromMean_TAB"], df["Log2FromMean_Nanopore"])
        
        return df
    
    def calcPearson(self):
        df = self.dfWithLogCols(False)

        return stats.pearsonr(df["Log2FromMean_TAB"],
                              df["Log2FromMean_Nanopore"])
    
    def calcSpearman(self):
        df = self.dfWithLogCols(False)

        return stats.spearmanr(df["Log2FromMean_TAB"],
                               df["Log2FromMean_Nanopore"])
    
    def calcMannWhitney(self, alternative="two-sided"):
        """
        Performs a Mann-Whitney(-Wilcoxon) U Test on the null hypothesis that the enrichment values from bisulphite (X) are from the same distribution as nanopore (Y). 
        Other values for the 'alternative' include 'greater' or 'less'. 
        """
        df = self.dfWithLogCols(False)

        return stats.mannwhitneyu(df["Log2FromMean_TAB"],
                                  df["Log2FromMean_Nanopore"], 
                                  alternative=alternative)
    
    def makeHist(self, stat, ax=None, cax=None):
        df = self.dfWithLogCols(False)

        if not ax:
            fig, ax = plt.subplots()

        hist = sns.histplot(df, x="Log2FromMean_TAB", y="Log2FromMean_Nanopore", cbar=True, cbar_ax=cax, stat=stat, ax=ax)

        return hist
    
class FeatureAndGene(GroupedDF):
    """
    Dataframe-like objects where CpG sites are grouped by gene, genomic feature, or CpG island. 
    """
    def __init__(self, df, cpg_threshold=None):
        super().__init__(df)
        self.cpg_threshold = cpg_threshold

    def asLongDf(self):
        """
        Converts the DF from a wide-form to a longer form. 
        """
        cdf = self.df.copy()

        if "percentMeth_TAB_5hmC" in cdf.columns:
            cdf.rename(columns={"percentMeth_TAB_5hmC" : "percentMeth_TAB"}, inplace=True)

        stubs = ["percentMeth"]
        indices = ["Name", "group_start", "group_end"]
        
        return pd.wide_to_long(cdf, stubs, indices, "method", sep="_", suffix="\D+").reset_index()
    
    def makeLineplot(self, ax=None):
        df = self.asLongDf()

        if not ax:
            fig, ax = plt.subplots()
        
        lineplot = sns.lineplot(df, x="feature_type", y="percentMeth", hue="method", errorbar=("pi", 50), estimator="median", ax=ax)
        return lineplot
    
    def makeBarplot(self, ax=None):
        df = self.asLongDf()

        if not ax:
            fig, ax = plt.subplots()
        
        barplot = sns.barplot(df, x="feature_type", y="percentMeth", hue="method", errorbar=("pi", 50),  estimator="median", capsize=0.1, errwidth=1, ax=ax)
        return barplot

    def makeBoxplots(self, ax=None):
        df = self.asLongDf()

        if not ax:
            fig, ax = plt.subplots()
        
        boxplot = sns.boxplot(df, x="feature_type", y="percentMeth", hue="method", width=0.8, fliersize=0.005, ax=ax)
        return boxplot

class tiledGroup(GroupedDF):
    """
    Dataframe-like objects where CpG sites are grouped by genomic window or tile. 
    """
    def __init__(self, df, cpg_threshold):
        super().__init__(df, cpg_threshold)

    def positiveControlGroupDF(self, number_target_tiles):
        """
        Returns the positive control group - entries enriched in both methods.
        """
        df = super().methodComparison()
        df = df.nlargest(50, "Average")
        df = df.nsmallest(number_target_tiles, "Difference")
        return tiledGroup(df, self.cpg_threshold)
    
    def NegativeControlGroupDF(self, number_target_tiles):
        """
        Returns the negative control group - entries not enriched in either method.
        """
        df = self.methodComparison()
        df = df.nsmallest(50, "Average")
        df = df.nsmallest(number_target_tiles, "Difference")
        return tiledGroup(df, self.cpg_threshold)

    def NanoporePositiveGroupDF(self, number_target_tiles):
        """
        Returns the Nanopore positive test group - entries enriched only in Nanopore.
        """
        df = self.methodComparison()
        df = df.loc[df["Log2FromMean_TAB"] <= 0]
        df = df.nlargest(number_target_tiles, "Log2FromMean_Nanopore")
        return tiledGroup(df, self.cpg_threshold)
    
    def TabPositiveGroupDF(self, number_target_tiles):
        """
        Returns the TAB positive test group - entries enriched only in TAB.
        """
        df = self.methodComparison()
        df = df = df.loc[df["Log2FromMean_Nanopore"] <= 0]
        df = df.nlargest(number_target_tiles, "Log2FromMean_TAB")
        return tiledGroup(df, self.cpg_threshold)
    
    def getSequences(self):
        pyrange = pr.PyRanges(self.df)
        ref_fa = './data/Reference_data/mm39.fa'

        sequences = pr.get_sequence(pyrange, ref_fa)
        sequences.name = "sequence"

        pyrange = pyrange.insert(sequences)       
        return tiledGroup(pyrange.as_df(), self.cpg_threshold)
    
    def getIGVcoords(self):
        df = self.df
        coordinates = []
        for line, value in df.iterrows():
            chrom = value[0]
            s = str(value[1])
            e = str(value[2])
            first_half = ":".join([chrom, s])
            coord = "-".join([first_half, e])
            coordinates.append(coord)
        df["coordinates"] = coordinates
        return tiledGroup(df, self.cpg_threshold)
    
    def __reorderDF(self):
        df = self.df
        return df[["coordinates", "CpG_count", "percentMeth_Nanopore", "percentMeth_TAB", "Log2FromMean_TAB", "Log2FromMean_Nanopore", "Average", "Difference",  "sequence"]]

    def exportTests(self, number_target_tiles):
        wr = pd.ExcelWriter('/u/n/doh28/Documents/Nanopore_HMC/primer_design_regions.xlsx')

        list_of_groups = [self.NanoporePositiveGroupDF(number_target_tiles), self.TabPositiveGroupDF(number_target_tiles), self.positiveControlGroupDF(number_target_tiles), self.NegativeControlGroupDF(number_target_tiles)]
        processed_groups = [group.getSequences().getIGVcoords().__reorderDF() for group in list_of_groups]
        list_of_names = ["Nanopore_positive", "TAB_positive", "Positive_ctrl", "Negative_ctrl"]

        for i in np.arange(0, 4):
            processed_groups[i].to_excel(wr, list_of_names[i], index=False)

        wr.close()
        return 
    
    def tileWithLogCols(self):
        df = super().dfWithLogCols(False)
        return tiledGroup(df, self.cpg_threshold)
    
    def asCpGIntersect(self):
        df = self.df
        return CpGIntersects(df)
