import pandas as pd
import pyranges as pr
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from FeatureReferences import *
from scipy import stats

def geneRefPyRange():
    gene_ref_path = './feature_references/revised/GENCODE_Basic_mm39_Genes_merged.bed'
    df = pd.read_csv(gene_ref_path, sep="\t", names=["Chromosome", "Start", "End", "Name", "Strand"])
    return pr.PyRanges(df).unstrand()

def barplotRefPyRange():
    gene_feature_list = subprocess.check_output(["ls", "./feature_references/revised/gene_features/name_adjusted/"]).decode("utf-8").split("\n") 
    gene_feature_list.pop(-1)

    df_list = []
    for file in gene_feature_list:
        path = "./feature_references/revised/gene_features/name_adjusted/" + file
        feature_tsv = Features(path)
        feature_df = feature_tsv.toDF()
        df_list.append(feature_df)

    feature_reference_df = pd.concat(df_list).drop(columns=["Score"])
    return pr.PyRanges(feature_reference_df)

def CGIrefPyRange():
    cgi_feature_list = subprocess.check_output(["ls", "./feature_references/revised/cgi/named/"]).decode("utf-8").split("\n") 
    cgi_feature_list.pop(-1)

    cgi_df_list = []
    for file in cgi_feature_list:
        path = "./feature_references/revised/cgi/named/" + file
        cgi_tsv = CGIs(path)
        cgi_df = cgi_tsv.toDF()
        cgi_df_list.append(cgi_df)

    cgi_reference_df = pd.concat(cgi_df_list)
    return pr.PyRanges(cgi_reference_df)

class CpGIntersects(pr.PyRanges):
    """
    Main class for feature/gene level comparison. Inherits from PyRange. 
    """
    def __init__(self, df):
        super().__init__(df, df)

    def intersectGenes(self):
        """
        Intersects CpGs with genes. Based on gene start/end coordinates in the GENCODE Basic reference build. 
        """
        gene_ref_pr = geneRefPyRange()
        df_with_genes = self.join(gene_ref_pr, slack=0).as_df()
        df_with_genes["feature_type"] = "Gene"
        return  df_with_genes
    
    def intersectFeatures(self):
        """
        Intersects CpGs with genomic features. Output is a dataframe-type object. 
        """
        feature_ref = barplotRefPyRange().unstrand()
        df_with_features = self.join(feature_ref, slack=0).as_df()
        categories = ["Intergenic", "Promoter", "5UTR", "TSS", "Intron", "Exon", "3UTR", "TTS"]
        df_with_features["feature_type"] = pd.Categorical(df_with_features["feature_type"], categories)

        return  df_with_features
    
    def intersectCpGIslands(self):
        """
        Intersects CpGs with CpG islands. Islands are broken into island feature (i.e.: shelf, shore). 
        """
        cgi_ref = CGIrefPyRange().unstrand()
        df_with_cgis = self.join(cgi_ref, slack=0).as_df()
        categories = ["Open sea", "Upstream shelf", "Upstream shore", "CGI", "Downstream shore", "Downstream shelf"]
        df_with_cgis["feature_type"] = pd.Categorical(df_with_cgis["feature_type"], categories)

        return  df_with_cgis
    
    def groupByGenomicWindow(self, window_size):
        """
        Groups CpGs based according to 1kb windows ("tiles"), using the average (mean) hydroxymethlyation of CpGs within those windows. Output is distinct from the grouping function below as the chromosomal coordinates are actually what defines each cluster. 
        """
        tiled_pr = self.tile(window_size, strand=False).cluster(slack=-1, strand=False)
        cluster_pr = tiled_pr.merge(slack=-1, count=True, strand=False)
        cluster_pr = cluster_pr.insert(
            tiled_pr.apply(
            f=lambda df: df.groupby("Cluster")[["percentMeth_Nanopore_5hmC", "percentMeth_Bisulphite_5hmC"]].mean(), 
            as_pyranges=False, strand=False
            )
            )
        cluster_df = cluster_pr.as_df()
        return cluster_df.rename(columns={
            "percentMeth_Bisulphite_5hmC" : "percentMeth_TAB",
            "percentMeth_Nanopore_5hmC" : "percentMeth_Nanopore",
            "Count": "CpG_count"}
            )
    
    def group(self, intersect_with):
        """
        Groups CpGs based on intersects.
        """
        if intersect_with == "other":
            intersect_df = self.df
        elif intersect_with == "features" or intersect_with == "Features":
            intersect_df = self.intersectFeatures()
        elif intersect_with == "CGI" or intersect_with == "islands":
            intersect_df = self.intersectCpGIslands()
        elif intersect_with == "genes" or intersect_with == "Genes":
            intersect_df = self.intersectGenes()
        else: 
            raise ValueError("Please input appropriate element to intersect with. ['genes', 'features', or 'CGI']")
        
        groupby_df = intersect_df.groupby(["Name", "feature_type", "Start_b", "End_b"], 
                                          observed=True).agg(
            {"percentMeth_Nanopore_5hmC" : "mean",
             "percentMeth_Bisulphite_5hmC" : "mean",
             "Start" : "count"}
             ).reset_index()
        groupby_df.rename(columns={"Start" : "CpG_count",
                                   "Start_b" : "group_start",
                                   "End_b" : "group_end",
                                   "percentMeth_Bisulphite_5hmC" : "percentMeth_TAB",
                                   "percentMeth_Nanopore_5hmC" : "percentMeth_Nanopore"}, 
                                   inplace=True)

        return  groupby_df
       
    def calculateMethodMean(self):
        df = self.as_df()

        mean_dict = {"bisulphite_mean" : df["percentMeth_Bisulphite_5hmC"].mean(),
                     "nanopore_mean" : df["percentMeth_Nanopore_5hmC"].mean()}

        return mean_dict
    
class groupedDF:
    """
    Dataframe-type objects where CpG positions are grouped. Child classes contain additional functionality. Contains interface for relevant Seaborn plotting functions.  
    """
    def __init__(self, df, cpg_threshold=None):
        self.df = df
        self.cpg_threshold = cpg_threshold
      
    def dfWithLogCols(self):
        """
        Adds columns indicating the log2scale enrichment/depletion of grouped CpG sites relative to the mean of each method. 
        """
        filtered_df = self.df.copy()
        filtered_df = filtered_df.loc[filtered_df["CpG_count"].ge(self.cpg_threshold)] # filters out groups with fewer than 'cpg_threshold' CpGs present in either dataset. 

        with np.errstate(divide="ignore"):
            filtered_df["Log2FromMean_TAB"] = np.log2(
                np.divide(
                filtered_df["percentMeth_Bisulphite"],
                filtered_df["percentMeth_Bisulphite"].mean()
                )
                )

            filtered_df["Log2FromMean_Nanopore"] = np.log2(
                np.divide(
                filtered_df["percentMeth_Nanopore"],
                filtered_df["percentMeth_Nanopore"].mean()
                )
                )
        
        filtered_df = filtered_df.replace(-np.inf, np.nan).dropna()

        return filtered_df
    
    def methodComparison(self):
        """
        Adds "Average" and "Difference" to the dataframe, displaying the average level of enrichment and difference between method enrichment levels respectively.
        """
        df = self.dfWithLogCols()

        df["Average"] = df[["Log2FromMean_TAB", "Log2FromMean_Nanopore"]].mean(axis=1)
        df["Difference"] = np.subtract(df["Log2FromMean_TAB"], df["Log2FromMean_Nanopore"])
        
        return df
    
    def calcPearson(self):
        df = self.dfWithLogCols()

        return stats.pearsonr(df["Log2FromMean_TAB"],
                              df["Log2FromMean_Nanopore"])
    
    def calcSpearman(self):
        df = self.dfWithLogCols()

        return stats.spearmanr(df["Log2FromMean_TAB"],
                               df["Log2FromMean_Nanopore"])
    
    def calcMannWhitney(self, alternative="two-sided"):
        """
        Performs a Mann-Whitney(-Wilcoxon) U Test on the null hypothesis that the enrichment values from bisulphite (X) are from the same distribution as nanopore (Y). 
        Other values for the 'alternative' include 'greater' or 'less'. 
        """
        df = self.dfWithLogCols()

        return stats.mannwhitneyu(df["Log2FromMean_TAB"],
                                  df["Log2FromMean_Nanopore"], 
                                  alternative=alternative)
    
    def makeHist(self, stat, ax=None):
        df = self.dfWithLogCols()

        if not ax:
            fig, ax = plt.subplots()

        hist = sns.histplot(df, x="Log2FromMean_TAB", y="Log2FromMean_Nanopore", cbar=True, cbar_kws={"label" : f"{stat}".capitalize()}, stat=stat, ax=ax)

        return hist