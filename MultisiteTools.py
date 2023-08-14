import subprocess
import pandas as pd
import pyranges as pr
from FeatureReferences import Reference
import numpy as np
from typing import Literal

def geneRefPyRange():
    gene_ref_path = './feature_references/revised/GENCODE_Basic_mm39_Genes_merged.bed'
    df = pd.read_csv(gene_ref_path, sep="\t", names=["Chromosome", "Start", "End", "Name", "Strand"])
    return pr.PyRanges(df).unstrand()

def featureRefPyRange(dir_path: str):
    """
    Takes a directory of BED4 or BED6 files containing lists of features and feature coordinates to be used for annotation purposes. 
    """
    gene_feature_list = subprocess.check_output(["ls", dir_path]).decode("utf-8").split("\n") 
    gene_feature_list.pop(-1) # removes the current directory dot node 

    df_list = []
    for file in gene_feature_list:
        path = dir_path + file
        feature_tsv = Reference(path)
        feature_df = feature_tsv.as_dataframe()
        df_list.append(feature_df)

    feature_reference_df = pd.concat(df_list).drop(columns=["Score", "ThickStart", "ThickEnd"])
    return pr.PyRanges(feature_reference_df)

def repeatTypeRefPyRange():
    repeat_ref_path = "./feature_references/revised/repeats/UCSC_rmsk_mm39_Repeat_IDd.bed"
    df = pd.read_csv(repeat_ref_path, sep="\t", names=["Chromosome", "Start", "End", "feature_type", "Name"])
    return pr.PyRanges(df).unstrand()

class CpGRange(pr.PyRanges):
    """
    Initial class for feature/gene level comparison. Inherits from PyRanges. 
    """
    def __init__(self, df):
        # give minion and prom datasets the same name for downstream compatibility
        df = df.rename(columns={
            "percentMeth_5mC_Min" : "percentMeth_5mC_Nanopore", 
            "percentMeth_5hmC_Min" : "percentMeth_5hmC_Nanopore", 
            "percentMeth_5mC_Prom" : "percentMeth_5mC_Nanopore", 
            "percentMeth_5hmC_Prom" : "percentMeth_5hmC_Nanopore"}, 
            errors="ignore")
        super().__init__(df, df)

    def __annotate_with_multiple(self, annotation_dir_path):
        annotation_ref = featureRefPyRange(annotation_dir_path).unstrand()
        annotated_df = self.join(annotation_ref, slack=0).as_df()
        return annotated_df
    
    def __annotate_with_single(self, feature_path):
        annotation_ref = Reference(feature_path)
        annotation_pr = pr.PyRanges(annotation_ref)
        return self.join(annotation_pr, slack=0).as_df()
                 
    def group_by_tile(self, window_size):
        """
        Groups CpGs based according to `window_size` bp windows ("tiles"), using the average (mean) hydroxymethlyation of CpGs within those windows. Output is distinct from the grouping function below as the chromosomal coordinates are actually what defines each cluster. 
        """
        tiled_pr = self.tile(window_size, strand=False).cluster(slack=-1, strand=False)

        grouped_df = tiled_df.groupby(["Chromosome", "Start", "End"], observed=True).aggregate(
            {"percentMeth_5mC_Nanopore" : np.mean,
             "percentMeth_5mC_Bisulphite" : np.mean,
             "percentMeth_5hmC_Nanopore" : np.mean,
             "percentMeth_5hmC_Bisulphite" : np.mean,
             "Cluster" : "count"}
             ).reset_index()
           
        return grouped_df.rename(columns={"Cluster":"CpG_count"})
    
    def group_by(self, 
              intersect_with: Literal["features", "CGI", "genes", "repeats"],
              annotation_path: str
              ):
        """
        Groups CpGs based on intersects.
        """
        intersecting_on = str(intersect_with).lower()
        if intersecting_on == "genes":
            intersect_df = self.annotate_with_single(annotation_path)
        elif intersecting_on == "features" or intersecting_on == "cgi" or intersecting_on == "repeats":
            intersect_df = self.annotate_with_multiple(annotation_path)

        groupby_df = intersect_df.groupby(["Name", "feature_type", "Start_b", "End_b"], 
                                          observed=True).agg({
                                              "percentMeth_5mC_Nanopore" : np.mean,
                                              "percentMeth_5mC_Bisulphite" : np.mean,
                                              "percentMeth_5hmC_Nanopore" : np.mean,
                                              "percentMeth_5hmC_Bisulphite" : np.mean,
                                              "Start" : "count"}
                                              ).reset_index()

        return  Multisite(groupby_df.rename(columns={
            "Start" : "CpG_count", 
            "Start_b" : "group_start", 
            "End_b" : "group_end"}))
    
class Multisite:
    """
    Dataframe-type objects where CpG positions are grouped. Child classes contain additional functionality. Contains interface for relevant Seaborn plotting functions.  
    """
    def __init__(self, df, cpg_threshold=1):
        self.df = df.loc[df.loc[:, "CpG_count"].ge(cpg_threshold)]
    
    def __ratioToMean(self, column: str):
        new_col = self.df[column].divide(self.df[column].mean())
        return new_col
    
    def __log2RatioWrapper(self, column: str, include_zeros=False):
        if include_zeros:
            add = 1
        else: 
            add = 0
        
        with np.errstate(divide="ignore"):
            log2_col = np.log2(
                np.add(self.__ratioToMean(column), add))
        return log2_col

    def enrichmentComparison(self, include_zeros=False):
        """
        Provides additional columns with log_2 scale scores showing enrichment relative to the mean for 5mC and 5hmC. 

        :param bool include_zeros: Whether groups with an average CpG modification of zero are kept. Adds 1 to the ratio calculation to avoid zero division. 
        """
        df = self.df.copy()
        df = df.assign(
            log2enrichment_5mC_Nanopore=self.__log2RatioWrapper("percentMeth_5mC_Nanopore", include_zeros), 
            log2enrichment_5mC_Bisulphite=self.__log2RatioWrapper("percentMeth_5mC_Bisulphite", include_zeros), 
            log2enrichment_5hmC_Nanopore=self.__log2RatioWrapper("percentMeth_5hmC_Nanopore", include_zeros),            
            log2enrichment_5hmC_Bisulphite=self.__log2RatioWrapper("percentMeth_5hmC_Bisulphite", include_zeros)
            )
        
        return df
