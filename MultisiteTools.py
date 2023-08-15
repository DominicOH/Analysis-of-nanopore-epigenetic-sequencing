import subprocess
import pandas as pd
import pyranges as pr
from FeatureReferences import Reference
import numpy as np
from typing import Literal

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

    feature_reference_df = pd.concat(df_list).drop(columns=["Score", "ThickStart", "ThickEnd"], errors="ignore")
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
        self._raw_means = None

    @property
    def __raw_means(self):
        """
        This property is intended to be used for downstream applications - saving the mean modification % per CpG. 
        """
        list_of_cols = ["percentMeth_5mC_Bisulphite", "percentMeth_5hmC_Bisulphite", "percentMeth_5mC_Nanopore", "percentMeth_5hmC_Nanopore"]
        raw_means = {}
        for col in list_of_cols:
            raw_means.update({
                col : self.df[col].mean()
            })
        return raw_means

    def __annotate_with_multiple(self, annotation_dir_path):
        annotation_ref = featureRefPyRange(annotation_dir_path).unstrand()
        annotated_df = annotation_ref.join(self, False, "right", suffix="_CpG", apply_strand_suffix=False).as_df()

        return annotated_df
    
    def __annotate_with_single(self, feature_path):
        annotation_ref = Reference(feature_path).as_dataframe().drop(
            columns=["Score", "ThickStart", "ThickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"],
            errors="ignore")
        annotation_pr = pr.PyRanges(annotation_ref).unstrand()
        annotated_df = annotation_pr.join(self, False, "right", suffix="_CpG", apply_strand_suffix=False).as_df()
    
        return annotated_df
                 
    def group_by_tile(self, window_size):
        """
        Groups CpGs based according to `window_size` bp windows ("tiles"), using the average (mean) hydroxymethlyation of CpGs within those windows. Outputs a Pandas DataFrame. 
        """
        tiled_df = self.tile(window_size, strand=False).cluster(slack=-1, strand=False).as_df()

        grouped_df = tiled_df.groupby(["Chromosome", "Start", "End"], observed=True).aggregate(
            {"percentMeth_5mC_Nanopore" : np.mean,
             "percentMeth_5mC_Bisulphite" : np.mean,
             "percentMeth_5hmC_Nanopore" : np.mean,
             "percentMeth_5hmC_Bisulphite" : np.mean,
             "Cluster" : "count"}
             ).reset_index().rename(columns={"Cluster":"CpG_count"})
           
        return Multisite(grouped_df, raw_means=self.__raw_means)
    
    def group_by_annotation(self, 
              intersect_with: str,
              annotation_path: str
              ):
        """
        Groups CpGs based on intersecting annotations. Outputs a Pandas DataFrame.

        :param str intersect_with: Type of object to be annotated by (available: "genes", "features", "CGI", or "repeats")
        :param str annotation_path: Path to BED file or directory of BED files. BEDs must be in standard BED4, BED6, BED8, or BED12 format. 
        """
        intersecting_on = str(intersect_with).lower()
        if intersecting_on == "genes":
            intersect_df = self.__annotate_with_single(annotation_path).replace(-1, "Intergenic")
        elif intersecting_on in ["features", "cgi", "repeats"]:
            intersect_df = self.__annotate_with_multiple(annotation_path)
        else: 
            raise ValueError("Choose appropriate annotation type.")

        groupby_df = intersect_df.groupby(["Name", "feature_type", "Start", "End"]).agg({
            "percentMeth_5mC_Nanopore" : np.mean,
            "percentMeth_5mC_Bisulphite" : np.mean,
            "percentMeth_5hmC_Nanopore" : np.mean,
            "percentMeth_5hmC_Bisulphite" : np.mean,
            "Start_CpG" : "count"}).reset_index()
        
        if intersecting_on in ["genes", "features"]: 
            output_df = groupby_df.replace("-1", "Intergenic")
        elif intersecting_on == "cgi":
            output_df = groupby_df.replace("-1", "Open sea")
        elif intersecting_on == "repeats":
            output_df = groupby_df.replace("-1", None).dropna()

        return  Multisite(output_df.rename(columns={"Start_CpG" : "CpG_count"}), raw_means=self.__raw_means)
    
class Multisite:
    """
    PyRange objects where CpG positions are grouped by feature, gene, or genomic window.
    """
    def __init__(self, df, cpg_threshold=1, raw_means=None):
        self._df = df.loc[df.loc[:, "CpG_count"].ge(cpg_threshold)]
        self._raw_means = raw_means

    @property
    def df(self):
        return self._df
    
    @property
    def raw_means(self):
        return self._raw_means
    
    def __calculate_ratio_to_mean(self, column: str):
        """
        Calculates ratio difference between the grouped percentage modification and the ORIGINAL CPG MODIFICATION RATE.
        """
        new_col = self.df[column].divide(self.raw_means[column])
        return new_col
    
    def __calculate_log2_difference(self, column: str, include_zeros=False):
        if include_zeros:
            add = 1
        else: 
            add = 0
        
        with np.errstate(divide="ignore"):
            log2_col = np.log2(
                np.add(self.__calculate_ratio_to_mean(column), add))
        return log2_col

    def enrichment_over_mean(self, include_zeros=False):
        """
        Provides additional columns with log_2 scale scores showing enrichment relative to the mean for 5mC and 5hmC. 

        :param bool include_zeros: Whether groups with an average CpG modification of zero are kept. Adds 1 to the ratio calculation to avoid zero division. 
        """
        df = self.df.copy()
        df = df.assign(
            log2enrichment_5mC_Nanopore=self.__calculate_log2_difference("percentMeth_5mC_Nanopore", include_zeros), 
            log2enrichment_5mC_Bisulphite=self.__calculate_log2_difference("percentMeth_5mC_Bisulphite", include_zeros), 
            log2enrichment_5hmC_Nanopore=self.__calculate_log2_difference("percentMeth_5hmC_Nanopore", include_zeros),            
            log2enrichment_5hmC_Bisulphite=self.__calculate_log2_difference("percentMeth_5hmC_Bisulphite", include_zeros)
            )
        
        return Multisite(df, raw_means=self.raw_means)
    

