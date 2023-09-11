import pandas as pd
import pyranges as pr
import FeatureReferences
import numpy as np
import warnings
from common import readBismarkZeroCov, readModbam2bed, asPyRangesDecorator

@asPyRangesDecorator
def Modbam2Pr(path):
    return readModbam2bed(path, 10, True)

@asPyRangesDecorator
def Bismark2Pr(path, mod):
    return readBismarkZeroCov(path, mod, 10, True)

def makeCpGRange(nanopore_path, tab_path):
    nanopore_pr = Modbam2Pr(nanopore_path)
    tab_pr = Bismark2Pr(tab_path, "5hmC")

    merged_df = nanopore_pr.join(tab_pr, False, suffix="_TAB").as_df()

    merged_df = merged_df.rename(columns={
        "percentMeth_5hmC" : "percentMeth_5hmC_Nanopore"
    })
    return CpGRange(merged_df)

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
            "percentMeth_5hmC_Prom" : "percentMeth_5hmC_Nanopore",
            "percentMeth_5hmC_Bisulphite" : "percentMeth_5hmC_TAB"}, 
            errors="ignore")
        
        super().__init__(df, df)
    
    @property
    def raw_means(self):
        """
        This is intended to be used for downstream applications - saving the mean modification % per CpG. 
        """
        list_of_cols = ["percentMeth_5hmC_TAB", "percentMeth_5hmC_Nanopore"]
        raw_means = {}
        df = self.df
        for col in list_of_cols:
            try:
                raw_means.update({
                    col : df[col].mean()
                })
            except: 
                pass
        return raw_means

    def __annotate_with_multiple(self, annotation_dir_path):
        annotation_ref = FeatureReferences.featureRefPyRange(annotation_dir_path).unstrand()
        with warnings.catch_warnings():
            warnings.simplefilter(action="ignore", category=FutureWarning)
            annotated_df = annotation_ref.join(self, False, "right", suffix="_CpG", apply_strand_suffix=False).as_df()

        return annotated_df
    
    def __annotate_with_single(self, feature_path):
        annotation_ref = FeatureReferences.Reference(feature_path).df
        # drops redundant columns 
        annotation_ref = annotation_ref.drop(
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
            {"percentMeth_5hmC_Nanopore" : np.mean,
             "percentMeth_5hmC_TAB" : np.mean,
             "Cluster" : "count"}
             ).reset_index().rename(columns={"Cluster":"CpG_count"})
           
        return Multisite(grouped_df, raw_means=self.raw_means)
    
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

        groupby_df = intersect_df.groupby(["feature_type", "Start", "End"]).agg({
            "percentMeth_5hmC_Nanopore" : np.mean,
            "percentMeth_5hmC_TAB" : np.mean,
            "Start_CpG" : "count"}).reset_index()
        
        if intersecting_on in ["genes", "features"]: 
            output_df = groupby_df.replace("-1", "Intergenic")
        elif intersecting_on == "cgi":
            output_df = groupby_df.replace("-1", "Open sea")
        elif intersecting_on == "repeats":
            output_df = groupby_df.replace("-1", None).dropna()

        return  Multisite(output_df.rename(columns={"Start_CpG" : "CpG_count"}), raw_means=self.raw_means)
    
class Multisite:
    """
    PyRange objects where CpG positions are grouped by feature, gene, or genomic window.

    Note: Not currently built to accommodate CpG 5mC.
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
            log2enrichment_5hmC_Nanopore=self.__calculate_log2_difference("percentMeth_5hmC_Nanopore", include_zeros),            
            log2enrichment_5hmC_TAB=self.__calculate_log2_difference("percentMeth_5hmC_TAB", include_zeros)
            )
        
        return Multisite(df, raw_means=self.raw_means)