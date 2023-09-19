import pandas as pd
import pyranges as pr
import FeatureReferences
import numpy as np
import warnings
from common import *

def makeCpGRange(nanopore_path, bis_path, mod="5hmC"):
    nanopore_pr = Modbam2Pr(nanopore_path)
    bis_path = Bismark2Pr(bis_path, mod)

    if mod == "5hmC":
        merged_df = nanopore_pr.join(bis_path, False, suffix="_TAB").as_df()
    else: 
        merged_df = nanopore_pr.join(bis_path, False, suffix="_oxBS").as_df()

    merged_df = merged_df.rename(columns={
        f"percentMeth_{mod}" : f"percentMeth_{mod}_Nanopore"
    })
    return CpGRange(merged_df)

def repeatTypeRefPyRange():
    repeat_ref_path = "./feature_references/revised/repeats/UCSC_rmsk_mm39_Repeat_IDd.bed"
    df = pd.read_csv(repeat_ref_path, sep="\t", names=["Chromosome", "Start", "End", "feature_type", "Name"])
    return pr.PyRanges(df).unstrand()

def annotationPivot(df):
    pivoted_df = pd.wide_to_long(df, 
        stubnames=["percentMeth_5mC", "log2enrichment_5mC", "percentMeth_5hmC", "log2enrichment_5hmC"], 
        i=["feature_type", "Chromosome", "Start", "End"], 
        j="method", sep="_", suffix="\D+")
    return pivoted_df.reset_index()

class CpGRange(pr.PyRanges):
    """
    Initial class for feature/gene level comparison. Inherits from PyRanges. 
    """
    def __init__(self, df):
        # give minion and prom datasets the same name for downstream compatibility
        df = df.rename(columns={
            "percentMeth_5hmC_Bisulphite" : "percentMeth_5hmC_TAB",
            "percentMeth_5mC_Bisulphite" : "percentMeth_5mC_oxBS"
            }, 
            errors="ignore")
        
        super().__init__(df, df)

    @property
    def __percent_cols(self):
        # Note: would be good to handle dataframes that have both
        df = self.df

        defaults = ["percentMeth_5mC_Nanopore", "percentMeth_5mC_oxBS", "percentMeth_5hmC_Nanopore", "percentMeth_5hmC_TAB"]
        cols_present = []

        for col in defaults: 
            if col in df.columns:
                cols_present.append(col)
        return cols_present
        
    @property
    def raw_means(self):
        """
        This is intended to be used for downstream applications - saving the mean modification % per CpG. 
        """
        list_of_cols = ["percentMeth_5hmC_TAB", "percentMeth_5hmC_Nanopore", "percentMeth_5mC_oxBS", "percentMeth_5mC_Nanopore"]
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
    
    def __groupby_function(self, df_with_groups):

        percent_cols = self.__percent_cols
        print(f"Aggregating all of {percent_cols}")

        if "feature_type" in df_with_groups.columns:
            groups = ["feature_type", "Chromosome", "Start", "End"]
        else: 
            groups = ["Chromosome", "Start", "End"]    

        if len(percent_cols) == 2:
            grouped_df = df_with_groups.groupby(groups, observed=True).aggregate(
                {percent_cols[0] : np.mean,
                percent_cols[1] : np.mean,
                "Cluster" : "count"}
                ).reset_index().rename(columns={"Cluster" : "CpG_count"})
        elif len(percent_cols) == 4:
            grouped_df = df_with_groups.groupby(groups, observed=True).aggregate(
                {percent_cols[0] : np.mean,
                percent_cols[1] : np.mean,
                percent_cols[2] : np.mean,
                percent_cols[3] : np.mean,
                "Cluster" : "count"}
                ).reset_index().rename(columns={"Cluster" : "CpG_count"})
        return grouped_df
                 
    def group_by_tile(self, window_size):
        """
        Groups CpGs based according to `window_size` bp windows ("tiles"), using the average (mean) hydroxymethlyation of CpGs within those windows. 
        
        :returns Multisite grouped_tiles:  
        """
        tiled_df = self.tile(window_size, strand=False).cluster(slack=-1, strand=False).as_df()

        grouped_df = self.__groupby_function(tiled_df)
        
        grouped_tiles = Multisite(grouped_df, raw_means=self.raw_means, percent_cols=self.__percent_cols)
           
        return grouped_tiles
    
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

        intersect_df = intersect_df.rename(columns={"Start_CpG" : "Cluster"})

        grouped_df = self.__groupby_function(intersect_df)
        
        if intersecting_on in ["genes", "features"]: 
            output_df = grouped_df.replace("-1", "Intergenic")
        elif intersecting_on == "cgi":
            output_df = grouped_df.replace("-1", "Open sea")
        elif intersecting_on == "repeats":
            output_df = grouped_df.replace("-1", None).dropna()

        return  Multisite(output_df.rename(columns={"Cluster" : "CpG_count"}), raw_means=self.raw_means, percent_cols=self.__percent_cols)
    
class Multisite:
    """
    PyRange objects where CpG positions are grouped by feature, gene, or genomic window.

    Note: Not currently built to accommodate CpG 5mC.
    """
    def __init__(self, df, cpg_threshold=1, raw_means=None, percent_cols=None):
        self._df = df.loc[df.loc[:, "CpG_count"].ge(cpg_threshold)]
        self._raw_means = raw_means
        self.__percent_cols = percent_cols

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

        percent_cols = self.__percent_cols
        new_cols = []

        for col in percent_cols:
            new_col_name = col.replace("percentMeth", "log2enrichment")
            new_col = self.__calculate_log2_difference(col, include_zeros)
            new_cols.append(pd.Series(new_col, name=new_col_name))

        for col in new_cols:
            df = df.assign(**{col.name : col.values})
        
        return Multisite(df, raw_means=self.raw_means)