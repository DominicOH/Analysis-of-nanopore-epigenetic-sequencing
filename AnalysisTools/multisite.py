import pandas as pd
import pyranges as pr
from AnalysisTools import annotation_features
from AnalysisTools.read_modbed import readModkit, readBismarkZeroCov
import numpy as np
import warnings
from AnalysisTools.common import *

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
def Bismark2Pr(path, mod, min_depth=10, max_depth=False, keep_raw=False):
    return readBismarkZeroCov(path, mod, min_depth, max_depth, keep_raw)

def loadChromSize():
    path = "./feature_references/mm39.chrom.sizes"

    df = pd.read_csv(path, sep="\t", names=["Chromosome", "End"])
    df["Start"] = 0

    return df[["Chromosome", "Start", "End"]]

def makeCpGRange(nanopore_path, bis_path, mod="5hmC"):
    nanopore_pr = Modkit2Pr(nanopore_path)
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
        try: 
            df = df.rename(columns={
                "percentMeth_5hmC_Bisulphite" : "percentMeth_5hmC_TAB",
                "percentMeth_5mC_Bisulphite" : "percentMeth_5mC_oxBS"
                }, 
                errors="ignore")
        except:
            Warning("This is already a PyRange.")
        
        super().__init__(df, df)

    @property
    def __percent_cols(self):
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
        annotation_ref = annotation_features.featureRefPyRange(annotation_dir_path).unstrand()

        with warnings.catch_warnings():
             warnings.simplefilter(action="ignore", category=FutureWarning)
             annotated_df = annotation_ref.join(self, False, "right", suffix="_CpG", apply_strand_suffix=False).as_df()

        return annotated_df
    
    def __annotate_with_single(self, feature_path):
        annotation_ref = annotation_features.Reference(feature_path).df

        # drops redundant columns 
        annotation_ref = annotation_ref.drop(
            columns=["Score", "ThickStart", "ThickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"],
            errors="ignore")
        annotation_pr = pr.PyRanges(annotation_ref).unstrand()

        annotated_df = annotation_pr.join(self, False, "right", suffix="_CpG", apply_strand_suffix=False).as_df()
    
        return annotated_df
    
    def group_by_tile(self, window_size, agg_cols_funcs):
        """
        Groups CpGs based according to `window_size` bp windows ("tiles"), using the average (mean) hydroxymethlyation of CpGs within those windows. 
        
        :returns Multisite grouped_tiles:  
        """
        tiled_df = (self
                    .tile(window_size, strand=False)
                    .cluster(slack=-1, strand=False, count=True)
                    .as_df()
                    )
        grouped_df = (tiled_df
                      .groupby(["Chromosome", "Start", "End", "Cluster"], observed=True)
                      .agg(agg_cols_funcs)
                      .reset_index()
                      )
        
        grouped_tiles = Multisite(grouped_df, raw_means=self.raw_means, percent_cols=self.__percent_cols)
           
        return grouped_tiles
    
    def group_by_annotation(self, 
              intersect_with: str,
              annotation_path: str,
              agg_cols_funcs: dict,
              replace_gaps: dict
              ):
        """
        Groups CpGs based on intersecting annotations. Outputs a Pandas DataFrame.

        :param str intersect_with: Type of object to be annotated by (available: "genes", "features", "CGI", or "repeats")
        :param str annotation_path: Path to BED file or directory of BED files. BEDs must be in standard BED4, BED6, BED8, or BED12 format. 
        """
        intersecting_on = str(intersect_with).lower()
        if intersecting_on == "genes":
            intersect_df = self.__annotate_with_single(annotation_path)
        elif intersecting_on in ["features", "cgi", "repeats"]:
            intersect_df = self.__annotate_with_multiple(annotation_path)
        else: 
            raise ValueError("Choose appropriate annotation type.")

        grouped_df = (intersect_df
                      .groupby(["Chromosome", "Start", "End", "feature_type"], observed=True)
                      .agg(agg_cols_funcs)
                      .replace(replace_gaps)
                      )
        
        grouped_df = grouped_df

        return  Multisite(grouped_df, 
                          raw_means=self.raw_means, 
                          percent_cols=self.__percent_cols)
    
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
    
    def __calculate_ratio_to_mean(self, column: str, native=False):
        """
        Calculates ratio difference between the grouped percentage modification and the ORIGINAL CPG MODIFICATION RATE.
        """
        if not native:
            epsilon = 1
        else: 
            epsilon = 0

        x = self.df[column].add(epsilon)
        x_bar = self.raw_means[column] + epsilon
        
        ratio = x.divide(x_bar)
        return ratio
    
    def __calculate_log2_difference(self, column: str, native=False):
        with np.errstate(divide="ignore"):
            log2_col = np.log2(
                self.__calculate_ratio_to_mean(column, native))
        return log2_col

    def enrichment_over_mean(self, native=False):
        """
        Provides additional columns with log_2 scale scores showing enrichment relative to the mean for 5mC and 5hmC. 

        :param bool native: Whether the log transformation of the ratio is done as is. Features with average modification of 0% are lost. 
        """
        df = self.df.copy()

        percent_cols = self.__percent_cols
        new_cols = []

        for col in percent_cols:
            new_col_name = col.replace("percentMeth", "log2enrichment")
            new_col = self.__calculate_log2_difference(col, native)
            new_cols.append(pd.Series(new_col, name=new_col_name))

        for col in new_cols:
            df = df.assign(**{col.name : col.values})
        
        return Multisite(df, raw_means=self.raw_means)