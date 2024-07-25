import pandas as pd
import pyranges as pr
import warnings
import numpy as np

ctcf_chip = pr.PyRanges(pd.read_table("data/ctcf/ENCSR000CBN/ChIP2MACS2/MACS2/ENCSR000CBN_peaks.narrowPeak", 
                                      names=["Chromosome", "Start", "End", "Name", "Pileup", 
                                             "Strand", "FoldDifference", "pValue", "qValue", 
                                             "Peak"],
                                      usecols=["Chromosome", "Start", "End", "FoldDifference", 
                                               "pValue", "qValue"]), 
                                      int64=True)

ctcf_motif = pr.PyRanges(pd.read_table("data/ctcf/MA0139.1.tsv", 
                                       names=["Chromosome", "Start", "End", "Name", "Score", "Score2", "Strand"],
                                       usecols=["Chromosome", "Start", "End", "Strand"]),
                                       int64=True)

print(len(ctcf_motif), "CTCF motifs")

ctcf_summits = pr.PyRanges(pd.read_table("data/ctcf/ENCSR000CBN/ChIP2MACS2/MACS2/ENCSR000CBN_summits.bed", 
                                         names=["Chromosome", "Start", "End", "Name", "Score"]),
                           int64=True).merge()

print(len(ctcf_summits), "CTCF summits")


def read_duplex_modbed(path, test_run=True, replicate=None):
    if test_run:
        nrows = 1000000
    else:
        nrows=None
    
    pattern_df = (pd.read_table(path, sep="\t", 
                            names=["Chromosome", "Start", "End", "Pattern", "readCount", "D0", "D1", "D2", "D3", "D4", "percentPattern", "N_Pattern", "N_Canonical", "N_Other", "D5", "D6", "D7", "D8"],
                            usecols=["Chromosome", "Start", "End", "Pattern", "readCount", "N_Pattern"],
                            nrows=nrows)
                            .pivot_table(values="N_Pattern", columns="Pattern", index=["Chromosome", "Start", "End"], fill_value=0)
                            .reset_index())
    print(f"Found {len(pattern_df)} sites in {path}")

    if replicate:
        pattern_df["Replicate"] = replicate

    return pattern_df

def merge_pattern(df):
        df = (df.rename(columns={"-,-,C" : "CC",
                                 "m,m,C" : "MM",
                                 "h,h,C" : "HH",
                                 "-,m,C" : "CM",
                                 "m,-,C" : "MC",
                                 "h,-,C" : "HC",
                                 "-,h,C" : "CH",
                                 "m,h,C" : "MH",
                                 "h,m,C" : "HM"})
                .eval("readCount = CC + MM + HH + CM + MC + HC + CH + MH + HM"))
        return PatternFrame(df)

def read_merge(path, test_run=True, replicate=None):
    return merge_pattern(read_duplex_modbed(path, test_run, replicate))

def eval_basecall_proportions(df: pd.DataFrame):
    count_sum = df["Count"].sum()
    df.eval("Proportion = Count / @count_sum", inplace=True)
    
    return df

class PatternFrame:
    def __init__(self, pattern_df: pd.DataFrame):
        self.pattern_df = pattern_df
        self.motif_patterns = motif_join(pattern_df)
        self.chip_patterns = chip_join(pattern_df)

    def melt_patterns(self, min_depth):
        df = (self.pattern_df
              .query(f"readCount > {min_depth}"))
        
        df = df.melt(id_vars=["Chromosome", "Start", "End"], 
                     value_vars=["CC", "MM", "HH", "CM" , "MC" , "HC" , "CH" , "MH" , "HM"],
                     var_name="Pattern", value_name="Count").query("Count > 0")
        return MeltDF(df, min_depth)
            
    def merge_chip_patterns(self, min_depth):
        df = (self.chip_patterns
              .query(f"readCount > {min_depth}"))
        
        df = df.melt(id_vars=["Chromosome", "Start", "End", "Strand_ChIP"], 
                     value_vars=["CC", "MM", "HH", "CM" , "MC" , "HC" , "CH" , "MH" , "HM"],
                     var_name="Pattern", value_name="Count")

        return MeltDF(df, min_depth)
    
    def merge_motif_patterns(self, min_depth):
        df = (self.motif_patterns
              .query(f"readCount > {min_depth}"))
        
        df = df.melt(id_vars=["Chromosome", "Start", "End", "Strand_Motif"], 
                     value_vars=["CC", "MM", "HH", "CM" , "MC" , "HC" , "CH" , "MH" , "HM"],
                     var_name="Pattern", value_name="Count")

        return MeltDF(df, min_depth)
    
    def site_summit_distances(self, min_depth, filter_distance):
        return self.melt_patterns(min_depth).distance_to_summit(filter_distance)
           
class MeltDF(PatternFrame):
    def __init__(self, pattern_df: pd.DataFrame, min_depth):
        super().__init__(pattern_df)
        self.min_depth = min_depth

    def replace_modnames(self):
        df = self.pattern_df.replace(["CM", "MC", "CH", "HC", "MH", "HM", "CC", "MM", "HH"], 
                                     ["C:5mC", "C:5mC", "C:5hmC", "C:5hmC", "5mC:5hmC", "5mC:5hmC", "C:C", "5mC:5mC", "5hmC:5hmC"])
        return MeltDF(df, self.min_depth)
        
    def piechart_data(self):
        """
        Outputs a dataframe counting all constitutive-modification states. Hetero-modification and hemi-methylation states are grouped. 
        """
        df = self.replace_modnames().pattern_df
        pie_data = df.groupby("Pattern")["Count"].sum().reset_index()

        return pie_data
    
    def distance_to_summit(self, filter_distance=np.inf):
        df = self.replace_modnames().pattern_df

        distances = (pr.PyRanges(df, int64=True)
                     .k_nearest(ctcf_summits, ties="first", suffix="_Summit")
                     .as_df())
        
        distances = (distances.groupby(["Distance", "Pattern"], observed=True).agg({"Count" : sum})
                     .reset_index())
        
        distances["abs"] = abs(distances["Distance"])
        distances = distances.loc[distances["abs"] <= filter_distance]        
        
        return distances
    
class DistDF:
    def __init__(self, df: pd.DataFrame):
        self.df = df
            
    def explode_reads(self) -> pd.DataFrame:
        df = self.df

        df =  df.loc[df.index.repeat(df['Count']), ("Pattern", "Distance", "abs")]
        return df
                
# removal of overlapping binding sites to prevent duplication downstream 
# binding site with least qValue chosen
with warnings.catch_warnings():
    warnings.simplefilter("ignore", FutureWarning)
    # will find the intersect of sites that are present on CTCF motifs, 
    # that overlap with summits from confident peaks
    chip_merge = pr.PyRanges(ctcf_motif
                             .join(ctcf_summits, apply_strand_suffix=False, suffix="_Summit")
                             .join(ctcf_chip, apply_strand_suffix=False, suffix="_Peak")
                             .as_df()
                             .sort_values("qValue", ascending=True)
                             .groupby(["Chromosome", "Start", "End"])
                             .head(1)
                             .reset_index(drop=True),
                             int64=True)
    
print(len(chip_merge), "CTCF motifs that overlap summits and peaks")

def motif_join(df: pd.DataFrame) -> pd.DataFrame:
    """
    Intersects with CTCF motifs - not necessarily those within ChIP summit peaks. 
    Overlapping motifs are merged. 
    """
    # Loading JASPAR CTCF binding sites 
    intersect = pr.PyRanges(df, int64=True).join(ctcf_motif, strandedness=False, apply_strand_suffix=True, suffix="_Motif").as_df()

    assert intersect.loc[:, ("Chromosome", "Start", "End")].duplicated().all() == False # type: ignore
        
    return intersect
    
def chip_join(df: pd.DataFrame) -> pd.DataFrame:
    """
    Intersects with CTCF motifs present in ChIP summit peaks. 
    Note that this joins with CTCF summits to remove adjacent motifs that may not be bound.
    """
    # Loading JASPAR CTCF binding sites 
    intersect = pr.PyRanges(df, int64=True).join(chip_merge, strandedness=False, apply_strand_suffix=True, suffix="_ChIP").as_df()

    assert intersect.loc[:, ("Chromosome", "Start", "End")].duplicated().all() == False # type: ignore
        
    return intersect