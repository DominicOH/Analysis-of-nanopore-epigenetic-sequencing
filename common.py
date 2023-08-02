import pandas as pd
from math import sqrt

def filterDepth(df, 
                min_depth: int = 10, 
                apply_max_depth: bool = False):
    """
    Filters the dataframe to only rows within the minimum and maximum coverage depth. 
    """
    average = df["readCount"].mean()
    df = df.loc[df.loc[:, "readCount"] >= min_depth]

    if apply_max_depth:
        df = df.loc[df.loc[:, "readCount"] < (average + 3*sqrt(average))]

    return 

def readBismarkZeroCov(
        path: str, 
        mod: str, 
        min_depth: int = 10,
        apply_max_depth: bool = False):
    """
    Reads the output file of Bismark methylation extractor. Requires the bed format output produced with the options: -p --bedGraph --zero_based --comprehensive

    :param str path: Path to the bed file of modified bases.
    :param str mod: Type of target modification ["5mC", "5hmC"] 
    :param bool filter: Whether to filter the dataframe by depth [default=True]
    :return: Pandas DataFrame object 

    """
    if mod == "5mC":
        df = pd.read_csv(path, sep="\t", names=[
            "chromosome", "chromStart", "chromEnd", "percentMeth_5mC", "N_mod", "N_unmod"]
            ).assign(readCount = lambda row: row.N_mod + row.N_unmod)
    elif mod == "5hmC":
        df = pd.read_csv(path, sep="\t", names=[
            "chromosome", "chromStart", "chromEnd", "percentMeth_5hmC", "N_mod", "N_unmod"]
            ).assign(readCount = lambda row: row.N_mod + row.N_unmod)
        
    if min_depth:
        df = filterDepth(df, min_depth, apply_max_depth)

    return df.drop(columns=["N_mod", "N_unmod"])

def get_nanopore_twoMod(path):
    """
    CURRENTLY NOT IN USE.
    """
    # Requires updates to Modbam2bed 0.9.1 -e 
    # Requires rebasecalling two-modwise in Dorado 
    names = ["chromosome", "chromStart", "chromEnd", "modification_type", "score", "strand", "thickStart", "thickEnd", "RGB", "readCount", "percentMeth", "N_unmod", "N_mod", "N_filtered"]
    redundant_cols = ["score", "thickStart", "thickEnd", "RGB"]
    nanopore_mod_df = pd.read_csv(path, sep="\t", index_col=False, names=names).drop(columns=redundant_cols)
    nanopore_mod_df["trueReadCount"] = nanopore_mod_df["readCount"].sub(nanopore_mod_df["N_filtered"])

    rename_dict = {
    "otherMod_percentMeth" : "percentMeth",
    "trueReadCount" : "readCount",
    "otherMod_method" : "method",
    "otherMod_type" : "modification_type"}
    
    nanopore_mod_df["method"] = "Nanopore 2-state"
    nanopore_mod_df["modification_type"] = "5mC"
    out_mod_df = nanopore_mod_df[["chromosome", "chromStart", "chromEnd", "strand", "modification_type", "trueReadCount", "percentMeth", "method"]].rename(columns=rename_dict)


    return out_mod_df.dropna()

def readModbam2bed(path: str,
                   min_depth: int = 10, 
                   apply_max_depth: bool = False
                   ):
    """
    Opens Modbam2bed bedMethyl files in an appropriate format for this analysis. 
    
    Note: Requires modbam2bed extended output with options: -e --cpg -m 5mC (other mod is assumed to be 5hmC)
    """
    modbed = pd.read_csv(path, sep="\t", 
                names=["chromosome", "chromStart", "chromEnd", "mod_type", "score", "strand", "i1", "i2", "i3", "readCount", "percentMeth_mC", "N_C", "N_mC", "N_filt", "N_NA", "N_hmC"])
    modbed["readCount_T"] = modbed.loc[:, ("N_C", "N_mC", "N_hmC")].sum(axis="columns")

    modbed.drop(columns="readCount", inplace=True)
    modbed.rename(columns={"readCount_T" : "readCount"}, inplace=True)

    if min_depth:
        modbed = filterDepth(modbed, min_depth, apply_max_depth)

    modbed["percentMeth_C"] = modbed.loc[:, "N_C"].divide(modbed.loc[:, "readCount"]).multiply(100)
    modbed["percentMeth_5hmC"] = modbed.loc[:, "N_hmC"].divide(modbed.loc[:, "readCount"]).multiply(100)
    modbed["percentMeth_5mC"] = modbed.loc[:, "N_mC"].divide(modbed.loc[:, "readCount"]).multiply(100)

    return modbed.loc[:, ("chromosome", "chromStart", "chromEnd", "strand", "readCount", "percentMeth_C", "percentMeth_5mC", "percentMeth_5hmC")]

def changeColNamesForPR(df):
    assert "chromosome" and "chromStart" in df.columns

    out_df = df.copy()
    out_df.rename(columns={
        "chromosome" : "Chromosome", 
        "chromStart" : "Start",
        "percentMeth_TAB_5hmC" : "percentMeth_Bisulphite_5hmC"
    }, inplace=True)

    if not "chromEnd" in out_df.columns:
        out_df["End"] = out_df["Start"].add(1)
    else: out_df.rename(columns={
        "chromEnd" : "End"
    }, inplace=True)
    
    return out_df
    
def loadChromSize():
    path = "./feature_references/revised/mm39.chrom.sizes"

    df = pd.read_csv(path, sep="\t", names=["chromosome", "chromEnd"])
    df["chromStart"] = 0

    return df[["chromosome", "chromStart", "chromEnd"]]

def optimisedResample(merged_df, left, right):
    """
    Calculates the most common shared read count between two methods then extracts only CpG sites with that read count.
    :param str left: The column name for readcounts in the left dataframe.
    :param str right: The column name for readcounts in the right dataframe. 
    """
    resampling = merged_df.groupby([left])[right].value_counts().sort_values(ascending=False).reset_index(name="Count")
    maxid = resampling.loc[resampling.loc[:, left] == resampling.loc[:, right], "Count"].idxmax()
    most_common = resampling.iloc[maxid, :][left]
    count = resampling.iloc[maxid, :]["Count"]

    print(f"Most common readcount is {most_common} with {count} CpGs.")

    return merged_df.loc[(merged_df.loc[:, left] == most_common) & (merged_df.loc[:, right] == most_common)]
