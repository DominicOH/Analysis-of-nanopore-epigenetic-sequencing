import pandas as pd
from math import sqrt

def get_bismark(path, mod):
    names = ["chromosome", "chromStart", "chromEnd", "percentMeth", "modified_reads", "unmodified_reads"]
    df = pd.read_csv(path, 
                         sep="\t", 
                         names=names)
    df["readCount"] = df[["modified_reads", "unmodified_reads"]].sum(axis=1)
    
    if mod == "5hmC":
        df["method"] = "TAB"
        df["modification_type"] = "5hmC"
    elif mod == "5mC":
        df["method"] = "oxBS"
        df["modification_type"] = "5mC"
    else: 
        raise ValueError("Please enter one of the compatible modification types: '5hmC' or '5mC'")
    
    df = df[["chromosome", "chromStart", "chromEnd", "modification_type", "readCount", "percentMeth", "method"]]

    return df

def get_nanopore_twoMod(path):
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

def get_nanopore_threeMod(path):
    names = ["chromosome", "chromStart", "chromEnd", "modification_type", "score", "strand", "thickStart", "thickEnd", "RGB", "readCount", "percentMeth_raw", "N_unmod", "N_mod", "N_filtered", "N_noCall", "N_other"]
    redundant_cols = ["score", "thickStart", "thickEnd", "RGB"]

    nanopore_mod_df = pd.read_csv(path, sep="\t", index_col=False, names=names).drop(columns=redundant_cols)
    nanopore_mod_df["trueReadCount"] = nanopore_mod_df["readCount"].sub(nanopore_mod_df["N_filtered"].add(nanopore_mod_df["N_noCall"]))
    nanopore_mod_df["otherMod_percentMeth"] = nanopore_mod_df["N_other"].divide(nanopore_mod_df["trueReadCount"], axis=0, fill_value=None).round(4).multiply(100)

    nanopore_mod_df["method"], nanopore_mod_df["otherMod_method"] = "Nanopore 5mC", "Nanopore 5hmC"
    nanopore_mod_df["modification_type"], nanopore_mod_df["otherMod_type"] = "5mC", "5hmC"

    rename_dict = {
    "otherMod_percentMeth" : "percentMeth",
    "trueReadCount" : "readCount",
    "otherMod_method" : "method",
    "otherMod_type" : "modification_type"}
    
    out_mod_df1 = nanopore_mod_df[["chromosome", "chromStart", "chromEnd", "strand", "modification_type", "trueReadCount", "percentMeth", "method"]].rename(columns=rename_dict)
    out_mod_df2 = nanopore_mod_df[["chromosome", "chromStart", "chromEnd", "strand", "otherMod_type", "trueReadCount", "otherMod_percentMeth", "otherMod_method"]].rename(columns=rename_dict)

    return out_mod_df1.dropna(), out_mod_df2.dropna()

def get_wgbs(path): 
    names=["chromosome", "chromStart", "chromEnd", "readCount", "percentMeth", "strand"]
    wgbs_df = pd.read_csv(path, 
                        names=names, 
                        sep="\t")
    wgbs_df["method"] = "WGBS"
    wgbs_df["modification_type"] = "5mC"
    wgbs_df = wgbs_df[["chromosome", "chromStart", "chromEnd", "strand", "modification_type", "readCount", "percentMeth", "method"]] 
    return wgbs_df

def filterDepth(df):
    average = df["readCount"].mean()
    df = df[df["readCount"].ge(10)]
    df = df[df["readCount"].le(average + 3*sqrt(average))]
    # df = df.loc[df["readCount"] == 15] # Testing constant readcount

    return df

def filterReadCountNanopore(df):
    mean = df.loc[:, "readCount_T"].mean() 
    max_depth = mean + 3*sqrt(mean) 
    filtered = df.loc[(df.loc[:, "readCount_T"] >= 10) & (df.loc[:, "readCount_T"] <= max_depth)]
    return filtered

def readModbam2bedTernary(path): 
    modbed = pd.read_csv(path, sep="\t", 
                names=["chromosome", "chromStart", "chromEnd", "mod_type", "score", "strand", "i1", "i2", "i3", "readCount", "percentMeth_mC", "N_C", "N_mC", "N_filt", "N_NA", "N_hmC"])
    modbed["readCount_T"] = modbed.loc[:, ("N_C", "N_mC", "N_hmC")].sum(axis="columns")
    modbed = filterReadCountNanopore(modbed)
    modbed["percentMeth_C"] = modbed.loc[:, "N_C"].divide(modbed.loc[:, "readCount_T"]).multiply(100)
    modbed["percentMeth_hmC"] = modbed.loc[:, "N_hmC"].divide(modbed.loc[:, "readCount_T"]).multiply(100)
    modbed["percentMeth_mC"] = modbed.loc[:, "N_mC"].divide(modbed.loc[:, "readCount_T"]).multiply(100)

    return modbed.loc[:, ("chromosome", "chromStart", "chromEnd", "strand", "readCount_T", "percentMeth_C", "percentMeth_mC", "percentMeth_hmC")]

def changeColNamesForPR(df):
    assert "chromosome" and "chromStart" in df.columns

    out_df = df.copy()
    out_df.rename(columns={
        "chromosome" : "Chromosome", 
        "chromStart" : "Start"
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
