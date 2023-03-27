import pandas as pd

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
    redundant_cols = ["score", "strand", "thickStart", "thickEnd", "RGB"]
    nanopore_mod_df = pd.read_csv(path, sep="\t", index_col=False, names=names).drop(columns=redundant_cols)
    nanopore_mod_df["trueReadCount"] = nanopore_mod_df["readCount"].sub(nanopore_mod_df["N_filtered"])

    rename_dict = {
    "otherMod_percentMeth" : "percentMeth",
    "trueReadCount" : "readCount",
    "otherMod_method" : "method",
    "otherMod_type" : "modification_type"}
    
    nanopore_mod_df["method"] = "Nanopore 5mC"
    nanopore_mod_df["modification_type"] = "5mC"
    out_mod_df = nanopore_mod_df[["chromosome", "chromStart", "chromEnd", "modification_type", "trueReadCount", "percentMeth", "method"]].rename(columns=rename_dict)


    return out_mod_df.dropna()


def get_nanopore_threeMod(path):
    names = ["chromosome", "chromStart", "chromEnd", "modification_type", "score", "strand", "thickStart", "thickEnd", "RGB", "readCount", "percentMeth", "N_unmod", "N_mod", "N_filtered", "N_noCall", "N_other"]
    redundant_cols = ["score", "strand", "thickStart", "thickEnd", "RGB"]

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
    
    out_mod_df1 = nanopore_mod_df[["chromosome", "chromStart", "chromEnd", "modification_type", "trueReadCount", "percentMeth", "method"]].rename(columns=rename_dict)
    out_mod_df2 = nanopore_mod_df[["chromosome", "chromStart", "chromEnd", "otherMod_type", "trueReadCount", "otherMod_percentMeth", "otherMod_method"]].rename(columns=rename_dict)

    return out_mod_df1.dropna(), out_mod_df2.dropna()


def get_nanopore_threeMod_wStrand(path):
    names = ["chromosome", "chromStart", "chromEnd", "modification_type", "score", "strand", "thickStart", "thickEnd", "RGB", "readCount", "percentMeth", "N_unmod", "N_mod", "N_filtered", "N_noCall", "N_other"]
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
    names=["chromosome", "chromStart", "chromEnd", "strand", "readCount", "percentMeth"]
    wgbs_df = pd.read_csv(path, 
                        names=names, 
                        sep="\t")
    wgbs_df["method"] = "WGBS"
    wgbs_df["modification_type"] = "5mC"
    wgbs_df = wgbs_df[["chromosome", "chromStart", "chromEnd", "modification_type", "readCount", "percentMeth", "method"]] 
    return wgbs_df