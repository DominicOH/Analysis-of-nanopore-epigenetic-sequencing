import pandas as pd

def get_tab(path):
    names = ["chromosome", "chromStart", "chromEnd", "percentMeth", "modified_reads", "unmodified_reads"]
    tab_df = pd.read_csv(path, 
                         sep="\t", 
                         names=names)
    tab_df["readCount"] = tab_df[["modified_reads", "unmodified_reads"]].sum(axis=1)
    tab_df["method"] = "TAB"
    tab_df["modification_type"] = "5hmC"

    tab_df = tab_df[["chromosome", "chromStart", "chromEnd", "modification_type", "readCount", "percentMeth", "method"]]

    return tab_df

def get_nanopore_mods(path):
    redundant_cols = ["score", "strand", "thickStart", "thickEnd", "RGB"]
    nanopore_mod_df = pd.read_csv(path, sep="\t", index_col=False, names=names).drop(columns=redundant_cols)
    nanopore_mod_df["trueReadCount"] = nanopore_mod_df["readCount"].sub(nanopore_mod_df["N_filtered"].add(nanopore_mod_df["N_noCall"]))
    nanopore_mod_df["otherMod_percentMeth"] = nanopore_mod_df["N_other"].divide(nanopore_mod_df["trueReadCount"], axis=0, fill_value=None).round(4).multiply(100)

    rename_dict = {"otherMod_percentMeth" : "percentMeth",
                "trueReadCount" : "readCount"}
    out_mod_df1 = nanopore_mod_df[["chromosome", "chromStart", "chromEnd", "trueReadCount", "percentMeth"]].rename(columns=rename_dict)
    out_mod_df2 = nanopore_mod_df[["chromosome", "chromStart", "chromEnd", "trueReadCount", "otherMod_percentMeth"]].rename(columns=rename_dict)

    return out_mod_df1, out_mod_df2

def get_wgbs(path): 
    names=["chromosome", "chromStart", "chromEnd", "strand", "readCount", "percentMeth"]
    wgbs_df = pd.read_csv(path, 
                        names=names, 
                        sep="\t")
    wgbs_df["method"] = "WGBS"
    wgbs_df["modification_type"] = "5mC"
    wgbs_df = wgbs_df[["chromosome", "chromStart", "chromEnd", "modification_type", "readCount", "percentMeth", "method"]] 
    return wgbs_df