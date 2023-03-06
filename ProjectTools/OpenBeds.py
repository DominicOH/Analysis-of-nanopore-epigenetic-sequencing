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

  def get_nanopore_5mc(path):
    names = ["chromosome", "chromStart", "chromEnd", "modification_type", "score", "strand", "thickStart", "thickEnd", "RGB", "readCount", "percentMeth", "unmodifiedReads", "modifiedReads", "filteredReads"]
    nanopore_df = pd.read_csv(path, sep="\t", names=names, index_col=False)
    nanopore_df.drop(columns=["score", "thickStart", "thickEnd", "RGB", "strand", "unmodifiedReads", "modifiedReads", "filteredReads"], inplace=True)
    nanopore_df["method"] = "Nanopore 5mC"
    nanopore_df.dropna(inplace=True)
   
    return nanopore_df

def get_nanopore_5hmc(path):
    names = ["chromosome", "chromStart", "chromEnd", "modification_type", "score", "strand", "thickStart", "thickEnd", "RGB", "readCount", "percentMeth", "unmodified_reads", "modified_reads", "filtered_reads", "nocall_reads"]
    nanopore_hmc_df = pd.read_csv(path, sep="\t", names=names, index_col=False)

    nanopore_hmc_df.drop(columns=["score", "strand", "thickStart", "thickEnd", "RGB", "unmodified_reads", "modified_reads", "filtered_reads", "nocall_reads"], inplace=True)
    nanopore_hmc_df["method"] = "Nanopore 5hmC"

    nanopore_hmc_df.dropna(inplace=True)
    
    return nanopore_hmc_df

def get_wgbs(path): 
    names=["chromosome", "chromStart", "chromEnd", "strand", "readCount", "percentMeth"]
    wgbs_df = pd.read_csv(path, 
                        names=names, 
                        sep="\t")
    wgbs_df["method"] = "WGBS"
    wgbs_df["modification_type"] = "5mC"
    wgbs_df = wgbs_df[["chromosome", "chromStart", "chromEnd", "modification_type", "readCount", "percentMeth", "method"]] 
    return wgbs_df