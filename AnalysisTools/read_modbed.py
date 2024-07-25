#!/usr/bin/python3

"""
Functions to load and process the files produced using Oxford Nanopore Technologies Modkit and Bismark for downstream analysis.

"""

##### Module imports #####

import pandas as pd
import argparse
from math import sqrt
import numpy as np
import time
try: 
    from helpers import timer
except:
    from AnalysisTools.helpers import timer

##### Function definitions #####

def checkBisOrModkit(path):
    with open(path, "r") as file:
        fline = file.readline()
        line_len = len("".join(fline).split("\t"))
        if line_len == 18:
            return "Modkit"
        elif line_len == 6:
            return "Bismark"
        elif line_len == 7:
            return "Bismark CpG report"
        else:
            raise ValueError("Format not recognised.")

def filterDepth(df, 
                min_depth: int = 1, 
                apply_max_depth: bool = False):
    """
    Filters the dataframe to rows meeting the user-supplied minimum coverage depth. 
    Maximum depth filtering can also be applied to reduce influence of repetitive elements.
    """
    
    filtered_df = df.copy()
    average = filtered_df["readCount"].mean()
    filtered_df = filtered_df.loc[filtered_df.loc[:, "readCount"].ge(int(min_depth))]

    if apply_max_depth:
        filtered_df = filtered_df.loc[filtered_df.loc[:, "readCount"] < (average + 3*sqrt(average))]

    return filtered_df

def readModkit(
        path: str, 
        min_depth: int = 1,
        apply_max_depth: bool = False,
        include_raw_counts: bool = False
):
    """
    Reads the bedmMethyl output of Modkit pileup into a pd.DataFrame. 
    Note: It's important that modkit pileup is run with the --only-tabs flag. Otherwise important data columns are separated only by spaces and missed by tab-parsing. 

    :param str path: Filepath of the bedMethyl file. 
    :param int min_depth: The minimum readcount of CpG sites.
    :param bool apply_max_depth: Whether to filter out modbases with a depth greater than d + 3*sqrt(d); where d is the mean depth.
    :param bool incl_raw_counts: Whether the raw count of modified basecalls should be kept in the resulting dataframe.

    """
    colnames = ["Chromosome", "Start", "End", "modBase", "modScore", "Strand", "rem1", "rem2", "rem3", "readCount", "percentMeth", "N_mod", "N_canonical", "N_other", "N_delete", "N_fail", "N_diff", "N_nocall"]
    df_init = pd.read_csv(path, 
                          sep="\t", 
                          names=colnames,
                          dtype={
                              "readCount" : np.int16,
                              "N_mod" : np.int16, 
                              "N_canonical" : np.int16, 
                              "N_other" : np.int16, 
                              "N_delete" : np.int16, 
                              "N_fail" : np.int16, 
                              "N_diff" : np.int16, 
                              "N_nocall" : np.int16
                                 }
)
    
    df_filtered = filterDepth(df_init, min_depth, apply_max_depth)

    if include_raw_counts:
        pivot_cols =  ["readCount", "percentMeth", "N_mod", "N_canonical", "N_other"]
    else:
        pivot_cols = ["readCount", "percentMeth"]

    df_pivot = df_filtered.pivot(index=["Chromosome", "Start", "End", "Strand"], columns="modBase", values=pivot_cols)
    df_pivot.columns = df_pivot.columns.to_flat_index()
    df_pivot = df_pivot.reset_index()
    
    df_pivot = df_pivot.rename(columns={
        ('readCount', 'h') : "readCount",
        "h" : "percentMeth_5hmC",
        "m" : "percentMeth_5mC",
        ("percentMeth", "h") : "percentMeth_5hmC",
        ("percentMeth", "m") : "percentMeth_5mC",
        ("N_mod", "h") : "N_hmC",
        ("N_mod", "m") : "N_mC",
        ("N_canonical", "h") : "N_C"
    }, errors="ignore")

    df_pivot = df_pivot.drop(columns=[
        ("readCount", "m"),
        ("N_canonical", "m"),
        ("N_other", "h"),
        ("N_other", "m")
    ], errors="ignore")

    return df_pivot 

def readBismarkZeroCov(
        path: str, 
        modbase: str = None, 
        min_depth: int = 1,
        apply_max_depth: bool = False,
        include_raw_counts: bool = False):
    
    """
    Reads the output file of Bismark methylation extractor. Requires the bed format output produced with the options: -p --bedGraph --zero_based --comprehensive

    :param str path: Path to the bed file of modified bases.
    :param str mod: Type of target modification ["5mC", "5hmC"] 
    :param bool filter: Whether to filter the dataframe by depth [default=True]
    :return: Pandas DataFrame object 

    """
    if not modbase:
        meth_col = "percentMeth_mod"
        mod_count = "N_mod"
    elif modbase == "5mC":
        meth_col = "percentMeth_5mC"
        mod_count = "N_mC"
    elif modbase == "5hmC":
        meth_col = "percentMeth_5hmC"
        mod_count = "N_hmC"
    else:
        raise ValueError("Mod_types entered must equal: '5mC' or '5hmC'")

    df = (pd.read_csv(path, sep="\t", names=[
            "Chromosome", "Start", "End", meth_col, mod_count, "N_unmod"],
            dtype={
                "mod_count" : np.int16,
                "N_unmod" : np.int16
            })
            .assign(readCount = lambda row: row[mod_count] + row["N_unmod"]))
        
    if min_depth:
        df = filterDepth(df, min_depth, apply_max_depth)

    df["readCount"] = df["readCount"].astype(np.int16)

    if include_raw_counts:
        return df
    else: 
        return df.drop(columns=["N_mC", "N_hmC", "N_mod", "N_unmod"], errors="ignore")

def open_single_file(file, min_depth=1, modbase=None, quiet=False, include_raw_counts=False):
    s = time.perf_counter()
    if checkBisOrModkit(file) == "Bismark":        
        mod_df = readBismarkZeroCov(file, modbase, min_depth, include_raw_counts=include_raw_counts)
    elif checkBisOrModkit(file) == "Modkit": 
        mod_df = readModkit(file, min_depth, include_raw_counts=include_raw_counts)
    e = time.perf_counter()

    if not quiet:
        print(f"Read {file} in {e-s} seconds.")
    return mod_df

@timer
def read_modbed(path, outpath, min_depth=1, include_raw_counts=False, quiet=False):
    if type(path) == list:
        for path, outpath in zip(path, outpath):
            print(f"Reading from {path}")

            df = open_single_file(path, min_depth, include_raw_counts=include_raw_counts, quiet=quiet)
            print(f"Success. Saving to {outpath}")
            df.to_csv(f"{outpath}", sep="\t", header=True, index=False)
        return
    elif type(path) == str:
        print(f"Reading from {path}")
        df = open_single_file(path, min_depth, include_raw_counts=include_raw_counts, quiet=quiet)
        print(f"Success. Saving to {outpath}")
        return df.to_csv(f"{outpath}", sep="\t", header=True, index=False)
        
##### Main function #####

if __name__=="__main__":
    parser = argparse.ArgumentParser(
                        prog = "read_modbed",
                        description = "Reads a bedfile containing modified base information, such as that produced by ONT's modkit or Bismark.")
    parser.add_argument("filenames", action="store", nargs="+")
    parser.add_argument("-o", "--outpath", action="store", nargs="+", type=str, required=True) 
    parser.add_argument("--include_raw_counts", action="store_true", default=False)
    parser.add_argument("--min_depth", action="store", dest="min_depth", default=1, type=int) 
    parser.add_argument("-q", "--quiet", action="store_true", default=False) 

    args = parser.parse_args()
    read_modbed(args.filenames, args.outpath,  min_depth=args.min_depth, include_raw_counts=args.include_raw_counts, quiet=args.quiet)
    print("Done")