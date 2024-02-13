#!/usr/bin/python3

"""
Functions to load and process the files produced using Oxford Nanopore Technologies Modkit and Bismark for downstream analysis.

Developed using mod_kit 0.1.13 and Bismark Version: v0.24.0. 
"""

##### Module imports #####

import pandas as pd
import argparse
import subprocess
from math import sqrt
import concurrent.futures

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
                min_depth: int = 10, 
                apply_max_depth: bool = False):
    """
    Filters the dataframe to rows meeting the user-supplied minimum coverage depth. 
    Maximum depth filtering can also be applied to reduce influence of repetitive elements.
    """
    
    filtered_df = df.copy()
    average = filtered_df["readCount"].mean()
    filtered_df = filtered_df.loc[filtered_df.loc[:, "readCount"] >= min_depth]

    if apply_max_depth:
        filtered_df = filtered_df.loc[filtered_df.loc[:, "readCount"] < (average + 3*sqrt(average))]

    return filtered_df

def readModkit(
        path: str, 
        min_depth: int = 1,
        apply_max_depth: bool = False,
        incl_raw_counts: bool = False
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
                          names=colnames)
    
    df_filtered = filterDepth(df_init, min_depth, apply_max_depth)

    if incl_raw_counts:
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
        modbase: str, 
        min_depth: int = 1,
        apply_max_depth: bool = False,
        incl_raw_counts: bool = False):
    
    """
    Reads the output file of Bismark methylation extractor. Requires the bed format output produced with the options: -p --bedGraph --zero_based --comprehensive

    :param str path: Path to the bed file of modified bases.
    :param str mod: Type of target modification ["5mC", "5hmC"] 
    :param bool filter: Whether to filter the dataframe by depth [default=True]
    :return: Pandas DataFrame object 

    """

    if modbase == "5mC":
        meth_col = "percentMeth_5mC"
        mod_count = "N_mC"
    elif modbase == "5hmC":
        meth_col = "percentMeth_5hmC"
        mod_count = "N_hmC"
    else:
        raise ValueError("Please enter a mod type: '5mC' or '5hmC'")

    df = pd.read_csv(path, sep="\t", names=[
            "Chromosome", "Start", "End", meth_col, mod_count, "N_unmod"]
            ).assign(readCount = lambda row: row[mod_count] + row["N_unmod"])
        
    if min_depth:
        df = filterDepth(df, min_depth, apply_max_depth)

    if incl_raw_counts:
        return df
    else: 
        return df.drop(columns=["N_mC", "N_hmC", "N_unmod"], errors="ignore")
    
def openReps(reps_path, select_cols=None, insert_cols=None, modbase=None, **kwargs):
    """
    Opens a directory of modkit '.bedMethyl' files or bismark_methylation_extractor '.zero.cov' files into a pd.Dataframe. 
    """
    if type(reps_path) == str:
        prep_iter = subprocess.check_output(["ls", reps_path]).split()
        rep_iter = []
        for file in prep_iter:
            rep_iter.append(reps_path + bytes.decode(file))
    else:
        rep_iter = reps_path

    rep_ls = []
    for index, file in enumerate(rep_iter): 
        if checkBisOrModkit(file) == "Bismark":
            if not modbase:
                raise ValueError("Bismark inputs requires modbase specification in args. ['5mC'|'5hmC']")
            
            mod_df = readBismarkZeroCov(file, modbase, **kwargs)
            mod_df = mod_df.assign(Replicate = f"Rep. {index + 1}")
            rep_ls.append(mod_df)

        elif checkBisOrModkit(file) == "Modkit": 
            mod_df = readModkit(file, **kwargs)
            mod_df = mod_df.assign(Replicate = f"Rep. {index + 1}")
            rep_ls.append(mod_df)
    
    df = pd.concat(rep_ls)

    if select_cols:
        df = df.loc[:, select_cols]

    if insert_cols:
        for col in insert_cols:
            df[col] = insert_cols[col]

    return df

def openReps_Parallel(reps_path, select_cols=None, insert_cols=None, modbase=None, **kwargs):
    """
    Opens a directory of modkit '.bedMethyl' files or bismark_methylation_extractor '.zero.cov' files into a pd.Dataframe. 
    """
    if type(reps_path) == str:
        prep_iter = subprocess.check_output(["ls", reps_path]).split()
        rep_iter = []
        for file in prep_iter:
            rep_iter.append(reps_path + bytes.decode(file))
    else:
        rep_iter = reps_path
        
    def open_single_rep(file, **kwargs):
        if checkBisOrModkit(file) == "Bismark":
            if not modbase:
                raise ValueError("Bismark input requires modbase specification in args. ['5mC'|'5hmC']")
            
            mod_df = readBismarkZeroCov(file, modbase, **kwargs)
        elif checkBisOrModkit(file) == "Modkit": 
            mod_df = readModkit(file, **kwargs)
        
        return mod_df

    with concurrent.futures.ThreadPoolExecutor() as tpe:
        futures = [tpe.submit(open_single_rep, file, **kwargs).result() for file in rep_iter]
        future_dfs = [future.assign(Replicate = f"Rep. {index + 1}") for index, future in enumerate(futures)]
    
    df = pd.concat(future_dfs)

    if select_cols:
        df = df.loc[:, select_cols]

    if insert_cols:
        for col in insert_cols:
            df[col] = insert_cols[col]

    return df

##### Main function #####

if __name__=="__main__":
    parser = argparse.ArgumentParser(
                        prog = "read_modbed",
                        description = "Reads a bedfile containing modified base information, such as that produced by ONT's modkit or Bismark.")

    parser.add_argument("-n ", "--nanopore_path", action="store", dest="nano_path", required=False) 
    parser.add_argument("-b ", "--BS_path", action="store", dest="BS_path", required=False) 
    parser.add_argument("-t ", "--TAB_path", action="store", dest="TAB_path", required=False) 
    parser.add_argument("-o ", "--oxBS_path", action="store", dest="oxBS_path", required=False) 
    parser.add_argument("-p ", "--out_prefix", action="store", dest="out_prefix", type=str, required=True) 
    parser.add_argument("--min_depth", action="store", dest="min_depth", default=10, type=int, required=False) 
    parser.add_argument("--max_depth", action="store", dest="max_depth", default=False, type=bool, required=False) 

    args = parser.parse_args()

    if not args.nano_path and not args.TAB_path and not args.oxBS_path:
        raise NameError("No valid file paths provided.")

    if args.nano_path:
        print("Nanopore data file found...")
        cpg_df = readModkit(args.nano_path, args.min_depth, args.max_depth).assign(method = "Nanopore")
        file_out_name = args.out_prefix + "_nanopore.tsv"
        cpg_df.to_csv(file_out_name, sep="\t", index=False)
        print(f"Processed file output as {file_out_name}")
   
    if args.TAB_path:
        print("TAB data file found...")
        cpg_df = readBismarkZeroCov(args.TAB_path, "5hmC", args.min_depth, args.max_depth).assign(method = "TAB-seq")
        file_out_name = args.out_prefix + "_tab.tsv"
        cpg_df.to_csv(file_out_name, sep="\t", index=False)
        print(f"Processed file output as {file_out_name}")
        
    if args.oxBS_path:
        print("oxBS data file found...")
        cpg_df = readBismarkZeroCov(args.oxBS_path, "5mC", args.min_depth, args.max_depth).assign(method = "oxBS-seq")
        file_out_name = args.out_prefix + "_oxbs.tsv"
        cpg_df.to_csv(file_out_name, sep="\t", index=False)
        print(f"Processed file output as {file_out_name}")

