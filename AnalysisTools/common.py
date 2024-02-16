import pandas as pd
from concurrent import futures
import read_modbed
import time
from functools import wraps
import os 
import subprocess
import concurrent.futures

def onlyAutosomal(df):
    df = df.loc[df.loc[:, "Chromosome"].str.match("^(chr)\d+$")]
    return df

def fetch_data(dry_run: bool, split_biorep=False, **kwargs):
    """
    Fetches modkit/bismark data. Order of return is: cbm2, cbm3, tab, oxbs. 

    ::bool dry_run:: Whether to return data for internal testing.  
    """

    if dry_run:
        cbm2_path = "data/dryruns/cbm2/"
        cbm3_path = "data/dryruns/cbm3/" 

        oxbs_path = "data/dryruns/oxbs/"
        tab_path = "data/dryruns/tab/"

    else:
        cbm2_path = "data/modbases/nanopore/cbm2/"
        cbm3_path = "data/modbases/nanopore/cbm3/"

        oxbs_path = "data/modbases/public/CRR008808_oxBS/masked/"
        tab_path = "data/modbases/public/CRR008807_TAB/masked/"

    if split_biorep:
        bio_reps = ["Nanopore 1", "Nanopore 2"]
    else:
        bio_reps = ["Nanopore", "Nanopore"]

    with futures.ThreadPoolExecutor(4) as ppe:
        all_futures = [ppe.submit(read_modbed.openReps, path, insert_cols={"Technique" : bio_rep}, **kwargs) for path, bio_rep in zip([cbm2_path, cbm3_path], bio_reps)]
        all_futures.append(ppe.submit(read_modbed.openReps, tab_path, insert_cols={"Technique" : "TAB"}, modbase="5hmC", **kwargs))
        all_futures.append(ppe.submit(read_modbed.openReps, oxbs_path, insert_cols={"Technique" : "oxBS"},  modbase="5mC", **kwargs))
        future_dfs = [onlyAutosomal(future.result()) for future in all_futures]

    return future_dfs

def load_controls(dry_run: bool, **kwargs):
    """
    Fetches modkit/bismark data. Order of return is: zymo_modified, zymo_unmodified. 

    ::bool dry_run:: Whether to return data for internal testing.  
    """

    if not dry_run:
        zymo_modified = "data/modbases/merged_reps/zymo_methylated.merged.bed"
        zymo_unmodified = "data/modbases/merged_reps/zymo_unmodified.merged.bed"

    else:
        zymo_modified = "data/dryruns/merged/zymo_methylated.merged.bed.head"
        zymo_unmodified = "data/dryruns/merged/zymo_unmodified.merged.bed.head"

    paths = [zymo_modified, zymo_unmodified]

    with futures.ProcessPoolExecutor(2) as ppe:
        all_futures = [ppe.submit(pd.read_table, path, sep="\t", 
                                  names=["Chromosome", "Start", "End", "readCount", "N_mC", "N_hmC", "percentMeth_5mC", "percentMeth_5hmC"],
                                  **kwargs) for path in paths]
        all_dfs = [future.result() for future in all_futures]

    return all_dfs
    
def openReps(reps_path, select_cols=None, insert_cols=None, modbase=None, min_depth=1, quiet=True, **kwargs):
    """
    Opens a directory of modkit '.bedMethyl' files or bismark_methylation_extractor '.zero.cov' files into a pd.Dataframe. 
    """
    if not quiet:
        print(f"Reading from {reps_path}")

    if type(reps_path) == list:
        rep_iter = [filepath for filepath in reps_path]
    elif os.path.isdir(reps_path):
        dir_ls = subprocess.check_output(["ls", reps_path]).split()
        rep_iter = [reps_path + bytes.decode(filepath) for filepath in dir_ls]
    else:
        raise ValueError("Supplied path must be a dir or list of files.")

    with concurrent.futures.ThreadPoolExecutor() as tpe:
        futures = [tpe.submit(read_modbed.open_single_file, file, min_depth, modbase, **kwargs).result() for file in rep_iter]
        future_dfs = [future.assign(Replicate = f"Rep. {index + 1}") for index, future in enumerate(futures)]
    
    df = pd.concat(future_dfs)

    if insert_cols:
        for col in insert_cols:
            df[col] = insert_cols[col]

    if select_cols:
        df = df.loc[:, select_cols]

    return df


