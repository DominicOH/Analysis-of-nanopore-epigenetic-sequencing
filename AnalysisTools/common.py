import pandas as pd
from concurrent import futures
from AnalysisTools import read_modbed
import os 
import subprocess
import concurrent.futures

def read_table(path, usecols: list=None, test_run: bool=False):
    default_usecols = ["Chromosome", "Start", "End"]

    if test_run:
        nrows=100000
    else: 
        nrows=None

    if usecols:
        if type(usecols) == list:
            default_usecols.extend(usecols)
        else: 
            default_usecols.append(usecols)

    df = pd.read_table(path, sep="\t", 
                       usecols=default_usecols,
                       nrows=nrows)
    return df

def fetch_modbeds(dirpath, usecols=None, test_run: bool = False):
    path_ls = subprocess.check_output(["ls", dirpath]).decode("utf-8").split("\n")
    path_ls.pop(-1)

    with concurrent.futures.ThreadPoolExecutor(len(path_ls)) as read_executor:
        tables = read_executor.map(lambda path: read_table(dirpath + path, usecols, test_run), path_ls)
        tables = [table.rename(columns={"N_mC" : "N_5mC", "N_hmC" : "N_5hmC"}) for table in tables]
    
    return tables

def fetch_controls(usecols, test_run=False):
    zymo_controls = ["data/modbases/controls/zymo_wga_methylated_rep1.sorted.modbed", "data/modbases/controls/zymo_wga_methylated_rep2.sorted.modbed",
                     "data/modbases/controls/zymo_wga_unmodified_rep1.sorted.modbed", "data/modbases/controls/zymo_wga_unmodified_rep2.sorted.modbed"]

    with concurrent.futures.ProcessPoolExecutor(4) as ppe:
        futures = [ppe.submit(read_table, path, usecols, test_run) for path in zymo_controls]
        ctrl_dfs = [future.result() for future in futures]

        [df.rename(columns={"N_mC" : "N_5mC", "N_hmC" : "N_5hmC"}, inplace=True) for df in ctrl_dfs]

    return ctrl_dfs

def calculate_percentages(df: pd.DataFrame, cols) -> pd.DataFrame:
    if type(cols) == list:
        for col in cols:
            df.eval(f"percentMeth_{col.split('_')[1]} = ({col}/readCount)*100", inplace=True)
    else:
        df.eval(f"percentMeth_{cols.split('_')[1]} = ({cols}/readCount)*100", inplace=True)

    return df

def merge_positions(dfs, calculate_percentage=False, cols=None, drop=None) -> pd.DataFrame:
    merged = (pd.concat(dfs)
              .groupby(["Chromosome", "Start", "End"], observed=False, as_index=False, sort=False)
              .sum(numeric_only=True))

    # Replace with two functions/methods that do this
    if calculate_percentage:
        merged = calculate_percentages(merged, cols)

    if drop:
        merged = merged.drop(columns=drop)

    return merged
