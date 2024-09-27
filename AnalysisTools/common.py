import pandas as pd
from concurrent import futures
from AnalysisTools import read_modbed
import os 
import subprocess
import concurrent.futures

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

    with concurrent.futures.ThreadPoolExecutor(4) as tpe:
        futures = [tpe.submit(read_modbed.open_single_file, file, min_depth, modbase, quiet=quiet, **kwargs).result() for file in rep_iter]
        future_dfs = [future.assign(Replicate = f"Rep. {index + 1}") for index, future in enumerate(futures)]
    
    df = pd.concat(future_dfs)

    if insert_cols:
        for col in insert_cols:
            df[col] = insert_cols[col]

    if select_cols:
        df = df.loc[:, select_cols]

    return df

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

def fetch_modbeds(dirpath, usecols=None):
    path_ls = subprocess.check_output(["ls", dirpath]).decode("utf-8").split("\n")
    path_ls.pop(-1)

    with concurrent.futures.ThreadPoolExecutor(len(path_ls)) as read_executor:
        tables = read_executor.map(lambda path: read_table(dirpath + path, usecols), path_ls)
        tables = [table.rename(columns={"N_mC" : "N_5mC", "N_hmC" : "N_5hmC"}) for table in tables]
    
    return tables
    
def fetch_oxbs(usecols, dryrun=True):    
    if dryrun:
        root_path = "data/dryruns/modbeds/"

        oxbs_1_path = root_path + "CRD018546.gz_val_1_bismark_bt2_pe.deduplicated_mapq10_sorted.bismark.cov.gz.cpg_only.bed.CpG_report.txt.bed.masked.modbed.head"
        oxbs_2_path = root_path + "CRD018548.gz_val_1_bismark_bt2_pe.deduplicated_mapq10_sorted.bismark.cov.gz.cpg_only.bed.CpG_report.txt.bed.masked.modbed.head"

    else:
        root_path = "data/modbases/modbeds/"

        oxbs_1_path = root_path + "CRD018546.gz_val_1_bismark_bt2_pe.deduplicated_mapq10_sorted.bismark.cov.gz.cpg_only.bed.CpG_report.txt.bed.masked.modbed"
        oxbs_2_path = root_path + "CRD018548.gz_val_1_bismark_bt2_pe.deduplicated_mapq10_sorted.bismark.cov.gz.cpg_only.bed.CpG_report.txt.bed.masked.modbed"

    return map(lambda path: read_table(path, usecols).rename(columns={"N_mod" : "N_5mC"}), [oxbs_1_path, oxbs_2_path])

def fetch_tab(usecols, dryrun=True):
    if dryrun:
        root_path = "data/dryruns/modbeds/"

        tab_1_path = root_path + "CRD018526-8.merged.sorted.q10.bismark.cov.gz.cpg_only.bed.CpG_report.txt.bed.masked.modbed.head"
        tab_2_path = root_path + "CRD018542.gz_val_1_bismark_bt2_pe.q10.bismark.cov.gz.cpg_only.bed.CpG_report.txt.bed.masked.modbed.head"
        tab_3_path = root_path + "CRD018544.gz_val_1_bismark_bt2_pe.q10.bismark.cov.gz.cpg_only.bed.CpG_report.txt.bed.masked.modbed.head"

    else:
        root_path = "data/modbases/modbeds/"

        tab_1_path = root_path + "CRD018526-8.merged.sorted.q10.bismark.cov.gz.cpg_only.bed.CpG_report.txt.bed.masked.modbed"
        tab_2_path = root_path + "CRD018542.gz_val_1_bismark_bt2_pe.q10.bismark.cov.gz.cpg_only.bed.CpG_report.txt.bed.masked.modbed"
        tab_3_path = root_path + "CRD018544.gz_val_1_bismark_bt2_pe.q10.bismark.cov.gz.cpg_only.bed.CpG_report.txt.bed.masked.modbed"

    return map(lambda path: read_table(path, usecols).rename(columns={"N_mod" : "N_5hmC"}), [tab_1_path, tab_2_path, tab_3_path])

def fetch_controls(usecols, test_run=False):
    zymo_controls = ["data/modbases/controls/zymo_wga_methylated_rep1.sorted.modbed", "data/modbases/controls/zymo_wga_methylated_rep2.sorted.modbed",
                     "data/modbases/controls/zymo_wga_unmodified_rep1.sorted.modbed", "data/modbases/controls/zymo_wga_unmodified_rep2.sorted.modbed"]

    with concurrent.futures.ProcessPoolExecutor(4) as ppe:
        futures = [ppe.submit(read_table, path, usecols, test_run) for path in zymo_controls]
        ctrl_dfs = [future.result() for future in futures]

        [df.rename(columns={"N_mC" : "N_5mC", "N_hmC" : "N_5hmC"}, inplace=True) for df in ctrl_dfs]

    return ctrl_dfs

def calculate_percentages(df, cols):
    if type(cols) == list:
        for col in cols:
            df.eval(f"percentMeth_{col.split('_')[1]} = ({col}/readCount)*100", inplace=True)
    else:
        df.eval(f"percentMeth_{cols.split('_')[1]} = ({cols}/readCount)*100", inplace=True)

    return df

def merge_positions(dfs, calculate_percentage=False, cols=None, drop=None):
    merged = (pd.concat(dfs)
              .groupby(["Chromosome", "Start", "End"], observed=False, as_index=False, sort=False)
              .sum(numeric_only=True))

    # Replace with two functions/methods that do this
    if calculate_percentage:
        merged = calculate_percentages(merged, cols)

    if drop:
        merged = merged.drop(columns=drop)

    return merged
