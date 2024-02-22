import pandas as pd
from concurrent import futures
from AnalysisTools import read_modbed
import os 
import subprocess
import concurrent.futures

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

def read_table(path, usecols):
    default_usecols = ["Chromosome", "Start", "End"]

    if type(usecols) == list:
        default_usecols.extend(usecols)
    else: 
        default_usecols.append(usecols)

    df = pd.read_table(path, sep="\t", 
                        usecols=default_usecols)
    return df

def fetch_nanopore(usecols, dryrun=True):
    if dryrun:
        root_path = "data/dryruns/modbeds/"

        cbm2_1_path = root_path + "CBM_2_rep1.masked.bed.modbed"
        cbm2_2_path = root_path + "CBM_2_rep2.masked.bed.modbed"

        cbm3_1_path = root_path + "CBM_3_rep2.masked.bed.modbed"
        cbm3_2_path = root_path + "CBM_3_rep2.masked.bed.modbed"

    else:
        root_path = "data/modbases/modbeds/"

        cbm2_1_path = root_path + "CBM_2_rep1.modbed"
        cbm2_2_path = root_path + "CBM_2_rep2.modbed"

        cbm3_1_path = root_path + "CBM_3_rep1.modbed"
        cbm3_2_path = root_path + "CBM_3_rep2.modbed"
    
    return map(lambda path: read_table(path, usecols).rename(columns={"N_mC" : "N_5mC", "N_hmC" : "N_5hmC"}), [cbm2_1_path, cbm2_2_path, cbm3_1_path, cbm3_2_path])

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

def fetch_controls(usecols, dryrun=True):
    if dryrun:
        root_path = "data/dryruns/modbeds/"
        
        zymo_m1 = root_path + "zymo_wga_methylated_rep1.masked.bed.modbed"
        zymo_m2 = root_path + "zymo_wga_methylated_rep2.masked.bed.modbed"

        zymo_u1 = root_path + "zymo_wga_unmodified_rep1.masked.bed.modbed"
        zymo_u2 = root_path + "zymo_wga_unmodified_rep2.masked.bed.modbed"

    else:
        root_path = "data/modbases/modbeds/"
        
        zymo_m1 = root_path + "zymo_methylated_rep1.modbed"
        zymo_m2 = root_path + "zymo_methylated_rep2.modbed"

        zymo_u1 = root_path + "zymo_unmodified_rep1.modbed"
        zymo_u2 = root_path + "zymo_unmodified_rep2.modbed"

    pos_controls = map(lambda path: read_table(path, usecols).rename(columns={"N_mC" : "N_5mC", "N_hmC" : "N_5hmC"}), [zymo_m1, zymo_m2])
    neg_controls = map(lambda path: read_table(path, usecols).rename(columns={"N_mC" : "N_5mC", "N_hmC" : "N_5hmC"}), [zymo_u1, zymo_u2])

    return pos_controls, neg_controls

def merge_positions(dfs, cols, drop=True):
    merged = pd.concat(dfs).groupby(["Chromosome", "Start", "End"]).sum(numeric_only=True)

    if type(cols) == list:
        for col in cols:
            merged[f"percentMeth_{col.split('_')[1]}"] = (merged[col]/merged["readCount"])*100
        if drop:
            merged = merged.drop(columns=["readCount", *cols])
        return merged
    else:
        merged[f"percentMeth_{cols.split('_')[1]}"] = (merged[cols]/merged["readCount"])*100
        if drop:
            merged = merged.drop(columns=["readCount", cols])
        return merged
