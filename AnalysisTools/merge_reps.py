import argparse
import pandas as pd
from AnalysisTools import read_modbed, common
from AnalysisTools.common import timer
import time

def load_data(path, modbase=None, dry_run=False, min_depth=1):
    if dry_run:
        nrows=10
    else:
        nrows=None

    print(f"Loading {path}.")
    df = read_modbed.openReps_Parallel(path, modbase=modbase, min_depth=min_depth, include_raw_counts=True, nrows=nrows)
    print("Removing non-autosomal CpG sites.")
    df = common.onlyAutosomal(df)
    
    return df

def group_position(df): 
    assert type(df) == pd.DataFrame, "Not a dataframe."

    if "N_mC" in df.columns and "N_hmC" in df.columns:
        aggfuncs = {
            "readCount" : sum,
            "N_mC" : sum, 
            "N_hmC" : sum
        }
        percent_cols = {
            "percentMeth_5mC" : lambda df: df.apply(lambda i: (i["N_mC"]/i["readCount"])*100, axis=1),
            "percentMeth_5hmC" : lambda df: df.apply(lambda i: (i["N_hmC"]/i["readCount"])*100, axis=1)
        }
    elif "N_mC" in df.columns and not "N_hmC" in df.columns:
        aggfuncs = {
            "readCount" : sum,
            "N_mC" : sum
        }
        percent_cols = {
            "percentMeth_5mC" : lambda df: df.apply(lambda i: (i["N_mC"]/i["readCount"])*100, axis=1)
        }
    elif "N_hmC" in df.columns and not "N_mC" in df.columns:
        aggfuncs = {
            "readCount" : sum,
            "N_hmC" : sum
        }
        percent_cols = {
            "percentMeth_5hmC" : lambda df: df.apply(lambda i: (i["N_hmC"]/i["readCount"])*100, axis=1)
        }
    else:
        raise ValueError("Can't determine type of data.")
    
    print("Grouping CpG positions.")
    print(f"Aggregating on {aggfuncs.keys}.")

    grouped_positions = df.groupby(["Chromosome", "Start", "End"], observed=False).agg(aggfuncs).reset_index()
    grouped_positions = grouped_positions.assign(**percent_cols)
    return grouped_positions

@timer
def merge_reps(inpath, outpath, modbase=None, dry_run=False, min_depth=1):
    df = load_data(inpath, modbase, dry_run, min_depth)
    grouped_positions = group_position(df)

    print(f"Success. Writing to {outpath}.")

    return grouped_positions.to_csv(outpath, sep="\t", header=False, index=False)


##### main function ##### 

if  __name__=="__main__":
    parser = argparse.ArgumentParser(
                        prog = "merge_reps",
                        description = "Merges two bedMethyl or bismark files.")
    parser.add_argument("filenames", action="append", nargs="+", required=True)
    parser.add_argument("-o", "--outpath", action="store", required=True)
    parser.add_argument("--mod_type", action="store", required=False)
    parser.add_argument("-d ", "--dryrun", action="store_true", dest="dryrun", required=False, help="Whether a test output is produced.") 
    parser.add_argument("--min_depth", action="store", required=False, default=1)

    args = parser.parse_args()


    merge_reps(args.filenames, args.outpath, args.mod_type, args.dryrun, args.min_depth)
    print(f"Done.")  