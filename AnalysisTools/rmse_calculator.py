from sklearn.metrics import root_mean_squared_error, median_absolute_error
import pandas as pd
import concurrent.futures
from AnalysisTools import common
import numpy as np
import argparse
from AnalysisTools.helpers import timer

"""
This script is used to calculate Root Mean Square Deviation between various permutations of the Nanopore and Bisulphite data. 

The dataframe inputs used are produced using AnalysisTools/read_modbed using a Modkit bedMethyl or Bismark zero.cov output. 
"""

def levelled_rmsd(merged_dataframe: pd.DataFrame, mod, depth_levels: range, intraassay: bool = False):
    level_dicts = []
    comparison = ""
    for level in depth_levels:
        level_df = merged_dataframe.loc[merged_dataframe.eval(f"readCount_x >= {level} & readCount_y >= {level}")]

        if len(level_df) != 0:
            if not intraassay:
                x = level_df[f"percentMeth_{mod}"].to_numpy()
                y = level_df["percentMeth_mod"].to_numpy()
                rmse = root_mean_squared_error(x, y)
                mad = median_absolute_error(x, y)
            else:
                x = level_df[f"percentMeth_{mod}_x"].to_numpy()
                y = level_df[f"percentMeth_{mod}_y"].to_numpy()
                rmse = root_mean_squared_error(x, y)
                mad = median_absolute_error(x, y)
            
        else:
            pass
        
        level_dict = {
            "Depth" : level,
            "Size" : len(level_df),
            "RMSD" : rmse, 
            "MAD" : mad,
            "Mean_X" : x.mean(),
            "Mean_Y" : y.mean(),
        }   
        if intraassay:
            if not comparison:
                comparison = merged_dataframe.get("Comparison")[0]
            else:
                level_dict.update({
                    "Comparison" : comparison
                })

        del level_df
        level_dicts.append(level_dict)
        
        if debug:
            print(level_dict)

    out_df = pd.DataFrame.from_records(level_dicts)

    return out_df

def one_vs_one(nanopore_df: pd.DataFrame, other_df: pd.DataFrame, mod: str, depth_levels: range, intrassay: bool = False):
    merge = (nanopore_df.merge(other_df, 
                               on=["Chromosome", "Start", "End"], 
                               how="inner")
                        .drop(columns=["Chromosome", "Start", "End", "N_mod", "N_5hmC", "N_5mC"], errors="ignore"))
    
    ovo_result = levelled_rmsd(merge, mod, depth_levels, intrassay)

    if debug:
        print(ovo_result)

    return  ovo_result  

def one_vs_many(nano_df: pd.DataFrame, other_dfs: list[pd.DataFrame], mod: str, depth_levels: range, other_rep_names, intraassay: bool = False):
    """
    Aggregates the result of several one-vs-one RMSD calculations. 
    """
    if type(other_rep_names) == list[int]:
        other_rep_names = [i+1 for i in other_rep_names]
    with concurrent.futures.ThreadPoolExecutor(len(other_dfs)) as ovo_tpe:
        futures = [ovo_tpe.submit(one_vs_one, nano_df, other_df, mod, depth_levels, intraassay) for other_df in other_dfs]
        one_vs_many_result = pd.concat([rmse.result().assign(Other_rep = i) for i, rmse in zip(other_rep_names, futures)])

    return one_vs_many_result

def many_vs_many_intrassay(compare_dfs: list[pd.DataFrame], mod: str, depth_levels: range, rep_names: list = None):
    if rep_names is None:
        rep_names = [i+1 for i in range(len(compare_dfs))]

    if debug:
        processors = 1
    else: 
        processors = len(compare_dfs)

    with concurrent.futures.ThreadPoolExecutor(processors) as mvm_executor:
        mvm_futures = [mvm_executor.submit(one_vs_many, 
                                           df_a.assign(Comparison = i), compare_dfs, mod, depth_levels, rep_names, intraassay=True) 
                       for i, df_a in zip(rep_names, compare_dfs)]
        mvm_result = pd.concat([future.result() for future in mvm_futures])
    
    mvm_result = (mvm_result.replace({"RMSD" : 0}, np.nan)
                  .dropna()
                  .drop_duplicates(["Depth", "Size", "RMSD"]))
                  
    return mvm_result 

def many_vs_many_interassay(nano_dfs: list[pd.DataFrame], bs_dfs: list[pd.DataFrame], mod: str, depth_levels: list[int]):
    """
    Aggregates the result of several one-vs-many RMSD calculations. Each one-vs-many is numbered. 
    """
    rep_names = [i for i in range(len(bs_dfs))]
    with concurrent.futures.ThreadPoolExecutor(len(nano_dfs)) as mvm_executor:
        ovm_futures = [mvm_executor.submit(one_vs_many, nano_df, bs_dfs, mod, depth_levels, rep_names) for nano_df in nano_dfs]
        mvm_result = pd.concat([ovm.result().assign(Comparison = i) for i, ovm in zip(("CBM2_1", "CBM2_2", "CBM3_1", "CBM3_2"), ovm_futures)])

    return mvm_result

def load_and_calculate_percentages(path, cols: list[str]):
    modbeds = common.fetch_modbeds(path, cols, test)

    cols.remove("readCount")
    if "N_mC" in cols:
        cols.remove("N_mC")
        cols.append("N_5mC")
        cols.remove("N_hmC")
        cols.append("N_5hmC")

    with concurrent.futures.ThreadPoolExecutor(len(modbeds)) as tpe: 
        modbeds = [tpe.submit(common.calculate_percentages, df, cols) for df in modbeds]
        modbeds = [future.result() for future in modbeds]

    return modbeds

@timer
def main():
    nano_path = "../data/modbases/modbeds/nanopore_base5/"
    tab_path = "../data/modbases/modbeds/tab/"
    oxbs_path = "../data/modbases/modbeds/oxbs/"

    paths = [nano_path, tab_path, oxbs_path]
    cols = [["readCount", "N_mC", "N_hmC"], ["readCount", "N_mod"], ["readCount", "N_mod"]]

    print("Loading datasets")
    with concurrent.futures.ProcessPoolExecutor(3) as fetch_executor:
        data_future = [fetch_executor.submit(load_and_calculate_percentages, path, cols) for path, cols in zip(paths, cols)]
        nano_mb, tab_mb, oxbs_mb = [future.result() for future in data_future]
    
    print("Performing many-vs-many RMSD calculations for Nanopore vs. Nanopore")

    level_range = range(4, 25, 1)
    nano_iassay_mc = many_vs_many_intrassay(nano_mb, "5mC", level_range, ("CBM2_1", "CBM2_2", "CBM3_1", "CBM3_2"))
    nano_iassay_hmc = many_vs_many_intrassay(nano_mb, "5hmC", level_range, ("CBM2_1", "CBM2_2", "CBM3_1", "CBM3_2"))
    nano_iassay_mc.to_csv("data/rmsd/nanopore_intrassay_5mC.tsv", index=False, sep="\t")
    nano_iassay_hmc.to_csv("data/rmsd/nanopore_intrassay_5hmC.tsv", index=False, sep="\t")
    
    print("Performing many-vs-many RMSD calculations for Bisulphite vs. Bisulphite")
    tab_iassay = many_vs_many_intrassay(tab_mb, "mod", level_range)
    tab_iassay.to_csv("data/rmsd/tab_intrassay.tsv", index=False, sep="\t")
    ox_iassay = many_vs_many_intrassay(oxbs_mb, "mod", level_range)
    ox_iassay.to_csv("data/rmsd/oxbs_intrassay.tsv", index=False, sep="\t")

    print("Performing many-vs-many RMSD calculations for Nanopore vs. oxBS")
    ox_result = many_vs_many_interassay(nano_mb, oxbs_mb, '5mC', level_range)
    print(f"Average RMSD at 10x depth: {ox_result.loc[ox_result['Depth']==10, 'RMSD'].mean()}")
    ox_result.to_csv("data/rmsd/nanopore_vs_ox_rmsd.tsv", index=False, sep="\t")

    print("Performing many-vs-many RMSD calculations for Nanopore vs. TAB")
    tab_result = many_vs_many_interassay(nano_mb, tab_mb, '5hmC', level_range)
    print(f"Average RMSD at 10x depth: {tab_result.loc[tab_result['Depth']==10, 'RMSD'].mean()}")
    tab_result.to_csv("data/rmsd/nanopore_vs_tab_rmsd.tsv", index=False, sep="\t")
    
    return print("Done.") 

##### Main function #####

if __name__=="__main__":
    parser = argparse.ArgumentParser(
                        prog = "cpg_methylation_compare",
                        description = "Compares the coverage of the different datasets.")
    parser.add_argument("-t ", "--test_run", 
                        action="store_true", 
                        dest="test", 
                        required=False,
                        help="Whether a test output is produced.") 
    parser.add_argument("--debug", 
                        action="store_true", 
                        dest="debug", 
                        default=False,
                        help="Whether a test output is produced.") 

    args = parser.parse_args()
    global test, debug
    test, debug = args.test, args.debug

    main()    