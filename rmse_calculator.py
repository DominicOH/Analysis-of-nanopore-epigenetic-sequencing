from sklearn.metrics import root_mean_squared_error
import pandas as pd
import subprocess
import concurrent.futures
from AnalysisTools import common
import numpy as np
import argparse
from AnalysisTools.helpers import timer

def one_vs_one_rmsd(nanopore_df: pd.DataFrame, bis_df: pd.DataFrame, mod: str, depth_levels: list[int]):
    merge = (nanopore_df.merge(bis_df, 
                               on=["Chromosome", "Start", "End"], 
                               how="inner")
                        .drop(columns=["Chromosome", "Start", "End"]))
    
    rmses = []
    sizes = []
    means_n = []
    means_b = []

    for level in depth_levels:
        merge = merge.query(f"readCount_x >= {level} & readCount_y >= {level}")
        try: 
            rmse = root_mean_squared_error(merge[f"percentMeth_{mod}"], merge[f"percentMeth_mod"])
        except:
            if len(merge) == 0:
                # print("Too few matched positions.")
                pass

        rmses.append(rmse)
        sizes.append(len(merge))
        means_n.append(merge[f"percentMeth_{mod}"].mean())
        means_b.append(merge[f"percentMeth_mod"].mean())

    ovo_result = pd.DataFrame({
        "Depth" : depth_levels,
        "Size" : sizes,
        "RMSD" : rmses, 
        "Nanopore_mean" : means_n,
        "Bisulphite_mean" : means_b
    })

    return  ovo_result            

def one_vs_many_rmsd(nano_df: pd.DataFrame, bs_dfs: list[pd.DataFrame], mod: str, depth_levels: list[int]):
    """
    Aggregates the result of several one-vs-one RMSD calculations. 
    """
    with concurrent.futures.ThreadPoolExecutor(len(bs_dfs)) as ovo_tpe:
        futures = [ovo_tpe.submit(one_vs_one_rmsd, nano_df, bs_df, mod, depth_levels) for bs_df in bs_dfs]
        one_vs_many_result = pd.concat([rmse.result().assign(Bisulphite_rep = i+1) for i, rmse in enumerate(futures)])

    return one_vs_many_result

def many_vs_many_rmsd(nano_dfs: list[pd.DataFrame], bs_dfs: list[pd.DataFrame], mod: str, depth_levels: list[int]):
    """
    Aggregates the result of several one-vs-many RMSD calculations. Each one-vs-many is numbered. 
    """
    with concurrent.futures.ProcessPoolExecutor(len(nano_dfs)) as mvm_executor:
        ovm_futures = [mvm_executor.submit(one_vs_many_rmsd, nano_df, bs_dfs, mod, depth_levels) for nano_df in nano_dfs]
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
    nano_path = "data/modbases/modbeds/nanopore/"
    tab_path = "data/modbases/modbeds/tab/"
    oxbs_path = "data/modbases/modbeds/oxbs/"

    paths = [nano_path, tab_path, oxbs_path]
    cols = [["readCount", "N_mC", "N_hmC"], ["readCount", "N_mod"], ["readCount", "N_mod"]]

    with concurrent.futures.ProcessPoolExecutor(3) as fetch_executor:
        data_future = [fetch_executor.submit(load_and_calculate_percentages, path, cols) for path, cols in zip(paths, cols)]
        nano_mb, tab_mb, oxbs_mb = [future.result() for future in data_future]
    
    # print(f"Interassay RMSE between Nanopore 5mC: {nested_mean(rmse_comparer(nano_mb, nano_mb, '5mC', '5mC'))}")
    # print(f"Interassay RMSE between Nanopore 5hmC: {nested_mean(rmse_comparer(nano_mb, nano_mb, '5hmC', '5hmC'))}")

    print("Performing many-vs-many RMSD calculations for Nanopore vs. oxBS")
    ox_result = many_vs_many_rmsd(nano_mb, oxbs_mb, '5mC', [i for i in range(5, 25, 1)])
    print(f"Average RMSD at 10x depth: {ox_result.loc[ox_result['Depth']==10, 'RMSD'].mean()}")
    ox_result.to_csv("nanopore_vs_ox_rmsd.tsv", index=False, sep="\t")

    print("Performing many-vs-many RMSD calculations for Nanopore vs. TAB")
    tab_result = many_vs_many_rmsd(nano_mb, tab_mb, '5hmC', [i for i in range(5, 25, 1)])
    print(f"Average RMSD at 10x depth: {tab_result.loc[tab_result['Depth']==10, 'RMSD'].mean()}")
    tab_result.to_csv("nanopore_vs_tab_rmsd.tsv", index=False, sep="\t")
    # print(f"RMSE between Nanopore 5hmC and TAB: {}")
        
    # print(f"Interassay RMSE between oxBS 5mC: {nested_mean(rmse_comparer(oxbs_mb, oxbs_mb, 'mod', 'mod'))}")
    # print(f"Interassay RMSE between TAB 5hmC: {nested_mean(rmse_comparer(tab_mb, tab_mb, 'mod', 'mod'))}")
    
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

    args = parser.parse_args()

    global test
    if args.test:
        test = True
    else: 
        test = False

    main()    