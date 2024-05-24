from sklearn.metrics import mean_squared_error
import concurrent.futures
from AnalysisTools.common import fetch_modbeds, calculate_percentages
import numpy as np
import argparse
from AnalysisTools.helpers import timer

def inner_comparison(a_df, b_df, a_mod, b_mod):
        with concurrent.futures.ThreadPoolExecutor(2) as tpe:
                futures = [tpe.submit(calculate_percentages, df, col) for df, col in zip([a_df, b_df], [f"N_{a_mod}", f"N_{b_mod}"])]
                a_df, b_df = [future.result().drop(columns=["readCount"]) for future in futures]
                merge = (a_df.merge(b_df, 
                                    on=["Chromosome", "Start", "End"], 
                                    how="inner", 
                                    suffixes=["_A", "_B"])
                                .drop(columns=["Chromosome", "Start", "End", ]))
            # Assuming test_data and control_data are arrays of predicted and actual values
            # Calculate the error metrics
        if a_mod == b_mod:
            mse = mean_squared_error(merge[f"percentMeth_{a_mod}_A"], merge[f"percentMeth_{b_mod}_B"])
        else:
            mse = mean_squared_error(merge[f"percentMeth_{a_mod}"], merge[f"percentMeth_{b_mod}"])
        
        rmse = np.sqrt(mse)
        if rmse != 0: # Remove self comparisons
            return rmse
        else: 
            return np.nan

def rmse_comparer(a, b, a_mod, b_mod):
    rmse_values = []
    for a_df in a:
        with concurrent.futures.ProcessPoolExecutor(4) as outer_tpe:
            rep_futures = [outer_tpe.submit(inner_comparison, a_df, b_df, a_mod, b_mod) for b_df in b]
            rmse_values.append([rmse.result() for rmse in rep_futures])
            
    return rmse_values

def nested_mean(ls):
    return np.mean([np.nanmean(nested_ls) for nested_ls in ls])

@timer
def main(dryrun=True):
    if dryrun:
        nano_path = "data/dryruns/modbeds/nanopore/"
        tab_path = "data/dryruns/modbeds/tab/"
        oxbs_path = "data/dryruns/modbeds/oxbs/"

    else:
        nano_path = "data/modbases/modbeds/nanopore/"
        tab_path = "data/modbases/modbeds/tab/"
        oxbs_path = "data/modbases/modbeds/oxbs/"

    paths = [nano_path, tab_path, oxbs_path]
    cols = [["readCount", "N_mC", "N_hmC"], ["readCount", "N_mod"], ["readCount", "N_mod"]]

    with concurrent.futures.ProcessPoolExecutor(3) as fetch_executor:
        modbeds_future = [fetch_executor.submit(fetch_modbeds, path, cols) for path, cols in zip(paths, cols)]
        nano_mb, tab_mb, oxbs_mb = [future.result() for future in modbeds_future]
    
    # print(f"Interassay RMSE between Nanopore 5mC: {nested_mean(rmse_comparer(nano_mb, nano_mb, '5mC', '5mC'))}")
    # print(f"Interassay RMSE between Nanopore 5hmC: {nested_mean(rmse_comparer(nano_mb, nano_mb, '5hmC', '5hmC'))}")

    # print(f"RMSE between Nanopore 5mC and oxBS: {nested_mean(rmse_comparer(nano_mb, oxbs_mb, '5mC', 'mod'))}")
    # print(f"RMSE between Nanopore 5hmC and TAB: {nested_mean(rmse_comparer(nano_mb, tab_mb, '5hmC', 'mod'))}")
        
    print(f"Interassay RMSE between oxBS 5mC: {nested_mean(rmse_comparer(oxbs_mb, oxbs_mb, 'mod', 'mod'))}")
    print(f"Interassay RMSE between TAB 5hmC: {nested_mean(rmse_comparer(tab_mb, tab_mb, 'mod', 'mod'))}")
    
    return print("Done.") 

##### Main function #####

if __name__=="__main__":
    parser = argparse.ArgumentParser(
                        prog = "cpg_methylation_compare",
                        description = "Compares the coverage of the different datasets.")
    parser.add_argument("-d ", "--dryrun", 
                        action="store_true", 
                        dest="dryrun", 
                        required=False,
                        help="Whether a test output is produced.") 

    args = parser.parse_args()

    if args.dryrun:
        dryrun = True
    else: 
        dryrun = False

    main(dryrun)    