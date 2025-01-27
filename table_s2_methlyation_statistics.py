import argparse
from AnalysisTools.common import *
from AnalysisTools.helpers import timer
from scipy import stats
import concurrent.futures
import numpy as np
import pingouin as pg

def averages(modbeds, col):
        means = []
        for mb in modbeds:
            mean = round(mb[col].mean(), 2)
            means.append(mean)
        
        return means

def assign_replicates(modbeds, technique):
    modbeds = [df.assign(Technique = technique, Replicate = f"Rep. {i+1}") 
                                for i, df in enumerate(modbeds)]
    return modbeds

def count_mods(modbeds, col):
    return [df[col].sum()/df["readCount"].sum() for df in modbeds]

def count_hypo(modbeds, col):
    return [len(df.query(f"{col} == 0"))/len(df) for df in modbeds]

@timer
def stats_main(dryrun):
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
        modbeds = [future.result() for future in modbeds_future]

    cols = [["N_5mC", "N_5hmC"], ["N_mod"], ["N_mod"]]
    
    with concurrent.futures.ThreadPoolExecutor(3) as merge_executor:
        merge_futures = [merge_executor.submit(merge_positions, modbed, True, col) for modbed, col in zip(modbeds, cols)]
        nanopore_average, tab_average, ox_average = [future.result().assign(Method = method) for future, method in zip(merge_futures, ["Nanopore mean", "TAB mean", "oxBS mean"])]

    ox_average = ox_average.rename(columns={"percentMeth_mod" : "percentMeth_5mC"})
    tab_average = tab_average.rename(columns={"percentMeth_mod" : "percentMeth_5hmC"})  
    
    nano_modbeds, tab_modbeds, oxbs_modbeds = [assign_replicates(mb, technique) for mb, technique in zip(modbeds, ["Nanopore", "TAB", "oxBS"])]

    nano_modbeds = [calculate_percentages(df, ["N_5mC", "N_5hmC"]) for df in nano_modbeds]
    oxbs_modbeds = [calculate_percentages(df, "N_mod").rename(columns={"percentMeth_mod" : "percentMeth_5mC"}) for df in oxbs_modbeds]
    tab_modbeds = [calculate_percentages(df, "N_mod").rename(columns={"percentMeth_mod" : "percentMeth_5hmC"}) for df in tab_modbeds]
    
    x, y = averages(nano_modbeds, "percentMeth_5mC"), averages(oxbs_modbeds, "percentMeth_5mC")
    print("Nanopore normality:", pg.normality(x), "\noxBS normality:", pg.normality(y))
    test = pg.ttest(x, y, paired=False)
    print("5mC averages:", test)

    x, y = averages(nano_modbeds, "percentMeth_5hmC"), averages(tab_modbeds, "percentMeth_5hmC")
    print("Nanopore normality:", pg.normality(x), "\nTAB normality:", pg.normality(y))
    test = pg.ttest(x, y, paired=False)
    print("5hmC averages:", test)

    x, y = count_hypo(nano_modbeds, "percentMeth_5hmC"), count_hypo(tab_modbeds, "percentMeth_5hmC")
    print("Nanopore normality:", pg.normality(x), "\nTAB normality:", pg.normality(y))
    test = pg.ttest(x, y, paired=False)
    print("5hmC hypo:", test)

    nano_oxbs = pd.concat([nanopore_average, ox_average], join="inner").pivot_table(values="percentMeth_5mC", index=["Chromosome", "Start", "End"], columns="Method").dropna()
    nano_tab = pd.concat([nanopore_average, tab_average], join="inner").pivot_table(values="percentMeth_5hmC", index=["Chromosome", "Start", "End"], columns="Method").dropna()
    
    print(f"5mC cor in {len(nano_oxbs)} positions: ", stats.spearmanr(nano_oxbs["Nanopore mean"], nano_oxbs["oxBS mean"]))
    print(f"5hmC cor in {len(nano_tab)} positions: ", stats.spearmanr(nano_tab["Nanopore mean"], nano_tab["TAB mean"]))
    return ... 

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

    stats_main(dryrun)  