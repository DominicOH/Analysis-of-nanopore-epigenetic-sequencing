import argparse
import pandas as pd
import numpy as np
import gc
from imblearn.over_sampling import ADASYN
from sklearn.preprocessing import Binarizer
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from AnalysisTools.helpers import timer
from AnalysisTools.common import *
import string


def prep_dev_plot(df, bs_method, mod):
    df = (df.reset_index()
            .drop(columns=["Chromosome", "Start", "End"])
            .eval(f"diff = `{bs_method} mean` - `Nanopore mean`")
            .assign(Mod = mod))
    return df
@timer
def fig_main(dryrun=True):
    fig = plt.figure(figsize=(89/25.4, 89/25.4), dpi=600, layout="constrained")

    sns.set_style("ticks")
    mpl.rc('font', size=5)

    gs = GridSpec(1, 3, fig)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])

    # Initial data collection # 
    # Should aim to generalise fetch_modbed data loading

    print("Loading data")
    
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
    gc.collect()
    
    with concurrent.futures.ProcessPoolExecutor(3) as merge_executor:
        merge_futures = [merge_executor.submit(merge_positions, modbed, True, col) for modbed, col in zip(modbeds, cols)]
        nanopore_average, tab_average, ox_average = [future.result().assign(Method = method) for future, method in zip(merge_futures, 
                                                                                                                       ["Nanopore mean", "TAB mean", "oxBS mean"])]
    ox_average = ox_average.rename(columns={"percentMeth_mod" : "percentMeth_5mC"})
    tab_average = tab_average.rename(columns={"percentMeth_mod" : "percentMeth_5hmC"})        

    print("Loaded data")

    # General distribution comparison # 

    dfs = [nanopore_average, ox_average, nanopore_average, tab_average]
    xs = ["percentMeth_5mC", "percentMeth_5mC", "percentMeth_5hmC", "percentMeth_5hmC"]
    cs = [sns.color_palette("Greens_r", 4)[0], sns.color_palette("Greens_r", 4)[0], sns.color_palette("Blues_r", 4)[0], sns.color_palette("Blues_r", 4)[0]]
    lss = [":", "-", ":", "-"]
    labels = ["Nanopore 5mC", "oxBS 5mC", "Nanopore 5hmC", "TAB 5hmC"]
    axs = [ax1, ax1, ax2, ax2]

    print("Plotting ECDFs")

    with concurrent.futures.ThreadPoolExecutor(4) as ecdf_executor:
        ecdf_futures = [ecdf_executor.submit(sns.ecdfplot, data=df, x=x, c=c, lw=.8, ls=ls, label=label, ax=ax) 
                        for df, x, c, ls, label, ax in zip(dfs, xs, cs, lss, labels, axs)]
        [future.result() for future in ecdf_futures]

    ax1.set_xlabel("Site 5mC (%)")
    ax2.set_xlabel("Site 5hmC (%)")
    
    nano_oxbs = pd.concat([nanopore_average, ox_average], join="inner", copy=False).pivot_table(values="percentMeth_5mC", index=["Chromosome", "Start", "End"], columns="Method").dropna()
    del ox_average
    
    nano_tab = pd.concat([nanopore_average, tab_average], join="inner", copy=False).pivot_table(values="percentMeth_5hmC", index=["Chromosome", "Start", "End"], columns="Method").dropna()
    del tab_average, nanopore_average

    gc.collect()

    for ax in [ax1, ax2]:
        ax.legend() 
        ax.set_aspect(100)
        ax.set_ylabel("Cumulative proportion")

    sns.move_legend(ax1, "upper left", frameon=False)
    sns.move_legend(ax2, "lower right", frameon=False)

    print("Plotting deviation")
    with concurrent.futures.ProcessPoolExecutor(2) as dev_plot_executor:
        dev_plot_futures = [dev_plot_executor.submit(prep_dev_plot, df, bs_method, mod) for df, bs_method, mod, in zip([nano_oxbs, nano_tab], 
                                                                                                                       ["oxBS", "TAB"], 
                                                                                                                       ["5mC", "5hmC"])]
        nano_oxbs, nano_tab = [future.result() for future in dev_plot_futures]

    gc.collect()
    sns.histplot(pd.concat([nano_oxbs, nano_tab]),
                x="diff", 
                stat="proportion",
                binrange=(-50, 50), binwidth=2,
                element="step", fill=False, 
                lw=0.8,
                hue="Mod", palette="GnBu",
                ax=ax3)
    
    ax3.set_ylim((0, 0.1))
    ax3.set_xlim((-50, 50))
    ax3.set_aspect(1000)
    ax3.set_xlabel("Bisulphite % - Nanopore %")

    sns.move_legend(ax3, "upper right", frameon=False, title=None)
    sns.despine()

    print("Done. Saving")

    for index, ax in enumerate(fig.axes):
        ax.set_title(f"{string.ascii_lowercase[index]}", fontdict={"fontweight" : "bold"}, loc="left")

    if dryrun:
        outpaths = ["plots/tests/cpg_methylation_compare.png"]
    else:
        outpaths = ["plots/cpg_methylation_compare.png"]
        outpaths = ["plots/cpg_methylation_compare.svg"]

    return [fig.savefig(path) for path in outpaths]

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

    fig_main(dryrun)    