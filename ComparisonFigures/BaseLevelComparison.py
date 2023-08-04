##### Imports #####

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from common import *
import argparse

##### Plotting ##### 

def plotFigure():
    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
    sns.set_style("whitegrid")

    comparison = [["percentMeth_5mC_Prom", "percentMeth_5mC_Bisulphite"], ["percentMeth_5hmC_Prom", "percentMeth_5hmC_Bisulphite"]]

    for index, ax in enumerate(axes[0]):
        sns.ecdfplot(all_dfs.melt(
            id_vars=["chromosome", "chromStart", "chromEnd"], value_vars=comparison[index], var_name="method", value_name="percentMeth_5mC"),
        x="percentMeth_5mC", hue="method", ax=ax)
        ax.set_aspect(100)

    axes[0][0].set_xlabel("5mC beta-value")
    axes[0][1].set_xlabel("5hmC beta-value")

    for index, ax in enumerate(axes[1]):
        div = make_axes_locatable(axes[1, index])
        cax = div.append_axes("right", size="5%", pad=0.1)
        sns.histplot(all_dfs, x=comparison[index][1], y=comparison[index][0], stat="proportion", binwidth=5, pthresh=0.01, cbar=True, ax=ax, cbar_ax=cax)
        ax.set_aspect("equal")

    axes[1][0].set_xlabel("oxBS 5mC beta-value")
    axes[1][0].set_ylabel("Nanopore 5mC beta-value")

    axes[1][1].set_xlabel("TAB 5hmC beta-value")
    axes[1][1].set_ylabel("Nanopore 5hmC beta-value")

    return fig

##### Args #####

if __name__=="__main__":
    parser = argparse.ArgumentParser(
                        prog = "CompareBaseLevel",
                        description = "Compares the base level modified base calls from Nanopore, oxBS, and TAB-seq data.")
    parser.add_argument('nanopore_modbed_path', type=str, help='Filepaths to the Modbam2bed bedMethyl file.')
    parser.add_argument('oxbs_modbed_path', type=str, help='Filepaths to the Bismark zero.cov output file for oxBS.')
    parser.add_argument('tab_modbed_path', type=str, help='Filepaths to the Bismark zero.cov output file for TAB.')
    parser.add_argument("-o ", "--output_path", type=str, help="Figure save path.") 

    args = parser.parse_args()
    nanopore_modbed_path = args.nanopore_modbed_path
    oxbs_modbed_path = args.oxbs_modbed_path
    tab_modbed_path = args.tab_modbed_path

    out_path = args.out_path

    ##### Load Nanopore dataset ##### 
    nanopore_ternary_data = readModbam2bed(nanopore_modbed_path, min_depth=10, 
                                       apply_max_depth=False)

    ##### Load oxBS dataset ##### 
    oxbs_df = readBismarkZeroCov(oxbs_modbed_path, mod="5mC", 
                             min_depth=10, apply_max_depth=False)

    ##### Load TAB dataset #####
    tab_df = readBismarkZeroCov(tab_modbed_path,
                            mod="5hmC", min_depth=10, apply_max_depth=False)

    # merge datasets on position
    bis_all = tab_df.merge(oxbs_df, "inner", ["chromosome", "chromStart", "chromEnd"], suffixes=["_TAB", "_oxBS"])
    all_dfs = bis_all.merge(nanopore_ternary_data, "inner", on=["chromosome", "chromStart", "chromEnd"], suffixes=["_Bisulphite", "_Nanopore"])

    ##### Make fig #####
    fig = plotFigure()
    fig.savefig(out_path, dpi=600, format="png", bbox_inches="tight")