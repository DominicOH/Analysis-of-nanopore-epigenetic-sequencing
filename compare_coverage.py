"""
Will read the table files to compare genomic coverage depth of multiple datasets. 

Takes a while to run as it needs to generate new intermediate data files for each dataset. 
"""

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
from AnalysisTools import common, annotation_features
from AnalysisTools.helpers import timer
import concurrent.futures
import numpy as np
import pandas as pd
import string
import pyranges as pr
import argparse
import gc

def fetch_data(dry_run: bool):
    """
    Fetches modkit/bismark data. Order of return is: nanopore, tab, oxbs. 

    ::bool dry_run:: Whether to return data for internal testing.  
    """    
    if dry_run:
        nano_path = "data/dryruns/coverage_plot/nanopore/"
        oxbs_path = "data/dryruns/coverage_plot/oxbs/"
        tab_path = "data/dryruns/coverage_plot/tab/"

    else:
        nano_path = "data/modbases/coverage_figure_beds/nanopore/"
        oxbs_path = "data/modbases/coverage_figure_beds/oxbs/"
        tab_path = "data/modbases/coverage_figure_beds/tab/"

    with concurrent.futures.ThreadPoolExecutor(3) as load_executor:
        nanopore, oxbs, tab = load_executor.map(lambda path: common.fetch_modbeds(path, ["readCount"]), [nano_path, oxbs_path, tab_path])

    nanopore = pd.concat([df.assign(Technique = "Nanopore", Replicate = f"Rep. {i}") for i, df in enumerate(nanopore)], copy=False)
    oxbs = pd.concat([df.assign(Technique = "oxBS", Replicate = f"Rep. {i}") for i, df in enumerate(oxbs)], copy=False)
    tab = pd.concat([df.assign(Technique = "TAB", Replicate = f"Rep. {i}") for i, df in enumerate(tab)], copy=False)
    
    return nanopore, oxbs, tab

def onlyAutosomal(df):
    df = df.loc[df.loc[:, "Chromosome"].str.match("^(chr)\d+$")]
    return df

def make_histogram(df, ax, color, **kwargs):
    return sns.histplot(df, x="readCount", ax=ax, color=color, **kwargs)  

def depth_to_median(df):
    median = df["readCount"].median()
    df = df.eval("readCount_vs_avg = ((readCount - @median)/@median)*100", local_dict={"median" : median})
    df["readCount_vs_avg"] = df["readCount_vs_avg"].astype(np.int16)
    return df

def prepare_boxplot(df):
    df = (df.replace(["intergenic", "genes", "intron", "cds", "1kbPromoter"], 
                   ["Intergenic", "Genic", "Intron", "CDS", "Promoter"])
            .drop(columns=["Chromosome", "Start", "End", "Replicate", "Start_Feature", "End_Feature", "readCount", "Strand"], 
                errors="ignore"))
    
    return df

def annotated_boxplot(df):
    feature_dir_path = "/mnt/data1/doh28/analyses/mouse_hydroxymethylome_analysis/feature_references/coverage_comparison_features/"
    annotator = annotation_features.Annotator(feature_dir_path)

    df = depth_to_median(df)
    df = annotator.annotate(df)

    df = prepare_boxplot(df)
    return df

@timer
def fig_main(dry_run):
    sns.set_style("ticks")
    mpl.rc('font', size=5)

    fig = plt.figure(figsize=(89/25.4, 60/25.4), dpi=600, layout="constrained")
    gs = GridSpec(2, 3, fig)

    print("Loading datasets...")
    nanopore, oxbs, tab = fetch_data(dry_run)
        
    print("Counting number of CpG sites")
    autosomal_dataframes = pd.concat([nanopore, tab, oxbs], copy=False)
    LEN_UNION = len(pr.PyRanges(autosomal_dataframes).merge(False, slack=-1)) # len of union of covered sites saved for later

    # Unique coverage comparison # 

    coverage_by_rep = autosomal_dataframes.loc[:, ("Replicate", "Technique")].value_counts(dropna=True).reset_index(name="Count")
    del autosomal_dataframes

    coverage_by_rep = coverage_by_rep.assign(Proportion_of_union = lambda r: (r["Count"]/LEN_UNION)*100)
    ax4 = fig.add_subplot(gs[1, 0])

    print("Making barplot of coverage breadth by replicate...")
    sns.barplot(coverage_by_rep, 
                x="Technique", y="Proportion_of_union", 
                hue="Technique", palette=sns.color_palette("Greens", 4)[1:],
                order=["oxBS", "TAB", "Nanopore"], 
                hue_order=["oxBS", "TAB", "Nanopore"], 
                errorbar=("sd", 1), err_kws={
                    "linewidth" : 0.8, 
                    "solid_capstyle" : "butt"
                    },
                capsize=0.5,
                ax=ax4)
    print("Done...")
    
    ax4.set_ylabel("Percentage of all\ncovered sites")
    ax4.set_xlabel(None)
    ax4.tick_params("x", labelrotation=15)

    # Coverage depth histograms # 

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])

    print("Making histograms of coverage depth by replicate...")
    with concurrent.futures.ThreadPoolExecutor(4) as hist_executor:
        hist_futures = [hist_executor.submit(make_histogram, df, ax, color, stat="percent", multiple="dodge", binrange=(0, 32), binwidth=2) for df, ax, color in zip(
            [oxbs, tab, nanopore],
            [ax1, ax2, ax3],
            [sns.color_palette("Greens", 4)[1], sns.color_palette("Greens", 4)[2], sns.color_palette("Greens", 4)[3]])]
                
        for future in hist_futures:
            future.result()
    
    del hist_futures
    print("Done...")
        
    ax1.set_title("oxBS-seq")
    ax2.set_title("TAB-seq")
    ax3.set_title("Nanopore")

    for index, ax in enumerate([ax1, ax2, ax3]):
        ax.set_title(string.ascii_lowercase[index], fontdict={"fontweight" : "bold"}, loc="left")
        ax.set_xlabel("Depth at CpG site")
        ax.set_ylabel("Percentage of sites")

    # Coverage depth by annotation # 

    gc.collect()

    ax5 = fig.add_subplot(gs[1, 1:])

    print("Preparing annotated boxplot data...")
    # with concurrent.futures.ThreadPoolExecutor(3) as annotation_executor:
    annotated_boxplot_futures = map(annotated_boxplot, [nanopore, tab, oxbs])
    print("Done. Joining...")
        
    boxplot_df = pd.concat([df for df in annotated_boxplot_futures], copy=False).reset_index(drop=True)
    del nanopore, tab, oxbs, annotated_boxplot_futures
    
    gc.collect()    
    print(f"Done. Joined {len(boxplot_df)} rows.")
    
    print("Making boxplot of coverage by context...")

    ax5.axhline(0, ls="--", c="k", lw=0.8, label="Median")
    sns.boxplot(boxplot_df, 
                x="feature_type", 
                y="readCount_vs_avg", 
                hue="Technique",
                hue_order=["oxBS", "TAB", "Nanopore"],
                order=["Intergenic", "Genic", "Intron", "CDS", "Promoter", "CGI"],
                gap=0.2,
                showfliers=False,
                palette=sns.color_palette("Greens", 4)[1:], 
                linewidth=0.8,
                ax=ax5)
    
    print(boxplot_df.groupby(["Technique", "feature_type"])["readCount_vs_avg"].median())
    
    print("Done.")
    ax5.set_ylabel("Difference to median (%)")
    ax5.yaxis.set_ticks(np.arange(-100, 400, 100))

    ax5.set_xlabel(None)
    ax5.tick_params("x")

    sns.move_legend(ax5, "lower center", bbox_to_anchor=(.5, 1), 
                    ncols=3, title=None, frameon=False, alignment="center")
    
    for index, ax in enumerate([ax4, ax5]):
        ax.set_title(string.ascii_lowercase[index+3], fontdict={"fontweight" : "bold"}, loc="left")

    sns.despine()

    if dryrun:
        outpath = "plots/tests/compare_coverage.png"
    else:
        outpath = "plots/compare_coverage.png"
    return fig.savefig(outpath)

##### main function ##### 

if __name__=="__main__":
    parser = argparse.ArgumentParser(
                        prog = "compare_coverage",
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
    print("Completed.")