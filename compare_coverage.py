"""
Will read the table files to compare genomic coverage depth of multiple datasets. 

Takes a while to run as it needs to generate new intermediate data files for each dataset. 
"""

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
from AnalysisTools import read_modbed
from AnalysisTools import annotation_features 
from AnalysisTools.helpers import timer
from AnalysisTools import common  
import concurrent.futures
import numpy as np
import pandas as pd
import string
import pyranges as pr
import argparse

def fetch_data(dry_run: bool, split_biorep=False, **kwargs):
    """
    Fetches modkit/bismark data. Order of return is: cbm2, cbm3, tab, oxbs. 

    ::bool dry_run:: Whether to return data for internal testing.  
    """
    def onlyAutosomal(df):
        df = df.loc[df.loc[:, "Chromosome"].str.match("^(chr)\d+$")]
        return df
    
    if dry_run:
        cbm2_path = "data/dryruns/cbm2/"
        cbm3_path = "data/dryruns/cbm3/" 

        oxbs_path = "data/dryruns/oxbs/"
        tab_path = "data/dryruns/tab/"

    else:
        cbm2_path = "data/modbases/nanopore/cbm2/"
        cbm3_path = "data/modbases/nanopore/cbm3/"

        oxbs_path = "data/modbases/public/CRR008808_oxBS/masked/"
        tab_path = "data/modbases/public/CRR008807_TAB/masked/"

    if split_biorep:
        bio_reps = ["Nanopore 1", "Nanopore 2"]
    else:
        bio_reps = ["Nanopore", "Nanopore"]

    with concurrent.futures.ThreadPoolExecutor(4) as ppe:
        all_futures = [ppe.submit(common.openReps, path, insert_cols={"Technique" : bio_rep}, **kwargs) for path, bio_rep in zip([cbm2_path, cbm3_path], bio_reps)]
        all_futures.append(ppe.submit(common.openReps, tab_path, insert_cols={"Technique" : "TAB"}, modbase="5hmC", **kwargs))
        all_futures.append(ppe.submit(common.openReps, oxbs_path, insert_cols={"Technique" : "oxBS"},  modbase="5mC", **kwargs))
        future_dfs = [onlyAutosomal(future.result()) for future in all_futures]

    return future_dfs

def annotate(df):
    feature_pr = annotation_features.fetch_feature_PyRange("/mnt/data1/doh28/analyses/mouse_hydroxymethylome_analysis/feature_references/RefSeq_Select/")

    annotated_df = pr.PyRanges(df).join(feature_pr, strandedness=False, suffix="_Feature", apply_strand_suffix=False).as_df()
    annotated_df = annotated_df.assign(readCount_vs_avg = lambda df: np.log2(df["readCount"]/df["readCount"].mean()))

    return annotated_df

def make_histogram(df, ax, color, **kwargs):
    return sns.histplot(df, x="readCount", ax=ax, color=color, **kwargs)  

@timer
def fig_main(dry_run):
    sns.set_style("ticks")
    mpl.rc('font', size=5)

    fig = plt.figure(figsize=(89/25.4, 60/25.4), dpi=600, layout="constrained")
    gs = GridSpec(2, 3, fig)

    cbm2, cbm3, tab, oxbs = fetch_data(dry_run, select_cols=["Chromosome",
                                                                    "Start", "End", 
                                                                    "readCount", 
                                                                    "Replicate", 
                                                                    "Technique"])
    cbm3 = cbm3.replace(["Rep. 1", "Rep. 2"], ["Rep. 3", "Rep. 4"])
    
    autosomal_dataframes = pd.concat([cbm2, cbm3, tab, oxbs])
    LEN_UNION = len(pr.PyRanges(autosomal_dataframes).merge(False, slack=-1)) # len of union of covered sites saved for later

    # Unique coverage comparison # 

    coverage_by_rep = autosomal_dataframes.loc[:, ("Replicate", "Technique")].value_counts(dropna=True).reset_index(name="Count")

    del autosomal_dataframes

    coverage_by_rep = coverage_by_rep.assign(Proportion_of_union = lambda r: (r["Count"]/LEN_UNION)*100)
    coverage_by_rep = coverage_by_rep.replace(["Nanopore 1", "Nanopore 2"], ["Nanopore", "Nanopore"])

    ax4 = fig.add_subplot(gs[1, 0])

    sns.barplot(coverage_by_rep, 
                x="Technique", y="Proportion_of_union", 
                hue="Technique", palette=sns.color_palette("BuGn_r", 4)[1:], 
                order=["oxBS", "TAB", "Nanopore"], 
                errorbar=("sd", 1), err_kws={
                    "linewidth" : 0.8, 
                    "solid_capstyle" : "butt"
                    },
                capsize=0.5,
                ax=ax4)
    
    ax4.set_ylabel("Percentage of all\ncovered sites")
    ax4.set_xlabel(None)
    ax4.tick_params("x", labelrotation=15)

    # Coverage depth histograms # 

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])

    with concurrent.futures.ThreadPoolExecutor(3) as hist_executor:
        hist_futures = [hist_executor.submit(make_histogram, df, ax, color, stat="percent", multiple="dodge", binrange=(0, 32), binwidth=2) for df, ax, color in zip(
            [oxbs, tab, pd.concat([cbm2, cbm3])],
            [ax1, ax2, ax3],
            [sns.color_palette("BuGn", 4)[1], sns.color_palette("BuGn", 4)[2], sns.color_palette("BuGn", 4)[3]])]
        for future in hist_futures:
            future.result() 
  
    ax1.set_title("oxBS-seq")
    ax2.set_title("TAB-seq")
    ax3.set_title("Nanopore")

    for index, ax in enumerate([ax1, ax2, ax3]):
        ax.set_title(string.ascii_lowercase[index], fontdict={"fontweight" : "bold"}, loc="left")
        ax.set_xlabel("Depth at CpG site")
        ax.set_ylabel("Percentage of sites")

    # Coverage depth by annotation # 

    ax5 = fig.add_subplot(gs[1, 1:])

    with concurrent.futures.ProcessPoolExecutor(4) as annotation_executor:
        annotated_futures = annotation_executor.map(annotate, [cbm2, cbm3, tab, oxbs])
        annotated_df = pd.concat([future for future in annotated_futures])

    annotated_df = annotated_df.replace(
        ["intergenic", "genes", "intron", "cds", "1kbPromoter"], 
        ["Intergenic", "Whole gene", "Intron", "CDS", "1kb Promoter"])

    annotated_df["feature_type"] = pd.Categorical(annotated_df["feature_type"], 
                                                  categories=["Intergenic", "Intron", "CDS", "1kb Promoter", "Whole gene"], 
                                                  ordered=True)

    sns.barplot(annotated_df, 
                x="feature_type", 
                y="readCount_vs_avg", 
                hue="Technique",
                palette=sns.color_palette("BuGn_r", 4)[1:], 
                errorbar=("sd", 1),
                err_kws={"linewidth" : 0.8,
                         "solid_capstyle" : "butt"},
                capsize=0.5,
                ax=ax5)

    ax5.axhline(0, linewidth=0.8, color="k", ls=":")
    ax5.set_ylabel(f"Depth difference\nto mean (Log$_2$)")
    ax5.set_xlabel(None)
    ax5.tick_params("x", labelrotation=15)

    sns.move_legend(ax5, "lower center", bbox_to_anchor=(.5, 1), ncol=3, title=None, frameon=False)

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