"""
Will read the table files to compare genomic coverage depth of multiple datasets. 

Takes a while to run as it needs to generate new intermediate data files for each dataset. 
"""

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
from AnalysisTools import annotation_features
from AnalysisTools import common
import numpy as np
import pandas as pd
import string
import pyranges as pr
import argparse

def fig_main(dry_run):
    cbm2_auto, cbm3_auto, tab_auto, oxbs_auto = common.fetch_data_Parallel(dry_run,
                                                                           select_cols=["Chromosome", "Start", "End", "readCount", "Replicate"])
    all_techs = [df.assign(Technique = f"{tech}") for df, tech in zip([cbm2_auto, cbm3_auto, tab_auto, oxbs_auto], 
                                                                     ["Nanopore", "Nanopore", "TAB", "oxBS"])]

    # calculate all covered CpG sites in all replicates
    union = len(pr.PyRanges(pd.concat(all_techs)).merge(slack=-1))

    rep_ls = []
    for tech in all_techs:
        for replicate in tech["Replicate"].unique():
            rep_df = len(tech.query(f"Replicate == '{replicate}'"))
            stat_df = pd.DataFrame({"Technique" : [tech.Technique.values.all()], 
                                    "Replicate" : [replicate], 
                                    "N_positions" : [rep_df], 
                                    "Proportion_of_Union" : [(rep_df/union)*100]})
            rep_ls.append(stat_df)

    feature_pr = annotation_features.featureRefPyRange("/mnt/data1/doh28/analyses/mouse_hydroxymethylome_analysis/feature_references/RefSeq_Select/")
    
    def annotate(df):
        annotated_df = pr.PyRanges(df).join(feature_pr, strandedness=False, suffix="_Feature", apply_strand_suffix=False).as_df()
        annotated_df = annotated_df.assign(readCount_vs_avg = lambda df: np.log2(df["readCount"]/df["readCount"].mean()))
        return annotated_df

    annotated_ls = [annotate(tech) for tech in all_techs]

    feature_depth = pd.concat(annotated_ls).replace(
        ["intergenic", "genes", "intron", "cds", "1kbPromoter"], 
        ["Intergenic", "Whole gene", "Intron", "CDS", "1kb Promoter"]
    )
    feature_depth["feature_type"] = pd.Categorical(feature_depth["feature_type"], 
                                                categories=["Intergenic", "Intron", "CDS", "1kb Promoter", "Whole gene"], 
                                                ordered=True)
    
    # Fig # 

    fig = plt.figure(figsize=(89/25.4, 60/25.4), dpi=600, layout="constrained")

    sns.set_style("ticks")
    mpl.rc('font', size=5)

    gs = GridSpec(2, 3, fig)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])

    # OxBS coverage histogram # 

    sns.histplot(oxbs_auto,
                x="readCount", 
                stat="percent",
                binrange=(0, 32),
                binwidth=2, 
                color=sns.color_palette("BuGn", 4)[1],
                multiple="dodge",
                ax=ax1)
    
    ax1.set_title("oxBS-seq")

    sns.histplot(tab_auto,
                x="readCount", 
                stat="percent",
                binrange=(0, 32),
                binwidth=2, 
                color=sns.color_palette("BuGn", 4)[2],
                multiple="dodge",
                ax=ax2)
    
    ax2.set_title("TAB-seq")
    
    sns.histplot(pd.concat([cbm2_auto, cbm3_auto]),
                x="readCount", 
                stat="percent",
                binrange=(0, 32),
                binwidth=2, 
                color=sns.color_palette("BuGn", 4)[3],
                multiple="dodge",
                ax=ax3)

    ax3.set_title("Nanopore")

    # Formatting # 

    for index, ax in enumerate([ax1, ax2, ax3]):
        #sns.move_legend(ax, "upper right")

        ax.set_title(string.ascii_lowercase[index], fontdict={"fontweight" : "bold"}, loc="left")
        ax.set_xlabel("Depth at CpG site")
        ax.set_ylabel("Percentage of sites")

    # Percentage of all covered sites # 

    ax4 = fig.add_subplot(gs[1, 0])

    sns.barplot(pd.concat(rep_ls), 
                x="Technique", y="Proportion_of_Union", 
                hue="Technique",
                order=["oxBS", "TAB", "Nanopore"],
                palette=sns.color_palette("BuGn_r", 4)[1:],  
                errorbar=("sd", 1), 
                err_kws={"linewidth" : 0.8,
                         "solid_capstyle" : "butt"},
                capsize=0.5,
                ax=ax4)

    ax4.set_ylabel("Percentage of all\ncovered sites")
    ax4.set_xlabel(None)
    ax4.tick_params("x", labelrotation=15)

    ax5 = fig.add_subplot(gs[1, 1:])

    sns.barplot(feature_depth, 
                x="feature_type", 
                y="readCount_vs_avg", 
                hue="Technique",
                # order=["oxBS", "TAB", "Nanopore"],
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
    return fig.savefig("plots/compare_coverage.png")

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