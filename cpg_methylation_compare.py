import argparse
import pandas as pd
from sklearn.metrics import RocCurveDisplay, roc_auc_score
from sklearn.preprocessing import Binarizer
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from AnalysisTools.helpers import timer
from AnalysisTools.common import *
import string

@timer
def fig_main(dryrun=True):
  
    fig = plt.figure(figsize=(89/25.4, 120/25.4), dpi=600, layout="constrained")

    sns.set_style("ticks")
    mpl.rc('font', size=5)

    gs = GridSpec(3, 2, fig)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    ax5 = fig.add_subplot(gs[2, 0])

    # Initial data collection # 

    nanopore_average = merge_positions(fetch_nanopore(["readCount", "N_mC", "N_hmC"], dryrun=dryrun), 
                                       ["N_5mC", "N_5hmC"])
    nanopore_average["Method"] = "Nanopore mean"

    ox_average = merge_positions(fetch_oxbs(["readCount", "N_mod"]), 
                                 "N_5mC")
    ox_average["Method"] = "oxBS mean"

    tab_average = merge_positions(fetch_tab(["readCount", "N_mod"]), 
                                  "N_5hmC")
    tab_average["Method"] = "TAB mean"

    # General distribution comparison # 

    sns.ecdfplot(nanopore_average,
             x="percentMeth_5mC",
             c=sns.color_palette("BuGn_r", 4)[0],
             lw=0.8, ls=":",
             label="Nanopore 5mC",
             ax=ax1)

    sns.ecdfplot(ox_average,
                x="percentMeth_5mC",
                c=sns.color_palette("BuGn_r", 4)[0],
                lw=0.8, 
                label="oxBS 5mC",
                ax=ax1)

    ax1.set_xlabel("Site 5mC (%)")
    
    nano_oxbs = pd.concat([nanopore_average, ox_average], join="inner").pivot_table(values="percentMeth_5mC", index=["Chromosome", "Start", "End"], columns="Method").dropna()
    del ox_average

    sns.ecdfplot(nanopore_average,
                x="percentMeth_5hmC",
                c=sns.color_palette("BuGn_r", 4)[1],
                lw=0.8, ls=":",
                label="Nanopore 5hmC",
                ax=ax2)
    
    
    sns.ecdfplot(tab_average,
                x="percentMeth_5hmC",
                c=sns.color_palette("BuGn_r", 4)[1],
                lw=0.8, 
                label="TAB 5hmC",
                ax=ax2)
    
    ax2.set_xlabel("Site 5hmC (%)")

    nano_tab = pd.concat([nanopore_average, tab_average], join="inner").pivot_table(values="percentMeth_5hmC", index=["Chromosome", "Start", "End"], columns="Method").dropna()
    del tab_average, nanopore_average

    # General distribution comparison # 

    pos_ctrl, neg_ctrl = fetch_controls(["readCount", "N_mC"], dryrun=True)
    m_average = merge_positions(pos_ctrl, "N_5mC")
    u_average = merge_positions(neg_ctrl, "N_5mC")

    del pos_ctrl, neg_ctrl

    m_average["Truth"], u_average["Truth"] = 1, 0

    for ax in [ax1, ax2]:
        ax.legend() 
        ax.set_aspect(100)
        ax.set_ylabel("Cumulative proportion")

    sns.move_legend(ax1, "upper left", frameon=False)
    sns.move_legend(ax2, "lower right", frameon=False)

    # KDE site comparison # 

    sns.kdeplot(nano_oxbs,
            x="Nanopore mean", y="oxBS mean", 
            fill=True, color=sns.color_palette("BuGn", 5)[4],
            ax=ax3)

    ax3.set_xlabel("oxBS-seq site 5mC (%)")
    ax3.set_ylabel("Nanopore site 5mC (%)")

    sns.kdeplot(nano_tab,
                x="Nanopore mean", y="TAB mean", 
                fill=True, color=sns.color_palette("BuGn", 5)[4],
                ax=ax4)

    ax4.set_xlabel("TAB-seq site 5hmC (%)")
    ax4.set_ylabel("Nanopore site 5hmC (%)")

    sns.despine()

    for ax in [ax3, ax4]:
        ax.set_aspect("equal")

        ax.set_xticks(range(0, 120, 20))
        ax.set_yticks(range(0, 120, 20))
        sns.despine(ax=ax, top=False, right=False)
        
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)
                        
    # site comparison ROC # 

    binarizer = Binarizer(threshold=(2/3)*100)

    nano_oxbs["Truth"] = binarizer.fit_transform(nano_oxbs["oxBS mean"].to_numpy().reshape((-1, 1)))
    nano_tab["Truth"] = binarizer.fit_transform(nano_tab["TAB mean"].to_numpy().reshape((-1, 1)))

    control_array = pd.concat([m_average, u_average], ignore_index=True)

    RocCurveDisplay.from_predictions(control_array["Truth"], control_array["percentMeth_5mC"],
                                    c=sns.color_palette("BuGn_r", 4)[0], lw=0.8,
                                    label=f"5mC (Controls): AUC {round(roc_auc_score(control_array['Truth'], control_array['percentMeth_5mC']), 3)}",
                                    ax=ax5)
    del control_array

    RocCurveDisplay.from_predictions(nano_oxbs["Truth"], nano_oxbs["Nanopore mean"],
                                    c=sns.color_palette("BuGn_r", 4)[1], lw=0.8,
                                    label=f"5mC (oxBS): AUC {round(roc_auc_score(nano_oxbs['Truth'], nano_oxbs['Nanopore mean']), 3)}",
                                    ax=ax5)
    del nano_oxbs

    RocCurveDisplay.from_predictions(nano_tab["Truth"], nano_tab["Nanopore mean"],
                                    c=sns.color_palette("BuGn_r", 4)[2], lw=0.8,
                                    label=f"5hmC (TAB): AUC {round(roc_auc_score(nano_tab['Truth'], nano_tab['Nanopore mean']), 3)}",
                                    ax=ax5)
    
    ax5.set_aspect("equal")
    sns.move_legend(ax5, "lower right", frameon=False, bbox_to_anchor=(1.25, 0))
    ax5.get_legend().set_in_layout(False)

    for index, ax in enumerate([ax1, ax2, ax3, ax4, ax5]):
        ax.set_title(f"{string.ascii_lowercase[index]}", fontdict={"fontweight" : "bold"}, loc="left")

    if dryrun:
        outpath = "plots/tests/cpg_methylation_compare.png"
    else:
        outpath = "plots/cpg_methylation_compare.png"
    return fig.savefig(outpath)

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