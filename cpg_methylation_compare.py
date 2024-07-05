import argparse
import pandas as pd
import gc
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from AnalysisTools.helpers import timer
from AnalysisTools.common import *
import string
from compare_features import *
from AnalysisTools import annotation_features
from scipy import stats

def prep_dev_plot(df, bs_method, mod):
    df = (df.reset_index()
            .drop(columns=["Chromosome", "Start", "End"])
            .eval(f"diff = `{bs_method} mean` - `Nanopore mean`")
            .assign(Mod = mod))
    return df

def annotate_wrapper(df_list: list[pd.DataFrame], annotator: annotation_features.Annotator):
    with concurrent.futures.ThreadPoolExecutor(4) as executor:
        annotated_all = executor.map(annotator.annotate, df_list)
        annotated_all = pd.concat([(annotation.drop(columns=["Strand"], errors="ignore").assign(Replicate = i+1))
                                   for i, annotation in enumerate(annotated_all)])
        
    grouped_all = annotated_all.groupby(["Chromosome", "Start", "End", "feature_type"], 
                                        observed=True)
    # now taking sum of all 
    rep_summary = grouped_all.sum(numeric_only=True).reset_index()

    for col in ["5mC", "5hmC"]:
        rep_summary.eval(f"percentMeth_{col} = (N_{col}/readCount)", inplace=True)
        rep_summary[f"percentMeth_{col}"] = rep_summary[f"percentMeth_{col}"].where(rep_summary[f"percentMeth_{col}"] <= 1, 1)        
        rep_summary[f"asin_{col}"] = np.arcsin(rep_summary[f"percentMeth_{col}"])
        # zscore is calculated element-wise           
        rep_summary[f"zscore_{col}"] = stats.zscore(rep_summary[f"asin_{col}"])
        
    return rep_summary

@timer
def fig_main(figsize, fontsize, dryrun=True):
    fig = plt.figure(figsize=[a for a in map(lambda d: float(d)/25.4, figsize)], 
                     dpi=600, layout="constrained")

    sns.set_style("ticks")
    mpl.rc('font', size=fontsize)

    gs = GridSpec(2, 3, fig, height_ratios=(1, 1.3))

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])
    sgs = gs[1, :].subgridspec(1, 5)
    ax4 = fig.add_subplot(sgs[0, :2])
    ax5 = fig.add_subplot(sgs[0, 2])
    ax6 = fig.add_subplot(sgs[0, 3:])

    # Initial data collection # 

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
        modbeds_list = [fetch_executor.submit(fetch_modbeds, path, cols) for path, cols in zip(paths, cols)]
        modbeds_list = [future.result() for future in modbeds_list]

    print("Loaded data")

    [modbed.rename(columns={"N_mod" : "N_5hmC"}, inplace=True) for modbed in modbeds_list[1]]
    [modbed.rename(columns={"N_mod" : "N_5mC"}, inplace=True) for modbed in modbeds_list[2]]
    
    feature_annotator = annotation_features.Annotator("/mnt/data1/doh28/analyses/mouse_hydroxymethylome_analysis/feature_references/feature_comparison/")
    gene_annotator = annotation_features.Annotator("/mnt/data1/doh28/analyses/mouse_hydroxymethylome_analysis/feature_references/genic/")
    cgi_annotator = annotation_features.Annotator("/mnt/data1/doh28/analyses/mouse_hydroxymethylome_analysis/feature_references/cgi/")

    nanopore_annotated = (annotate_wrapper(modbeds_list[0], feature_annotator)
                          .reset_index()
                          .replace(["1kbPromoter", "allExons", "intron"],
                                   ["Promoter", "Exon", "Intron"])
                          .melt(id_vars=["Start_Feature", "End_Feature", "feature_type"], 
                                value_vars=["zscore_5mC", "zscore_5hmC"],
                                value_name="zscore", var_name="Mod")
                                .replace(["zscore_5mC", "zscore_5hmC"],
                                         ["5mC", "5hmC"]))

    tickorder = ["Promoter", "5UTR", "Intron", "Exon", "3UTR"]

    [ax.axhline(0, ls=":", c="grey", label="Genomic mean") for ax in [ax4, ax5, ax6]]
    
    nanopore_annotated["feature_type"] = pd.Categorical(nanopore_annotated["feature_type"],
                                                        categories=tickorder,
                                                        ordered=True)

    sns.barplot(nanopore_annotated, 
                x="feature_type", y="zscore",
                estimator="mean", errorbar=None,
                hue="Mod", palette="GnBu",
                hue_order=["5mC", "5hmC"],
                ax=ax4)

    sns.move_legend(ax4, loc="lower right", 
                    frameon=False, title=None, 
                    )
    # ax4.get_legend().set_in_layout(False)
                    
    nanopore_annotated = (annotate_wrapper(modbeds_list[0], gene_annotator)
                        .reset_index()
                        .replace(["genes", "intergenic"], ["Genic", "Intergenic"])
                        .melt(id_vars=["Start_Feature", "End_Feature", "feature_type"],
                            value_vars=["zscore_5mC", "zscore_5hmC"],
                            value_name="zscore", var_name="Mod")
                            .replace(["zscore_5mC", "zscore_5hmC"], 
                                     ["5mC", "5hmC"]))
    
    tickorder = ["Genic", "Intergenic"]

    nanopore_annotated["feature_type"] = pd.Categorical(nanopore_annotated["feature_type"],
                                                        categories=tickorder,
                                                        ordered=True)
    
    sns.barplot(nanopore_annotated, 
                x="feature_type", y="zscore",
                estimator="mean", errorbar=None,
                legend=None,
                hue="Mod", palette="GnBu",
                hue_order=["5mC", "5hmC"],
                ax=ax5)

    nanopore_annotated = (annotate_wrapper(modbeds_list[0], cgi_annotator)
                          .reset_index()
                          .replace(["OpenSea"], ["Sea"])
                          .melt(id_vars=["Start_Feature", "End_Feature", "feature_type"],
                                value_vars=["zscore_5mC", "zscore_5hmC"],
                                value_name="zscore", var_name="Mod")
                                .replace(["zscore_5mC", "zscore_5hmC"], 
                                         ["5mC", "5hmC"]))
    
    tickorder = ["Sea", "Shelf", "Shore", "CGI"]

    nanopore_annotated["feature_type"] = pd.Categorical(nanopore_annotated["feature_type"],
                                                        categories=tickorder,
                                                        ordered=True)
    
    sns.lineplot(nanopore_annotated, 
                x="feature_type", y="zscore",
                hue="Mod", style="Mod", palette="GnBu",
                hue_order=["5mC", "5hmC"],
                ax=ax6)
    
    sns.move_legend(ax6, "lower left", frameon=False, title=None)
    [ax.set_xlabel(None) for ax in [ax4, ax5, ax6]]
    [ax.set_ylabel(None) for ax in [ax5, ax6]]

    ax4.set_ylabel("Site modification (%) Z-Score")

    with concurrent.futures.ProcessPoolExecutor(3) as merge_executor:
        merge_futures = [merge_executor.submit(merge_positions, modbed, True, col) 
                         for modbed, col in zip(modbeds_list, [["N_5mC", "N_5hmC"], ["N_5hmC"], ["N_5mC"]])]
        nanopore_average, tab_average, ox_average = [future.result().assign(Method = method) 
                                                     for future, method in zip(merge_futures, 
                                                                               ["Nanopore mean", "TAB mean", "oxBS mean"])]

    dfs = [nanopore_average, ox_average, nanopore_average, tab_average]
    xs = ["percentMeth_5mC", "percentMeth_5mC", "percentMeth_5hmC", "percentMeth_5hmC"]
    cs = [sns.color_palette("Greens_r", 4)[0], sns.color_palette("Greens_r", 4)[0], sns.color_palette("Blues_r", 4)[0], sns.color_palette("Blues_r", 4)[0]]
    lss = ["-", ":", "-", ":"]

    labels = [f"Nanopore,\nn={len(nanopore_average):.2e}", f"oxBS,\nn={len(ox_average):.2e}", f"Nanopore,\nn={len(nanopore_average):.2e}", f"TAB,\nn={len(tab_average):.2e}"]
    axs = [ax1, ax1, ax2, ax2]

    print("Plotting ECDFs")

    with concurrent.futures.ThreadPoolExecutor(4) as ecdf_executor:
        ecdf_futures = [ecdf_executor.submit(sns.ecdfplot, data=df, x=x, c=c, lw=.8, ls=ls, label=label, ax=ax) 
                        for df, x, c, ls, label, ax in zip(dfs, xs, cs, lss, labels, axs)]
        [future.result() for future in ecdf_futures]

    ax1.set_xlabel("Site 5mC (%)")
    ax2.set_xlabel("Site 5hmC (%)")

    ax1.set_title("5mC", loc="center")
    ax2.set_title("5hmC", loc="center")

    nano_oxbs = (pd.concat([nanopore_average, ox_average], join="inner", copy=False)
                 .pivot_table(values="percentMeth_5mC", index=["Chromosome", "Start", "End"], columns="Method").dropna())
    del ox_average

    nano_tab = (pd.concat([nanopore_average, tab_average], join="inner", copy=False)
                .pivot_table(values="percentMeth_5hmC", index=["Chromosome", "Start", "End"], columns="Method").dropna())
    del tab_average, nanopore_average

    gc.collect()

    for ax in [ax1, ax2]:
        ax.legend() 
        ax.set_ylim(0, 1)
        ax.set_xlim(0, 100)
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
    
    ax3.axvline(0, ls=":", c="grey", lw=0.8)
    
    sns.histplot(pd.concat([nano_oxbs, nano_tab]),
                x="diff", 
                stat="proportion",
                binrange=(-50, 50), binwidth=2,
                element="step", fill=False, 
                legend=False,
                lw=0.8,
                hue="Mod", palette="GnBu",
                ax=ax3)
        
    ax3.lines[1].set_linestyle("--")
    legend = [Line2D([0], [0], lw=.8, ls="-", c=sns.color_palette("GnBu", 2)[0]),
              Line2D([0], [0], lw=.8, ls="--", c=sns.color_palette("GnBu", 2)[1])]
    ax3.legend(legend, ["5mC", "5hmC"])
    sns.move_legend(ax3, "upper left", frameon=False, title=None)

    ax3.set_ylim((0, 0.1))
    ax3.set_xlim((-50, 50))
    ax3.set_xlabel("Bisulphite % - Nanopore %")

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
                        default=False,
                        required=False,
                        help="Whether a test output is produced.") 
    parser.add_argument("--figsize", 
                        dest="figsize", 
                        nargs=2,
                        default=[120, 80],
                        required=False,
                        help="Size of the figure produced in mm: (w, h).") 
    parser.add_argument("--fontsize", 
                        dest="fontsize", 
                        default=5,
                        required=False,
                        help="Size of figure font.") 

    args = parser.parse_args()

    fig_main(args.figsize, args.fontsize, args.dryrun)    