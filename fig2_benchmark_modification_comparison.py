import argparse
import pandas as pd
import gc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from AnalysisTools.helpers import timer
from AnalysisTools.common import *
import string
from AnalysisTools import annotation_features
from scipy import stats

def prep_dev_plot(df, bs_method, mod):
    df = (df.reset_index()
            .drop(columns=["Chromosome", "Start", "End"])
            .eval(f"diff = `{bs_method} mean` - `Nanopore mean`")
            .assign(Mod = mod))
    return df


def feature_stats(annotated_dataset: pd.DataFrame):
    feauture_stats = annotated_dataset.groupby("feature_type", observed=True).aggregate(
        median_5mC_z = pd.NamedAgg("zscore_5mC", np.median),
        median_5hmC_z = pd.NamedAgg("zscore_5hmC", np.median),
        mean_5mC = pd.NamedAgg("percentMeth_5mC", np.mean),
        mean_5hmC = pd.NamedAgg("percentMeth_5hmC", np.mean),
        median_5mC = pd.NamedAgg("percentMeth_5mC", np.median),
        median_5hmC = pd.NamedAgg("percentMeth_5hmC", np.median),
        var_5mC = pd.NamedAgg("percentMeth_5mC", np.var),
        var_5hmC = pd.NamedAgg("percentMeth_5hmC", np.var),
    )
    return feauture_stats.multiply({
        "median_5mC_z" : 1,
        "median_5hmC_z" : 1,
        "mean_5mC" : 100, 
        "mean_5hmC" : 100,
        "median_5mC" : 100,
        "median_5hmC" : 100,
        "var_5mC" : 100,
        "var_5hmC" : 100
    })

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
    ax4 = fig.add_subplot(gs[1, :])

    # Initial data collection # 

    print("Loading data")
    
    nano_path = "data/modbases/modbeds/nanopore_base5/"
    tab_path = "data/modbases/modbeds/tab/"
    oxbs_path = "data/modbases/modbeds/oxbs/"
    
    paths = [nano_path, tab_path, oxbs_path]
    cols = [["readCount", "N_mC", "N_hmC"], ["readCount", "N_mod"], ["readCount", "N_mod"]]

    with concurrent.futures.ProcessPoolExecutor(3) as fetch_executor:
        modbeds_list = [fetch_executor.submit(fetch_modbeds, path, cols, dryrun) for path, cols in zip(paths, cols)]
        modbeds_list = [future.result() for future in modbeds_list]

    print("Loaded data")

    [modbed.rename(columns={"N_mod" : "N_5hmC"}, inplace=True) for modbed in modbeds_list[1]]
    [modbed.rename(columns={"N_mod" : "N_5mC"}, inplace=True) for modbed in modbeds_list[2]]
    
    feature_annotator = annotation_features.Annotator("feature_references/feature_comparison/")
    gene_annotator = annotation_features.Annotator("feature_references/genic/")
    cgi_annotator = annotation_features.Annotator("feature_references/cgi/")
    
    nano_mb, *_ = modbeds_list
    annotators = [feature_annotator, gene_annotator, cgi_annotator]

    print("Annotating")

    def annotate_wrapper(annotator: annotation_features.Annotator):
        with concurrent.futures.ThreadPoolExecutor(4) as executor:
            annotated_all = executor.map(annotator.annotate, nano_mb)
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
        
        return rep_summary.reset_index()
    
    annotated_all = map(annotate_wrapper, annotators)
    annotated_all = pd.concat(annotated_all).replace(
        ["1kbPromoter", "allExons", "intron", "genes", "intergenic", "OpenSea"],
        ["Promoter", "Exon", "Intron", "Genic", "Intergenic", "Sea"]
    )

    print("Outputting annotated statistics to data/context_stats.tsv")
    feature_stats(annotated_all).to_csv("data/context_stats.tsv", sep="\t")
    
    annotated_all = annotated_all.melt(id_vars=["Start_Feature", "End_Feature", "feature_type"], 
                            value_vars=["zscore_5mC", "zscore_5hmC"],
                            value_name="zscore", var_name="Mod").replace(["zscore_5mC", "zscore_5hmC"], ["5mC", "5hmC"])

    tickorder = ["Promoter", "5UTR", "Intron", "Exon", "3UTR", "Genic", "Intergenic", "Sea", "Shelf", "Shore", "CGI"]
       
    annotated_all["feature_type"] = pd.Categorical(values=annotated_all["feature_type"], 
                                                   categories=tickorder, 
                                                   ordered=True)
    ax4.axhline(0, c="grey", lw=.8, ls=":")

    print("Writing violin plot of annotated positions")
    annotated_all.to_csv('source_data/fig2d_violin.csv.gz')
    sns.violinplot(annotated_all, 
                   x="feature_type", y="zscore", 
                   cut=2, linewidth=.5, density_norm="width",
                   hue="Mod", palette="GnBu", 
                   ax=ax4)    
    gc.collect()
    
    sns.move_legend(ax4, "lower left", title=None, bbox_to_anchor=(.4, 1), ncol=2, frameon=False)
    ax4.axvline(4.5, c="k", lw=.8)
    ax4.axvline(6.5, c="k", lw=.8)
    
    ax4.set_ylabel("Site modification (%) Z-Score")
    ax4.set_xlabel(None)
    ax4.set_ybound(-3, 3)

    with concurrent.futures.ProcessPoolExecutor(3) as merge_executor:
        merge_futures = [merge_executor.submit(merge_positions, modbed, True, col) 
                         for modbed, col in zip(modbeds_list, [["N_5mC", "N_5hmC"], ["N_5hmC"], ["N_5mC"]])]
        nanopore_average, tab_average, ox_average = [future.result().assign(Method = method) 
                                                     for future, method in zip(merge_futures, 
                                                                               ["Nanopore mean", "TAB mean", "oxBS mean"])]

    if dryrun:
        n=100000
    else: 
        n=1000000

    dfs = [nanopore_average, ox_average, nanopore_average, tab_average]
    csv_names = ['nanopore_site_average.csv.gz', 'ox_site_average.csv.gz', 'tab_site_average.csv.gz']
    csv_names = ['source_data/' + name for name in csv_names]

    for df, name in zip([nanopore_average, ox_average, tab_average], csv_names):
        df.to_csv(name)

    dfs = [df.sample(n, random_state=42) for df in dfs]

    xs = ["percentMeth_5mC", "percentMeth_5mC", "percentMeth_5hmC", "percentMeth_5hmC"]
    cs = [sns.color_palette("Greens_r", 4)[0], sns.color_palette("Greens_r", 4)[0], sns.color_palette("Blues_r", 4)[0], sns.color_palette("Blues_r", 4)[0]]
    lss = ["-", ":", "-", ":"]

    labels = [f"Nanopore", f"oxBS", f"Nanopore", f"TAB"]
    axs = [ax1, ax1, ax2, ax2]

    print("Plotting densities")

    with concurrent.futures.ThreadPoolExecutor(4) as kde_executor:
        kde_futures = [kde_executor.submit(sns.kdeplot, data=df, x=x, c=c, lw=.8, ls=ls, label=label, ax=ax) 
                        for df, x, c, ls, label, ax in zip(dfs, xs, cs, lss, labels, axs)]
        [future.result() for future in kde_futures]

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

    sns.move_legend(ax1, "upper left", frameon=False)
    sns.move_legend(ax2, "upper right", frameon=False)

    print("Plotting deviation")
    with concurrent.futures.ProcessPoolExecutor(2) as dev_plot_executor:
        dev_plot_futures = [dev_plot_executor.submit(prep_dev_plot, df, bs_method, mod) for df, bs_method, mod, in zip([nano_oxbs, nano_tab], 
                                                                                                                       ["oxBS", "TAB"], 
                                                                                                                       ["5mC", "5hmC"])]
        nano_oxbs, nano_tab = [future.result() for future in dev_plot_futures]
    
    gc.collect()
    [print("Sites compared: ", len(df)) for df in [nano_oxbs, nano_tab]]
    ax3.axvline(0, ls=":", c="grey", lw=0.8)

    concat = pd.concat([nano_oxbs, nano_tab])
    del nano_oxbs, nano_tab

    concat.to_csv('fig2c_diff_hist.csv.gz')

    sns.histplot(concat,
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

    ax3.set_ylim(0, 0.1)
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
                        default=[180, 80],
                        required=False,
                        help="Size of the figure produced in mm: (w, h).") 
    parser.add_argument("--fontsize", 
                        dest="fontsize", 
                        default=5,
                        required=False,
                        help="Size of figure font.") 

    args = parser.parse_args()

    fig_main(args.figsize, args.fontsize, args.dryrun)    