import argparse
import pandas as pd
import gc
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from AnalysisTools.helpers import timer
from AnalysisTools.common import *
import string
from compare_features import *
from AnalysisTools import annotation_features

def prep_dev_plot(df, bs_method, mod):
    df = (df.reset_index()
            .drop(columns=["Chromosome", "Start", "End"])
            .eval(f"diff = `{bs_method} mean` - `Nanopore mean`")
            .assign(Mod = mod))
    return df

def annotate_wrapper(df_list: list[pd.DataFrame], annotator: annotation_features.Annotator):
    mean_dict = {}
    for col in ["5mC", "5hmC"]:
        try:
            readcount_tot = np.sum([df["readCount"].sum() for df in df_list])
            total = np.sum([df[f"N_{col}"].sum() for df in df_list])
            mean_value = (total/readcount_tot)*100

            mean_dict.update({col : mean_value})
        except:
            # print(col, "not found in", df_list[0].columns)
            pass

    def group_feature(annotated_df: pd.DataFrame):
        grouped = annotated_df.groupby(["Start_Feature", "End_Feature", "feature_type"], 
                                       observed=True, as_index=True)
        
        features = grouped.sum(numeric_only=True)
        features["CpG_count"] = grouped["Start"].count()
        features = features.loc[features.eval("CpG_count > 9")]

        for col in ["5mC", "5hmC"]:
            try:
                features.eval(f"percentMeth_{col} = (N_{col}/readCount)*100", inplace=True)
                features[f"enrichment_{col}"] = np.log2((features[f"percentMeth_{col}"]+1)/(mean_dict[col]+1))
            except:
                pass    
        return features

    with concurrent.futures.ThreadPoolExecutor(len(df_list)) as executor:
        annotated_all = executor.map(annotator.annotate, df_list)
        annotated_all = [(annotation.drop(columns=["Strand"], errors="ignore")) 
                         for annotation in annotated_all]
            
    with concurrent.futures.ThreadPoolExecutor(len(df_list)) as executor:
        grouped_all = pd.concat([df.assign(Replicate = i) 
                                 for i, df in enumerate(executor.map(group_feature, annotated_all))])
    
    grouped_all = grouped_all.groupby(["Replicate", "feature_type"])

    rep_summary = grouped_all.mean(numeric_only=True)

    return rep_summary

@timer
def fig_main(figsize, fontsize, dryrun=True):
    fig = plt.figure(figsize=[a for a in map(lambda d: float(d)/25.4, figsize)], 
                     dpi=600, layout="constrained")

    sns.set_style("ticks")
    mpl.rc('font', size=fontsize)

    gs = GridSpec(2, 3, fig)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])
    ax4 = fig.add_subplot(gs[1, :2])
    ax5 = fig.add_subplot(gs[1, 2])

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
        modbeds_list = [fetch_executor.submit(fetch_modbeds, path, cols) for path, cols in zip(paths, cols)]
        modbeds_list = [future.result() for future in modbeds_list]

    print("Loaded data")

    [modbed.rename(columns={"N_mod" : "N_5hmC"}, inplace=True) for modbed in modbeds_list[1]]
    [modbed.rename(columns={"N_mod" : "N_5mC"}, inplace=True) for modbed in modbeds_list[2]]
    
    feature_annotator = annotation_features.Annotator("/mnt/data1/doh28/analyses/mouse_hydroxymethylome_analysis/feature_references/feature_comparison/")
    cgi_annotator = annotation_features.Annotator("/mnt/data1/doh28/analyses/mouse_hydroxymethylome_analysis/feature_references/cgi/")

    # feature annotation 

    with concurrent.futures.ProcessPoolExecutor(3) as annotation_executor:
        feature_annotations = [annotation_executor.submit(annotate_wrapper, df_list, feature_annotator) 
                               for df_list in modbeds_list]
        feature_annotations = [future.result().assign(Method = method)
                               for future, method in zip(feature_annotations, ["Nanopore", "TAB", "oxBS"])]
            
    nanopore_tab = pd.concat(feature_annotations[:2]).reset_index().replace(
        ["1kbPromoter", "allExons", "intron", "genes"], 
        ["Promoter", "Exon", "Intron", "Genic"]
    )

    # nanopore_oxbs = pd.concat([feature_annotations[0], feature_annotations[2]])
    tickorder = ["Promoter", "5UTR", "Intron", "Exon", "3UTR", "Genic"]
    
    other_fig = plt.figure(figsize=(16.73/2.54, 6.47/2.54), layout="constrained")
    gs = other_fig.add_gridspec(1, 3)
    ax4 = other_fig.add_subplot(gs[0, :2])
    ax5 = other_fig.add_subplot(gs[0, 2])

    [ax.axhline(0, ls="--", c="grey", label="Genomic mean") for ax in [ax4, ax5]]
    
    sns.barplot(nanopore_tab, 
                x="feature_type", y="enrichment_5hmC",
                hue="Method", palette="Blues",
                errorbar=("sd", 1),
                err_kws={'linewidth': 1}, capsize=.5,
                order=tickorder,
                ax=ax4)
    
    feature_stats = nanopore_tab.groupby(["feature_type", "Method"])
    for feature in tickorder:
        nanopore = np.array(feature_stats.get_group((feature, "Nanopore"))["enrichment_5hmC"])
        tab = np.array(feature_stats.get_group((feature, "TAB"))["enrichment_5hmC"])

        print(feature, stats.ttest_ind(nanopore, tab, equal_var=False))

    sns.move_legend(ax4, loc="lower right", frameon=False, title=None)

    with concurrent.futures.ProcessPoolExecutor(3) as annotation_executor:
        cgi_annotations = [annotation_executor.submit(annotate_wrapper, df_list, cgi_annotator) 
                               for df_list in modbeds_list]
        cgi_annotations = [future.result().assign(Method = method)
                           for future, method in zip(cgi_annotations, ["Nanopore", "TAB", "oxBS"])]
    
    nanopore_tab = pd.concat(cgi_annotations[:2]).reset_index()
    # nanopore_oxbs = pd.concat([cgi_annotations[0], cgi_annotations[2]])
    tickorder = ["Shelf", "Shore", "CGI"]

    sns.barplot(nanopore_tab, 
                x="feature_type", y="enrichment_5hmC",
                hue="Method", palette="Blues",
                err_kws={'linewidth': 1}, capsize=.5,
                errorbar=("sd", 1),
                order=tickorder,
                ax=ax5)
    
    [ax.set_xlabel(None) for ax in [ax4, ax5]]
    
    feature_stats = nanopore_tab.groupby(["feature_type", "Method"])
    for feature in tickorder:
        nanopore = np.array(feature_stats.get_group((feature, "Nanopore"))["enrichment_5hmC"])
        tab = np.array(feature_stats.get_group((feature, "TAB"))["enrichment_5hmC"])

        print(feature, stats.ttest_ind(nanopore, tab, equal_var=False))

    sns.move_legend(ax5, loc="lower left", frameon=False, title=None)

    ax4.set_ylabel(f"Log$_{2}$ FC vs. genomic mean")
    ax4.set_ylim(-3.5, .5)
    ax5.set_ylabel(f"Log$_{2}$ FC vs. genomic mean")
    ax5.set_ylim(-3, .5)

    sns.despine(other_fig)
    other_fig.savefig("plots/presentation_features.svg", dpi=600)
    exit()

    with concurrent.futures.ProcessPoolExecutor(3) as merge_executor:
        merge_futures = [merge_executor.submit(merge_positions, modbed, True, col) 
                         for modbed, col in zip(modbeds_list, [["N_5mC", "N_5hmC"], ["N_5hmC"], ["N_5mC"]])]
        nanopore_average, tab_average, ox_average = [future.result().assign(Method = method) 
                                                     for future, method in zip(merge_futures, 
                                                                               ["Nanopore mean", "TAB mean", "oxBS mean"])]
     
    # General distribution comparison # 
        
    # other_fig = plt.figure(figsize=(14.97/2.54, 12.37/2.54), layout="constrained", dpi=600)
    # gs = other_fig.add_gridspec(2, 6)
    # ax1 = other_fig.add_subplot(gs[0, 1:3])
    # ax2 = other_fig.add_subplot(gs[0, 3:5])
    # ax3 = other_fig.add_subplot(gs[1, 2:4])

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

    ax1.set_title(f"5mC", loc="center")
    ax2.set_title(f"5hmC", loc="center")
    
    # for average in [nanopore_average, ox_average]:
    #     len(average.query("percentMeth_5mC > 0"))/len(average)

    nano_oxbs = pd.concat([nanopore_average, ox_average], join="inner", copy=False).pivot_table(values="percentMeth_5mC", index=["Chromosome", "Start", "End"], columns="Method").dropna()
    del ox_average
    
    # for average in [nanopore_average, tab_average]:
    #     len(average.query("percentMeth_5hmC > 0"))/len(average)

    nano_tab = pd.concat([nanopore_average, tab_average], join="inner", copy=False).pivot_table(values="percentMeth_5hmC", index=["Chromosome", "Start", "End"], columns="Method").dropna()
    del tab_average, nanopore_average

    gc.collect()

    for ax in [ax1, ax2]:
        ax.legend() 
        ax.set_ylim(0, 1)
        ax.set_xlim(0, 100)
        ax.set_ylabel("Cumulative proportion")
        # ax.get_legend().set_in_layout(False)

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
    ax3.set_xlabel("Bisulphite % - Nanopore %")

    sns.move_legend(ax3, "upper right", frameon=False, title=None)
    sns.despine()

    # other_fig.savefig("plots/presentation_cpg_compare.svg")
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
                        default=[89, 89],
                        required=False,
                        help="Size of the figure produced in mm: (w, h).") 
    parser.add_argument("--fontsize", 
                        dest="fontsize", 
                        default=5,
                        required=False,
                        help="Size of figure font.") 

    args = parser.parse_args()

    fig_main(args.figsize, args.fontsize, args.dryrun)    