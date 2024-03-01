import argparse
import pandas as pd
from AnalysisTools.common import *
from AnalysisTools import annotation_features
from AnalysisTools.helpers import timer
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
import string
from scipy import stats

def group_feature(annotated_df):
    grouped = annotated_df.groupby(["Start_Feature", "End_Feature", "feature_type", "Method"], observed=True, as_index=False).agg({
        "readCount" : np.sum,
        "N_5hmC" : np.sum,
        "Start" : "count"
    }).rename(columns={"Start" : "CpG_count"})
    return grouped

def assign_enrichment(grouped_df, mean):
    out_df = (grouped_df.assign(avg_5hmC = lambda r: (r["N_5hmC"]/r["readCount"])*100)
                        .assign(enrichment = lambda r: np.log2((r["avg_5hmC"]+1)/(mean+1))))
    return out_df

def enrichment_wrap(annotated_df, mean):
    grouped_df = group_feature(annotated_df)
    return assign_enrichment(grouped_df, mean)

def annotate_main(df_list, annotator):
    
    global means
    means = [df["percentMeth_5hmC"].mean() for df in df_list]

    with concurrent.futures.ProcessPoolExecutor(len(df_list)) as executor:
        annotation_futures = [executor.submit(annotator.annotate, df) for df in df_list]
        annotated_all = [future.result().drop(columns=["Strand"], errors="ignore") for future in annotation_futures]
        enrichment_dfs = [executor.submit(enrichment_wrap, df, mean).result() for df, mean in zip(annotated_all, means)]
        enrichment_dfs = [df.drop(columns=["readCount", "N_5hmC", "Method"]) for df in enrichment_dfs]

        enrichment_all = pd.merge(*enrichment_dfs, 
                          on=["Start_Feature", "End_Feature", "feature_type"],
                          suffixes=["_Nanopore", "_TAB"])
    return enrichment_all

def plot_kde(df, ax):
    df = df.query("CpG_count_Nanopore > 4 & CpG_count_TAB > 4")
    
    return sns.kdeplot(df, x="enrichment_Nanopore", y="enrichment_TAB", 
                       color=sns.color_palette("Greens", 5)[4], fill=True,
                       ax=ax)

@timer
def fig_main(dryrun):
    # Loading data # 

    nanopore_dfs = fetch_nanopore(["readCount", "N_hmC"], dryrun=dryrun)
    nanopore_df = merge_positions(nanopore_dfs, ["N_5hmC"], False).reset_index()
    nanopore_df["Method"] = "Nanopore"

    tab_dfs = fetch_tab(["readCount", "N_mod"], dryrun=dryrun)
    tab_df = merge_positions(tab_dfs, ["N_5hmC"], False).reset_index()
    tab_df["Method"] = "TAB"

    # Data processing # 

    all_dfs = [nanopore_df, tab_df]
    del nanopore_df, tab_df

    feature_annotator = annotation_features.Annotator("/mnt/data1/doh28/analyses/mouse_hydroxymethylome_analysis/feature_references/feature_comparison/")
    cgi_annotator = annotation_features.Annotator("/mnt/data1/doh28/analyses/mouse_hydroxymethylome_analysis/feature_references/cgi/")

    features_enrichment = annotate_main(all_dfs, feature_annotator)
    cgi_enrichment = annotate_main(all_dfs, cgi_annotator)

    del all_dfs

    features_enrichment = features_enrichment.replace(["intron", "cds", "intergenic", "allExons", "genes", "1kbPromoter"],
                                                  ["Intron", "Exon (coding)", "Intergenic", "Exon (all)", "Whole gene", "1kb promoter"])

    features_enrichment["feature_type"] = pd.Categorical(features_enrichment.feature_type, 
                                                        ["Intergenic", "1kb promoter", "5UTR", "Intron", "Exon (coding)", "Exon (all)", "Whole gene", "3UTR"], 
                                                        ordered=True)
    
    def prep_barplot(annotated_df, categories):
        barplot_df = (annotated_df.query("CpG_count_Nanopore > 4 & CpG_count_TAB > 4")
                                    .melt(id_vars=["Start_Feature", "End_Feature", "feature_type"], 
                                          value_vars=["avg_5hmC_Nanopore", "avg_5hmC_TAB"])
                                    .replace(["avg_5hmC_Nanopore", "avg_5hmC_TAB"], ["Nanopore", "TAB"]))
        
        barplot_df["feature_type"] = pd.Categorical(barplot_df.feature_type, categories, ordered=True)
        return barplot_df
    
    barplot_feature = (prep_barplot(features_enrichment, ["Intergenic", "1kb promoter", "5UTR", "Intron", "Exon (coding)", "3UTR"])
                     .replace(["Whole gene", "Exon (all)"], [None, None]))

    barplot_cgi = prep_barplot(cgi_enrichment, ["Shelf", "Shore", "CGI"])
    
    # Fig setup # 

    fig = plt.figure(figsize=(120/25.4, 100/25.4), 
                     dpi=600, 
                     layout="constrained")
    gs = GridSpec(3, 6, fig)

    sns.set_style("ticks")
    mpl.rc('font', size=5)

    ax0 = fig.add_subplot(gs[0, :4])

    sns.barplot(barplot_feature, x="feature_type", y="value", 
                hue="variable", 
                errorbar=("sd", 1), err_kws={"lw" : 0.8}, capsize=.25, 
                palette="Greens", hue_order=["TAB", "Nanopore"],
                ax=ax0)
    del barplot_feature

    ax0.set_ylim(0, 30)

    ax0.set_xlabel(None)
    ax0.set_ylabel("Mean feature 5hmC %")

    sns.move_legend(ax0, "upper left", ncol=2, title=None, frameon=False)

    ax1 = fig.add_subplot(gs[0, 4:])

    sns.barplot(barplot_cgi, x="feature_type", y="value", 
                hue="variable", 
                errorbar=("sd", 1), err_kws={"lw" : 0.8}, capsize=.25,
                palette="Greens", hue_order=["TAB", "Nanopore"],
                ax=ax1)
    del barplot_cgi

    ax1.set_ylim(0, 30)

    ax1.set_xlabel(None)
    ax1.set_ylabel(None)

    sns.move_legend(ax1, "upper right", ncol=2, title=None, frameon=False)

    ax00 = fig.add_subplot(gs[1, :2])
    ax01 = fig.add_subplot(gs[1, 2:4])
    ax02 = fig.add_subplot(gs[1, 4:])

    ax10 = fig.add_subplot(gs[2, :2])
    ax11 = fig.add_subplot(gs[2, 2:4])
    ax12 = fig.add_subplot(gs[2, 4:])

    kde_plots = pd.concat([features_enrichment.groupby("feature_type").get_group(feature) for feature in ["Intron", "Exon (coding)", "Exon (all)", "Whole gene", "1kb promoter"]])
    
    # Also interested to include CGI in this comparison
    kde_plots = pd.concat([kde_plots, cgi_enrichment.groupby("feature_type").get_group("CGI")])

    kde_features = kde_plots.groupby("feature_type")

    for feature, ax in zip(["Whole gene", "Intron", "Exon (coding)", "Exon (all)", "1kb promoter", "CGI"], 
                                [ax00, ax01, ax02, ax10, ax11, ax12]):
        df = kde_features.get_group(f"{feature}")

        plot_kde(df, ax)
        ax.axhline(0, ls=":", c="grey", lw=0.8)
        ax.axvline(0, ls=":", c="grey", lw=0.8)

        ax.set_title(f"{feature}\nn={len(df)}, r={round(stats.pearsonr(df.enrichment_Nanopore, df.enrichment_TAB).statistic, 3)}")

        ax.set_ylim((-6, 4))
        ax.set_xlim((-6, 4))

        ax.set_ylabel(None)
        ax.set_xlabel(None)

        ax.set_xticks(np.arange(-6, 4, 2))
        ax.set_yticks(np.arange(-6, 4, 2))

    fig.supylabel("5hmC enrichment (Nanopore)", y=.35, x=0.01, ha="center", size=7)
    fig.supxlabel("5hmC enrichment (TAB-seq)", size=7)

    for i, ax in enumerate(fig.axes):
        ax.set_title(f"{string.ascii_lowercase[i]}", fontdict={"fontweight" : "bold"}, loc="left")

    sns.despine()

    if dryrun:
        outpath = "plots/tests/compare_features.png"
    else:
        outpath = "plots/compare_features.png"

    return fig.savefig(outpath) 

##### main function ##### 

if __name__=="__main__":
    parser = argparse.ArgumentParser(
                        prog = "compare_features",
                        description = "Compares 5hmC in features of the different datasets.")
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