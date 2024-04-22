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
from sklearn import metrics
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
    means = [df["percentMeth_5hmC"].mean() for df in df_list]
    print(means)

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

def fetch_and_merge(path, cols):
    dfs = fetch_modbeds(path, cols)
    dfs = [df.rename(columns={"N_mod":"N_5hmC"}, errors="ignore") for df in dfs]
    return merge_positions(dfs, calculate_percentage=True, cols=["N_5hmC"])

def cohen(x, y):
    mean_diff = np.mean(x) - np.mean(y)
    std = np.sqrt((np.std(x, ddof=1)**2 + np.std(y, ddof=1)**2) / 2)
    return mean_diff / std

def calculate_stats(df):
    out = {}
    for feature_type in df["feature_type"].unique():
        x = (df.groupby(["variable", "feature_type"])
             .get_group(("TAB", feature_type))["value"]
             .to_numpy())
        y = (df.groupby(["variable", "feature_type"])
             .get_group(("Nanopore", feature_type))["value"]
             .to_numpy())
        mw = stats.mannwhitneyu(x, y)
        cohen_stat = cohen(x, y)
        rmsd = metrics.mean_squared_error(x, y, squared=False)
        out.update({feature_type : [round(n, 4) for n in [cohen_stat, rmsd, mw.pvalue]]})
    return out


@timer
def fig_main(dryrun):
    # Loading data # 
    if dryrun:
        nano_path = "data/dryruns/modbeds/nanopore/"
        tab_path = "data/dryruns/modbeds/tab/"

    else:
        nano_path = "data/modbases/modbeds/nanopore/"
        tab_path = "data/modbases/modbeds/tab/"

    paths = [nano_path, tab_path]
    cols = [["readCount", "N_hmC"], ["readCount", "N_mod"]]

    with concurrent.futures.ProcessPoolExecutor(2) as fetch_executor:
        modbeds_future = [fetch_executor.submit(fetch_and_merge, path, cols) for path, cols in zip(paths, cols)]
        dfs = [future.result().reset_index() for future in modbeds_future]
    dfs = [df.assign(Method = method) for df, method in zip(dfs, ["Nanopore", "TAB"])]

    # Data processing # 

    feature_annotator = annotation_features.Annotator("/mnt/data1/doh28/analyses/mouse_hydroxymethylome_analysis/feature_references/feature_comparison/")
    cgi_annotator = annotation_features.Annotator("/mnt/data1/doh28/analyses/mouse_hydroxymethylome_analysis/feature_references/cgi/")

    with concurrent.futures.ProcessPoolExecutor(2) as annotation_executor:
        annotation_futures = [annotation_executor.submit(annotate_main, dfs, annotator) for annotator in [feature_annotator, cgi_annotator]]
        features_enrichment, cgi_enrichment = [future.result() for future in annotation_futures]

    del dfs

    features_enrichment = features_enrichment.replace(["intron", "allExons", "genes", "1kbPromoter"],
                                                  ["Intron", "Exon", "Genic", "Promoter"])

    features_enrichment["feature_type"] = pd.Categorical(features_enrichment.feature_type, 
                                                         ["Promoter", "Intron", "Exon", "Whole gene"], 
                                                         ordered=True)
    cgi_enrichment["feature_type"] = pd.Categorical(cgi_enrichment.feature_type, 
                                                    ["Shelf", "Shore", "CGI"], 
                                                    ordered=True)
    
    # Fig setup # 

    fig = plt.figure(figsize=(120/25.4, 100/25.4), 
                     dpi=600, 
                     layout="constrained")
    gs = GridSpec(1, 6, fig)

    sns.set_style("ticks")
    mpl.rc('font', size=5)

    ax00 = fig.add_subplot(gs[0, :2])
    ax01 = fig.add_subplot(gs[0, 2:4])
    ax02 = fig.add_subplot(gs[0, 4:])
    ax03 = fig.add_subplot(gs[0, :2])

    # Also interested to include CGI in this comparison
    kde_plots = pd.concat([features_enrichment, cgi_enrichment])
    kde_features = kde_plots.groupby("feature_type")

    def kde_outer(feature, ax):
        df = kde_features.get_group(f"{feature}")

        plot_kde(df, ax)
        ax.axhline(0, ls=":", c="grey", lw=0.8)
        ax.axvline(0, ls=":", c="grey", lw=0.8)

        stat = stats.spearmanr(df.enrichment_Nanopore, df.enrichment_TAB)
        sign = "="
        pval = round(stat.pvalue, 3)
        if pval < 0.001:
            pval = 0.001
            sign = "<"

        ax.set_title(f"{feature} (n={len(df)})\n\N{GREEK SMALL LETTER RHO}={round(stat.statistic, 3)}, p{sign}{round(pval, 3)}")

        ax.set_ylim((-6, 4))
        ax.set_xlim((-6, 4))

        ax.set_ylabel(None)
        ax.set_xlabel(None)

        ax.set_xticks(np.arange(-6, 4, 2))
        ax.set_yticks(np.arange(-6, 4, 2))
        return print(f"Plotted {feature}.")
    
    kde_targets = ["Whole gene", "Intron", "CDS", "Exon (all)", "Promoter", "CGI"]
    axes = [ax00, ax01, ax02, ax10, ax11, ax12]

    with concurrent.futures.ThreadPoolExecutor(len(kde_targets)) as kde_executor:
        kde_futures = [kde_executor.submit(kde_outer, feature, ax) for feature, ax in zip(kde_targets, axes)]
        [kde_plot.result() for kde_plot in kde_futures]

    print("Done plotting.")

    def calculate_spearman(feature):
        df = kde_features.get_group(feature)
        stat = stats.spearmanr(df.enrichment_Nanopore, df.enrichment_TAB)
        return [feature, round(stat.statistic, 3), round(stat.pvalue, 3)]

    with concurrent.futures.ThreadPoolExecutor(len(kde_features.groups.keys())) as spearman_executor:
        spearman = spearman_executor.map(calculate_spearman, kde_features.groups.keys()) 
        [print(stat) for stat in spearman]       

    fig.supylabel("5hmC enrichment (Nanopore)", y=.35, x=0.01, ha="center", size=7)
    fig.supxlabel("5hmC enrichment (TAB-seq)", size=7)

    for i, ax in enumerate(fig.axes):
        ax.set_title(f"{string.ascii_lowercase[i]}", fontdict={"fontweight" : "bold"}, loc="left")

    sns.despine()

    if dryrun:
        fig.savefig("plots/tests/compare_features.png") 
    else:
        fig.savefig("plots/compare_features.png") 
        fig.savefig("plots/compare_features.svg") 

    return print("Completed")

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