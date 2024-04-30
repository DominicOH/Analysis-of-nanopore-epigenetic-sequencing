import argparse
import pandas as pd
from AnalysisTools.common import *
from AnalysisTools import annotation_features
from AnalysisTools.helpers import timer
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn import metrics
from scipy import stats

def annotate_wrapper(df_list: list[pd.DataFrame], annotator: annotation_features.Annotator):
    readcount_tot = np.sum([df["readCount"].sum() for df in df_list])
    hmc_tot = np.sum([df["N_5hmC"].sum() for df in df_list])

    mean_value = (hmc_tot/readcount_tot)*100
    
    def group_feature(annotated_df: pd.DataFrame):
        grouped = annotated_df.groupby(["Start_Feature", "End_Feature", "feature_type"], 
                                       observed=True, as_index=True)
        features = grouped.sum(numeric_only=True)
        features["CpG_count"] = grouped["Start"].count()

        features = features.loc[features.eval("CpG_count > 9")]

        features.eval(f"percentMeth_5hmC = (N_5hmC/readCount)*100", inplace=True)
        features[f"enrichment_5hmC"] = np.log2((features["percentMeth_5hmC"]+1)/(mean_value+1))
        return features

    with concurrent.futures.ThreadPoolExecutor(len(df_list)) as executor:
        annotated_all = executor.map(annotator.annotate, df_list)
        annotated_all = [(annotation.drop(columns=["Strand"], errors="ignore")) 
                         for annotation in annotated_all]
        
    grouped_all = group_feature(pd.concat(annotated_all))

    return grouped_all.drop(columns=["Start", "End", "N_5hmC", "readCount", "CpG_count"])

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
def fig_main(dryrun, fontsize=5, figsize=None):
    # Loading data # 
    if dryrun:
        nano_path = "data/dryruns/modbeds/nanopore/"
        tab_path = "data/dryruns/modbeds/tab/"

    else:
        nano_path = "data/modbases/modbeds/nanopore/"
        tab_path = "data/modbases/modbeds/tab/"

    paths = [nano_path, tab_path]
    cols = [["readCount", "N_hmC"], ["readCount", "N_mod"]]

    with concurrent.futures.ProcessPoolExecutor(3) as fetch_executor:
        modbeds_list = [fetch_executor.submit(fetch_modbeds, path, cols) for path, cols in zip(paths, cols)]
        modbeds_list = [future.result() for future in modbeds_list]

    [modbed.rename(columns={"N_mod" : "N_5hmC"}, inplace=True) for modbed in modbeds_list[1]]

    # Data processing # 

    feature_annotator = annotation_features.Annotator("/mnt/data1/doh28/analyses/mouse_hydroxymethylome_analysis/feature_references/feature_comparison/")

    with concurrent.futures.ThreadPoolExecutor(2) as annotation_executor:
        annotated_data = annotation_executor.map(lambda df: annotate_wrapper(df, feature_annotator), modbeds_list)
        nanopore_annotated, tab_annotated = annotated_data

    del annotated_data, modbeds_list

    features_enrichment = (pd.merge(nanopore_annotated, tab_annotated, 
                                    how="inner", 
                                    on=["Start_Feature", "End_Feature", "feature_type"],
                                    suffixes=["_Nanopore", "_TAB"])
                                    .reset_index()
                                    .replace(["intron", "allExons", "genes", "1kbPromoter"],
                                             ["Intron", "Exon", "Genic", "Promoter"]
                                             ))
    
    print(features_enrichment.groupby("feature_type").mean())
    
    # Fig setup # 
    
    sns.set_style("ticks")
    mpl.rc('font', size=fontsize)

    if figsize:
        figsize = float(figsize)/25.4
    else:
        figsize = 3

    fg = sns.FacetGrid(features_enrichment, col="feature_type", 
                       hue="feature_type", palette="PuBuGn",
                       col_wrap=2, col_order=["Intron", "Exon", "Genic", "Promoter"],
                       height=figsize)
    
    fg.map_dataframe(sns.kdeplot, 
                     x="enrichment_5hmC_TAB", 
                     y="enrichment_5hmC_Nanopore",
                     fill=True)
    
    fg.set_xlabels("")
    fg.set_ylabels("")
    
    fg.map(plt.axhline, y=0, ls=":", c="grey")
    fg.map(plt.axvline, x=0, ls=":", c="grey")

    fg.set_titles("{col_name}")

    fg.figure.supxlabel(f"Log$_{2}$ FC from mean (TAB)", y=-.05)
    fg.figure.supylabel(f"Log$_{2}$ FC from mean (Nanopore)", x=-.05)

    for feature in ["Intron", "Exon", "Genic", "Promoter"]:
        ax = fg.axes_dict[feature]

        test_df = features_enrichment.query(f"feature_type == '{feature}'")
        stat = stats.spearmanr(test_df["enrichment_5hmC_TAB"], 
                               test_df["enrichment_5hmC_Nanopore"])
        print(feature, stat.pvalue)
        
        if stat.pvalue < 0.001:
            star = "***"
        elif stat.pvalue < 0.01:
            star = "**"
        elif stat.pvalue < 0.05:
            star = "*"
        
        ax.annotate(f"$\\rho$={round(stat.statistic, 3)}$^{{{star}}}$", 
                    xy=(-3.5, 2))

    # def calculate_spearman(feature):
    #     df = kde_features.get_group(feature)
    #     stat = stats.spearmanr(df.enrichment_Nanopore, df.enrichment_TAB)
    #     return [feature, round(stat.statistic, 3), round(stat.pvalue, 3)]

    if dryrun:
        fg.savefig("plots/tests/compare_features.png", dpi=600) 
    else:
        fg.savefig("plots/compare_features.png", dpi=600) 
        fg.savefig("plots/compare_features.svg", dpi=600, transparent=True) 

    return

##### main function ##### 

if __name__=="__main__":
    parser = argparse.ArgumentParser(
                        prog = "compare_features",
                        description = "Compares 5hmC in features of the different datasets.")
    parser.add_argument("-d ", "--dryrun", 
                        action="store_true", 
                        dest="dryrun", 
                        default=False,
                        required=False,
                        help="Whether a test output is produced.") 
    parser.add_argument("--figsize", 
                        dest="figsize", 
                        default=None,
                        required=False,
                        help="Size of the figure produced in mm") 
    parser.add_argument("--fontsize", 
                        dest="fontsize", 
                        default=5,
                        required=False,
                        help="Size of figure font.") 

    args = parser.parse_args()

    fig_main(args.dryrun, args.fontsize, args.figsize)
