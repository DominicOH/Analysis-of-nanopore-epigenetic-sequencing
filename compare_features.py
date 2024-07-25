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
    def calculate_zscore(annotated_df: pd.DataFrame):
        """
        Z-Score is here calculated after annotation. Z-Score here calculated element-wise rather than CpG-wise. 
        """
        annotated_df.eval("percentMeth_5hmC = (N_5hmC/readCount)", inplace=True)
        annotated_df[f"asin_5hmC"] = np.arcsin(annotated_df[f"percentMeth_5hmC"])
        annotated_df["zscore_5hmC"] = stats.zscore(annotated_df["asin_5hmC"])
        
        return annotated_df   
     
    def group_feature(annotated_df: pd.DataFrame):
        grouped = annotated_df.groupby(["Start_Feature", "End_Feature", "feature_type"], 
                                       observed=True, as_index=True)
        
        features = grouped.sum(numeric_only=True)
        features["CpG_count"] = grouped["Start"].count()

        features = features.loc[features.eval("CpG_count > 9")]
        features = calculate_zscore(features)

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
def fig_main(dryrun, fontsize=5):
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
    gene_annotator = annotation_features.Annotator("/mnt/data1/doh28/analyses/mouse_hydroxymethylome_analysis/feature_references/genic/")

    with concurrent.futures.ThreadPoolExecutor(2) as annotation_executor:
        annotated_data = annotation_executor.map(lambda df: annotate_wrapper(df, feature_annotator), modbeds_list)
        nanopore_annotated, tab_annotated = annotated_data

    del annotated_data

    features_enrichment = (pd.merge(nanopore_annotated, tab_annotated, 
                                    how="inner", 
                                    on=["Start_Feature", "End_Feature", "feature_type"],
                                    suffixes=["_Nanopore", "_TAB"])
                                    .reset_index()
                                    .replace(["intron", "allExons", "1kbPromoter"],
                                             ["Intron", "Exon", "Promoter"]
                                             ))
    
    with concurrent.futures.ThreadPoolExecutor(2) as annotation_executor:
        annotated_data = annotation_executor.map(lambda df: annotate_wrapper(df, gene_annotator), modbeds_list)
        nanopore_annotated, tab_annotated = annotated_data

    gene_enrichment = (pd.merge(nanopore_annotated, tab_annotated, 
                                    how="inner", 
                                    on=["Start_Feature", "End_Feature", "feature_type"],
                                    suffixes=["_Nanopore", "_TAB"])
                                    .reset_index()
                                    .replace(["genes"], ["Gene body"]))
    gene_enrichment = gene_enrichment.loc[gene_enrichment["feature_type"] == "Gene body"]

    all_enrichment = pd.concat([features_enrichment, gene_enrichment])
        
    # Fig setup # 
    
    sns.set_style("ticks")
    mpl.rc('font', size=fontsize)


    fg = sns.FacetGrid(all_enrichment, col="feature_type", 
                       hue="feature_type", palette="PuBuGn",
                       col_wrap=3, col_order=["Gene body", "Promoter", "5UTR", "Intron", "Exon", "3UTR"],
                       xlim=(-3, 4), ylim=(-3, 4))
    
    fg.figure.set_size_inches(180/25.4, 120/25.4)
    fg.figure.set_constrained_layout(True)
    
    fg.map_dataframe(sns.kdeplot, 
                     x="zscore_5hmC_TAB", 
                     y="zscore_5hmC_Nanopore",
                     fill=True)
    
    fg.set_xlabels("")
    fg.set_ylabels("")
    
    fg.map(plt.axhline, y=0, ls=":", c="grey")
    fg.map(plt.axvline, x=0, ls=":", c="grey")

    fg.set_titles("{col_name}")

    fg.figure.supxlabel(f"Site modification (%) Z-Score (TAB)", y=-.01)
    fg.figure.supylabel(f"Site modification (%) Z-Score (Nanopore)", x=-.01)

    for feature in ["Intron", "Exon", "Gene body", "Promoter", "5UTR", "3UTR"]:
        ax = fg.axes_dict[feature]

        test_df = all_enrichment.query(f"feature_type == '{feature}'")
        stat = stats.spearmanr(test_df["zscore_5hmC_TAB"], 
                               test_df["zscore_5hmC_Nanopore"])
        print(feature, len(test_df), stat)
        
        if stat.pvalue < 0.0001:
            star = "****"
        elif stat.pvalue < 0.001:
            star = "***"
        elif stat.pvalue < 0.01:
            star = "**"
        elif stat.pvalue < 0.05:
            star = "*"
        
        ax.annotate(f"$\\rho$={round(stat.statistic, 3)}$^{{{star}}}$", 
                    xy=(-2, 3))

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
    parser.add_argument("--fontsize", 
                        dest="fontsize", 
                        default=5,
                        required=False,
                        help="Size of figure font.") 

    args = parser.parse_args()

    fig_main(args.dryrun, args.fontsize)
