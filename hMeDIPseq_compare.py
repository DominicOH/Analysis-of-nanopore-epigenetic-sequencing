import argparse
import pandas as pd
import pyranges as pr
import subprocess
from AnalysisTools import annotation_features, common, helpers, multisite
import numpy as np
import string
import concurrent.futures
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from scipy import stats

def tile_wgs(df, tile_size=500):
    tiles = (pr.PyRanges(df.reset_index())
             .tile(tile_size)
             .as_df()
             .groupby(["Chromosome", "Start", "End"], observed=True)
             .agg({
                 "readCount" : np.sum,
                 "N_5hmC" : np.sum,
                 "Start" : "count"
             }).rename(columns={"Start" : "CpG_count"})
             .reset_index())
    return tiles

def remove_low_count_tiles(tiled_df, threshold=5):
    tiled_df = tiled_df.loc[tiled_df.eval("CpG_count >= @threshold")]
    return tiled_df.drop(columns="CpG_count")

def assign_enrichment(grouped_df):
    grouped_df = grouped_df.eval("tile_5hmC = (N_5hmC/readCount)*100")
    mean = grouped_df["tile_5hmC"].mean()

    grouped_df = grouped_df.assign(enrichment = lambda r: np.log2((r["tile_5hmC"]+1)/(mean+1)))
    return grouped_df.drop(columns=["N_5hmC", "readCount"])

def get_basecall_stats(df):
    calls = df.agg({
        "N_C" : sum, 
        "N_5mC" : sum, 
        "N_5hmC" : sum
        })
    call_stats = pd.DataFrame(calls, columns=["count"])
    call_stats = call_stats.assign(proportion = lambda r: (r.div(r.sum()))*100)
    return call_stats

@multisite.asPyRangesDecorator
def chromSizePr():
    return multisite.loadChromSize()

def read_narrowPeak(path):
    df = pd.read_table(path, 
                       sep="\t",
                       names=["Chromosome", "Start", "End", "Name", "Score", ".", "fold_enrichment", "PVAL", "QVAL", "lenght"],
                       usecols=["Chromosome", "Start", "End", "fold_enrichment", "PVAL", "QVAL"])
    return df

def load_hmedip(dirpath):
    path_ls = subprocess.check_output([f"ls {dirpath}*.narrowPeak"], shell=True).decode("utf-8").split("\n")
    path_ls.pop(-1)

    with concurrent.futures.ThreadPoolExecutor(len(path_ls)) as read_executor:
        peak_files = read_executor.map(read_narrowPeak, path_ls)
        peak_files = [df.assign(Replicate = i) for i, df in enumerate(peak_files)]
    return peak_files

def prep_feature_type_counts(obj):
    hmedip_annotator = annotation_features.Annotator("feature_references/hmedip_comparison_features/")
    annotated = (hmedip_annotator.annotate(obj)
                 .groupby("Replicate")
                 .value_counts(["feature_type"], normalize=True)
                 .reset_index(name="proportion")
                 .replace(["1kbPromoter", "intron", "allExons", "intergenic"], 
                          ["Promoter", "Intron", "Exon", "Intergenic"]))
    return annotated

@helpers.timer
def fig_main(dryrun=True):
    # Formatting # 
    sns.set_style("ticks")
    mpl.rc('font', size=5)
    fig = plt.figure(figsize=(180/25.4, 89/25.4),   
                    dpi=600, layout="constrained")

    gs = GridSpec(2, 4, fig)

    ax1 = fig.add_subplot(gs[0, :2]) # pyGenomeTracks
    ax2 = fig.add_subplot(gs[1, :2])
    ax3 = fig.add_subplot(gs[:, 2])
    ax4 = fig.add_subplot(gs[0, 3])
    ax5 = fig.add_subplot(gs[1, 3])

    for i, ax in enumerate(fig.axes):
        ax.set_title(string.ascii_lowercase[i], fontdict={"weight" : "bold",
                                                        "fontsize" : 5},
                                                        loc="left")

    ax1.tick_params("both", bottom=False, left=False, 
                    labelbottom=False, labelleft=False)
                
    # Peaks vs. Nanopore WGS
        
    print("Loading MACS2 peak data.")
    nanopore_hmedip_dirpath = "/mnt/data1/doh28/data/nanopore_hmc_validation/hmedip_nanopore/summary/"
    narrow_peaks = pd.concat(load_hmedip(nanopore_hmedip_dirpath))

    blacklist = pr.read_bed("data/hmedip/blacklist/mm39-blacklist.v2_adj.bed")
    narrow_peaks_pr = pr.PyRanges(narrow_peaks)
                                                        
    hmedip_annotator = annotation_features.Annotator("feature_references/hmedip_comparison_features/")
    confident_peaks = hmedip_annotator.annotate(narrow_peaks_pr.intersect(blacklist, invert=True), report_overlap=True)
    confident_peaks = (confident_peaks.sort_values("Overlap", ascending=False)
                                    .groupby(["Chromosome", "Start", "End"], observed=True)
                                    .head(1)
                                    .replace(["1kbPromoter", "intron", "allExons", "intergenic"], ["Promoter", "Intron", "Exon", "Intergenic"]))

    blacklisted_peaks = narrow_peaks_pr.intersect(blacklist)

    # Plot peak confidence scatter


    sns.scatterplot(blacklisted_peaks.as_df(), 
                    x="fold_enrichment", y="QVAL",
                    marker="x", 
                    label="Blacklist",
                    s=10,
                    color=sns.color_palette("Paired", 6)[5],
                    ax=ax3)

    sns.scatterplot(confident_peaks.sort_values("feature_type"),
                    x="fold_enrichment", y="QVAL",
                    s=10, hue="feature_type", style="feature_type",
                    palette="Greens",
                    ax=ax3)

    ax3.set_ylabel(r"$-log_{10}$ q-value")
    ax3.set_xlabel("Fold enrichment over input")
    sns.move_legend(ax3, "upper left", title=None)

    # Plot proportion of C basecalls peak vs. Nanopore WGS # 
    print("Loading Nanopore WGS reference data...")

    if dryrun:
        nano_path = "data/dryruns/modbeds/nanopore/"
        tab_path = "data/dryruns/modbeds/tab/"

    else:
        nano_path = "data/modbases/modbeds/nanopore/"
        tab_path = "data/modbases/modbeds/tab/"
        
    nano_cols = ["readCount", "N_C", "N_mC", "N_hmC", "percentMeth_5hmC"]
    tab_cols = ["readCount", "N_mod", "percentMeth_mod"]

    with concurrent.futures.ThreadPoolExecutor(4) as fetch_executor:
        nano_modbeds_fut, tab_modbeds_fut = [fetch_executor.submit(common.fetch_modbeds, path, sel_cols) for path, sel_cols in zip([nano_path, tab_path], [nano_cols, tab_cols])]
        nano_modbeds = [future for future in nano_modbeds_fut.result()]
        tab_modbeds = [future.rename(columns={"percentMeth_mod" : "percentMeth_5hmC",
                                                "N_mod" : "N_5hmC"}) for future in tab_modbeds_fut.result()]

    with concurrent.futures.ThreadPoolExecutor(4) as merge_executor:
        nanopore_average, tab_average = merge_executor.map(common.merge_positions, [nano_modbeds, tab_modbeds])
        nanopore_average, tab_average = [df.reset_index() for df in [nanopore_average, tab_average]]

    hmedip_modbeds = common.fetch_modbeds("/mnt/data1/doh28/data/nanopore_hmc_validation/hmedip_nanopore/modkit_reports/modbeds/", ["N_C", "N_mC", "N_hmC"])

    print("Summarising all nanopore basecalls.")
    hmedip_call_stats = [get_basecall_stats(df).assign(Replicate = i) for i, df in enumerate(hmedip_modbeds)]
    wgs_call_stats = [get_basecall_stats(df).assign(Replicate = i) for i, df in enumerate(nano_modbeds)]

    callstats_bar_df = pd.concat([pd.concat(hmedip_call_stats).assign(Method="All hMeDIP"), 
                                  pd.concat(wgs_call_stats).assign(Method="WGS"),
                                  pd.concat(hmedip_call_stats).query("Replicate > 0").assign(Method="Monoclonal")]).reset_index().rename(columns={"index" : "mod"})
    callstats_bar_df = callstats_bar_df.replace(["N_C", "N_5mC", "N_5hmC"], ["C", "5mC", "5hmC"])
    
    def stat_arrays(df, groups, target):
        array = [df.groupby(groups)
        .get_group((mod, target))["proportion"]
        .to_numpy() for mod in ["C", "5mC", "5hmC"]]
        return array

    all_hmedips = stat_arrays(callstats_bar_df, ["mod", "Method"], "All hMeDIP")
    wgss = stat_arrays(callstats_bar_df, ["mod", "Method"], "WGS")
    monoclonal_only = stat_arrays(callstats_bar_df, ["mod", "Method"], "Monoclonal")

    print("Welch's T-Test: all hMeDIPs with WGS:")
    for mod, hmedip, wgs in zip(["C", "5mC", "5hmC"], all_hmedips, wgss):
        # Have to assume unequal variance due to the presence of the polyclonal antibody
        stat = stats.ttest_ind(hmedip, wgs, equal_var=False).pvalue
        print(f"{mod}: p={round(stat, 6)}")

    print("Unpaired T-Test: monoclonal hMeDIPS with WGS:")
    for mod, hmedip, wgs in zip(["C", "5mC", "5hmC"], monoclonal_only, wgss):
        # We can assume equal variance because all WGS are consistent replicates of littermate mice
        # And the hMeDIP are produced using the same antibody 
        stat = stats.ttest_ind(hmedip, wgs, equal_var=True).pvalue
        print(f"{mod}: p={stat}")
    
    ax5.axhline(0, c=(0.15, 0.15, 0.15, 1.0), lw=0.8)
    sns.barplot(callstats_bar_df.sort_values("Method"), 
                x="mod",
                hue="Method",
                y="proportion",
                palette="Greens",
                width=.6,
                order=["C", "5mC", "5hmC"], 
                hue_order=["WGS", "All hMeDIP", "Monoclonal"], 
                errorbar=("sd", 1), err_kws={"lw" : 0.8}, capsize=.25,            
                ax=ax5)
    
    # ax5.annotate("*", (2, 35), size=7)
    sns.move_legend(ax5, "lower center", ncols=2, title=None, frameon=False, reverse=True, bbox_to_anchor=(.5, 1))

    ax5.set_xlabel(None)
    ax5.set_ylabel("% of CpG-context basecalls")

    # hMeDIP vs. Nanopore WGS enriched regions
                    
    with concurrent.futures.ProcessPoolExecutor(2) as tile_executor:
        print("Tiling Nanopore WGS reference data.")
        tile_dfs = tile_executor.map(tile_wgs, [nanopore_average, tab_average])
        tile_dfs = map(remove_low_count_tiles, tile_dfs)

    nanopore_wgs_tiles, tab_tiles = [assign_enrichment(df) for df in tile_dfs]
    kde_background = pd.merge(nanopore_wgs_tiles, tab_tiles, on=["Chromosome", "Start", "End"], suffixes=["_Nanopore", "_TAB"])

    peak_overlay = (pr.PyRanges(kde_background)
                    .join(narrow_peaks_pr.intersect(blacklist, invert=True), report_overlap=True)
                    .as_df()
                    .sort_values("Overlap", ascending=False)
                    .groupby(["Chromosome", "Start", "End"], observed=True)
                    .head(1))

    sns.kdeplot(kde_background, 
                x="enrichment_TAB",
                y="enrichment_Nanopore",
                fill=True, 
                color=sns.color_palette("Greens", 5)[4],
                ax=ax2)

    sns.scatterplot(peak_overlay.sort_values("fold_enrichment"),
                    x="enrichment_TAB",
                    y="enrichment_Nanopore",
                    palette="rocket_r", 
                    size="fold_enrichment", sizes=(1, 15),
                    hue="fold_enrichment",
                    ax=ax2)

    sns.move_legend(ax2, "lower right", title="hMeDIP fold\nenrichment")

    ax2.set_xlabel("Mean TAB Enrichment")
    ax2.set_ylabel("Mean Nanopore Enrichment")

    ax2.set_xlim((-4, 4))
    ax2.set_ylim((-4, 4))

    ax2.axhline(0, ls=":", c="grey", lw=0.8)
    ax2.axvline(0, ls=":", c="grey", lw=0.8)

    public_hmedip_dirpath = "/mnt/data1/doh28/data/public/PRJNA134133_hMeDIP/macs2/ucsc_converted/"
    public_hmedip = pd.concat([table.assign(Replicate = i) for i, table in enumerate(load_hmedip(public_hmedip_dirpath))])

    g_tiles = pr.genomicfeatures.tile_genome(chromSizePr(), 500).as_df()
    g_tiles["Replicate"] = 0

    annotated_all = map(prep_feature_type_counts, [confident_peaks, g_tiles, public_hmedip])
    feature_barplot = pd.concat([df.assign(Method = method) for method, df in zip(["Nanopore", "Genome", "Ref."], annotated_all)])

    feature_barplot.eval("proportion = proportion * 100", inplace=True)
    sns.barplot(feature_barplot.reset_index(),
                x="feature_type", 
                y="proportion",
                errorbar=("sd", 1), err_kws={"lw" : 0.8}, capsize=.25,
                order=["Intergenic", "Intron", "Exon", "Promoter"],
                hue="Method", hue_order=["Genome", "Ref.", "Nanopore"],
                palette="Greens",
                ax=ax4)

    ax4.set_ylabel("Percentage of coverage")
    ax4.set_xlabel(None)

    sns.move_legend(ax4, loc="upper right", title=None, frameon=False)
        
    sns.despine()
    sns.despine(ax=ax1, left=True, bottom=True)

    print("Done. Saving...")
    
    if dryrun:
        fig.savefig("plots/tests/f4_canvas_w_plots.svg")
        fig.savefig("plots/tests/f4_canvas_w_plots.png")
    else:
        fig.savefig("plots/f4_canvas_w_plots.svg")
        fig.savefig("plots/f4_canvas_w_plots.png")

##### main function ##### 

if __name__=="__main__":
    parser = argparse.ArgumentParser(
                        prog = "compare_hmedip",
                        description = "Compares Nanopore hMeDIP to Nanopore WGS and standard hMeDIP.")
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
    print("Completed.")