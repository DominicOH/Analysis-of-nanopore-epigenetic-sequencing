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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pingouin as pg

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

def remove_low_count_tiles(tiled_df, threshold=10):
    tiled_df = tiled_df.query(f"CpG_count >= {threshold}")
    return tiled_df.drop(columns="CpG_count")

def calculate_zscore(grouped_df):
    grouped_df = grouped_df.eval("tile_5hmC = (N_5hmC/readCount)")
    grouped_df["asin"] = np.arcsin(grouped_df["tile_5hmC"])
    grouped_df["zscore"] = stats.zscore(grouped_df["asin"])
    return grouped_df.drop(columns=["N_5hmC", "readCount"])

def get_basecall_stats(df):
    calls = df.agg({
        "N_C" : sum, 
        "N_5mC" : sum, 
        "N_5hmC" : sum
        })
    call_stats = pd.DataFrame(calls, columns=["count"])
    call_stats = call_stats.assign(proportion = lambda r: (r.div(r.sum())))
    call_stats = call_stats.eval("percentage = proportion*100")
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
    annotated = hmedip_annotator.annotate(obj)
                 
    feature_counts = (annotated.groupby("Replicate")
                      .value_counts(["feature_type"], normalize=True, sort=False)
                      .reset_index(name="proportion")
                      .replace(["1kbPromoter", "intron", "allExons", "intergenic"], 
                               ["Promoter", "Intron", "Exon", "Intergenic"]))
    
    counts = annotated.groupby("Replicate", as_index=False).value_counts(["feature_type"], sort=False)
    feature_counts["count"] = counts["count"]
    return feature_counts

@helpers.timer
def fig_main(dryrun=True, fontsize=5, threshold=10):
    # Formatting # 
    sns.set_style("ticks")
    mpl.rc('font', size=fontsize)
    fig = plt.figure(figsize=(89/25.4, 140/25.4),   
                    dpi=600, layout="constrained")

    gs = GridSpec(4, 2, fig)

    ax1 = fig.add_subplot(gs[0, :]) # pyGenomeTracks
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1])
    ax4 = fig.add_subplot(gs[2:, :])

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

    blacklist = pr.read_bed("/mnt/data1/doh28/data/reference_genomes/mm39/mm39-blacklist.v2_sorted.bed")
    # these are regions not currently on the blacklist, inspected in igv, 
    # found to likely represent repetitive elements
    belongs_on_blacklist = pr.PyRanges(
        pd.DataFrame({"Chromosome" : ["chr13", "chr17", "chr17"], 
                      "Start" : [99221000, 4676277, 56506763],
                      "End" : [99223000, 4676915, 56507303]})
    )				
    narrow_peaks_pr = pr.PyRanges(narrow_peaks)

    hmedip_annotator = annotation_features.Annotator("feature_references/hmedip_comparison_features/")
    confident_peaks = narrow_peaks_pr.intersect(blacklist, invert=True).intersect(belongs_on_blacklist, invert=True)
    confident_peaks = hmedip_annotator.annotate(confident_peaks, 
                                                report_overlap=True)
    confident_peaks = (confident_peaks.sort_values("Overlap", ascending=False)
                                      .groupby(["Chromosome", "Start", "End"], observed=True)
                                      .head(1)
                                      .replace(["1kbPromoter", "intron", "allExons", "intergenic"], ["Promoter", "Intron", "Exon", "Intergenic"]))

    # Plot proportion of C basecalls peak vs. Nanopore WGS # 
    print("Loading Nanopore WGS reference data...")

    nano_path = "data/modbases/modbeds/nanopore_base5/"
    tab_path = "data/modbases/modbeds/tab/"
    
    nano_cols = ["readCount", "N_C", "N_mC", "N_hmC", "percentMeth_5hmC"]
    tab_cols = ["readCount", "N_mod", "percentMeth_mod"]

    with concurrent.futures.ThreadPoolExecutor(4) as fetch_executor:
        nano_modbeds_fut, tab_modbeds_fut = [fetch_executor.submit(common.fetch_modbeds, path, sel_cols, dryrun) for path, sel_cols in zip([nano_path, tab_path], [nano_cols, tab_cols])]
        nano_modbeds = [future for future in nano_modbeds_fut.result()]
        tab_modbeds = [future.rename(columns={"percentMeth_mod" : "percentMeth_5hmC",
                                              "N_mod" : "N_5hmC"}) for future in tab_modbeds_fut.result()]

    with concurrent.futures.ThreadPoolExecutor(4) as merge_executor:
        nanopore_average, tab_average = merge_executor.map(common.merge_positions, [nano_modbeds, tab_modbeds])
        nanopore_average, tab_average = [df.reset_index() for df in [nanopore_average, tab_average]]

    def extract_confident_peaks(df):
        pyrange = pr.PyRanges(df)
        confident_peaks_pr = pr.PyRanges(confident_peaks)

        return pyrange.intersect(confident_peaks_pr, strandedness=False, how= "first").as_df()
    
    hmedip_modbeds = common.fetch_modbeds("/mnt/data1/doh28/data/nanopore_hmc_validation/hmedip_nanopore/modkit_reports/modbeds/", ["N_C", "N_mC", "N_hmC"])
    hmedip_modbeds = map(extract_confident_peaks, hmedip_modbeds)

    print("Summarising all nanopore basecalls.")
    hmedip_call_stats = [get_basecall_stats(df).assign(Replicate = i) for i, df in enumerate(hmedip_modbeds)]
    wgs_call_stats = [get_basecall_stats(df).assign(Replicate = i) for i, df in enumerate(nano_modbeds)]

    callstats_bar_df = (pd.concat([pd.concat(hmedip_call_stats).assign(Method="hMeDIP"), 
                                   pd.concat(wgs_call_stats).assign(Method="WGS")])
                          .reset_index()
                          .rename(columns={"index" : "mod"}))

    callstats_bar_df = callstats_bar_df.replace(["N_C", "N_5mC", "N_5hmC"], ["C", "5mC", "5hmC"])
    
    wgs = pd.concat(wgs_call_stats).reset_index().rename(columns={"index" : "mod"})
    wgs_sum = wgs["count"].sum()
    wgs_mods = wgs.groupby("mod").sum().reset_index().replace(["N_C", "N_5mC", "N_5hmC"], ["C", "5mC", "5hmC"])
    wgs_mods["proportion"] = wgs_mods["count"].div(wgs_sum)

    ttest_df = callstats_bar_df.groupby(["mod", "Method"])

    for mod in ["C", "5mC", "5hmC"]:
        genome = ttest_df.get_group((mod, "WGS"))
        n_hmedip = ttest_df.get_group((mod, "hMeDIP"))
        print(mod, n_hmedip["proportion"].to_numpy())
        print("T-Test:\n", pg.ttest(genome["proportion"].to_numpy(), n_hmedip["proportion"].to_numpy(), paired=False))

    sns.barplot(callstats_bar_df.sort_values("Method"), 
                x="mod", y="percentage",
                hue="Method", palette="PuBuGn",
                width=.6,
                order=["C", "5mC", "5hmC"], 
                hue_order=["WGS", "hMeDIP"], 
                errorbar=("sd", 1), err_kws={"lw" : 1}, capsize=.25,            
                ax=ax3)
    
    sns.stripplot(callstats_bar_df.sort_values("Method"), 
                  x="mod", y="percentage", 
                  dodge=True, order=["C", "5mC", "5hmC"], 
                  hue="Method", hue_order=["WGS", "hMeDIP"], color="k",
                  legend=False, 
                  size=3,
                  ax=ax3)

    sns.move_legend(ax3, "lower center", ncols=2, title=None, frameon=False, bbox_to_anchor=(.5, 1))

    ax3.set_xlabel(None)
    ax3.set_ylabel("Percent of CpG-context basecalls")
                    
    with concurrent.futures.ProcessPoolExecutor(2) as tile_executor:
        print("Tiling Nanopore WGS reference data.")
        tile_dfs = tile_executor.map(tile_wgs, [nanopore_average, tab_average])
        tile_dfs = map(lambda df: remove_low_count_tiles(df, threshold), tile_dfs)

    nanopore_wgs_tiles, tab_tiles = [calculate_zscore(df) for df in tile_dfs]
    kde_background = pd.merge(nanopore_wgs_tiles, tab_tiles, on=["Chromosome", "Start", "End"], suffixes=["_Nanopore", "_TAB"])

    peak_overlay = (pr.PyRanges(kde_background)
                    .join(pr.PyRanges(confident_peaks), report_overlap=True, suffix="_Peak")
                    .as_df()
                    .sort_values("Overlap", ascending=False)
                    .groupby(["Chromosome", "Start_Peak", "End_Peak", "Replicate"], observed=True) 
                    .head(1)
                    .eval("Peak_length = End_Peak - Start_Peak")
                    .query("Overlap >= (Peak_length / 2)")) # At least half the peak must sit in the region
    
    sns.kdeplot(kde_background, 
                x="zscore_TAB",
                y="zscore_Nanopore",
                fill=True, 
                # alpha=.6, 
                color="#9ecae1",
                ax=ax4)

    ax4in = inset_axes(ax4, height="30%", width="30%", loc="upper left")

    sns.kdeplot(kde_background, 
                x="zscore_TAB",
                y="zscore_Nanopore",
                fill=True, 
                # alpha=.6, 
                color="#9ecae1",
                ax=ax4in)

    ax4in.tick_params("both", left=False, bottom=False, labelleft=False, labelbottom=False)
    ax4in.set_xlabel(None)
    ax4in.set_ylabel(None)

    sns.scatterplot(peak_overlay.sort_values("fold_enrichment"),
                    x="zscore_TAB",
                    y="zscore_Nanopore",
                    hue="fold_enrichment",
                    palette="OrRd", 
                    size="fold_enrichment", sizes=(1, 15),
                    ax=ax4)
    
    print("Peaks over Z0 (Nanopore):", len(peak_overlay.query("zscore_Nanopore > 0"))/len(peak_overlay))
    print("Peaks over Z0 (TAB):", len(peak_overlay.query("zscore_TAB > 0"))/len(peak_overlay))
    
    sns.move_legend(ax4, "lower right", title="hMeDIP fold\nenrichment")

    ax4.set_xlabel("5hmC Z-Score (TAB)")
    ax4.set_ylabel("5hmC Z-Score (Nanopore)")

    for ax in [ax4, ax4in]:
        ax.set_xlim((-3, 7))
        ax.set_ylim((-3, 7))

        ax.axhline(kde_background["zscore_TAB"].mean(), ls=":", c="grey", lw=0.8)
        ax.axvline(kde_background["zscore_Nanopore"].mean(), ls=":", c="grey", lw=0.8)

    sns.despine()

    public_hmedip_dirpath = "/mnt/data1/doh28/data/public/PRJNA134133_hMeDIP/ucsc_converted/"
    public_hmedip = pd.concat([table.assign(Replicate = i) for i, table in enumerate(load_hmedip(public_hmedip_dirpath))])

    g_tiles = pr.genomicfeatures.tile_genome(chromSizePr(), 500).as_df()
    g_tiles["Replicate"] = 0

    annotated_all = map(prep_feature_type_counts, [confident_peaks, g_tiles, public_hmedip])
    feature_barplot = pd.concat([df.assign(Method = method) for method, df in zip(["Direct hMeDIP", "Genome", "PCR hMeDIP"], annotated_all)])

    feature_barplot.eval("percentage = proportion * 100", inplace=True)
    sns.barplot(feature_barplot.reset_index(),
                x="feature_type", y="percentage",
                errorbar=("sd", 1), err_kws={"lw" : 0.8}, capsize=.25,
                order=["Intergenic", "Intron", "Exon", "Promoter"],
                hue="Method", hue_order=["Genome", "Direct hMeDIP", "PCR hMeDIP"],
                palette="BuPu",
                ax=ax2)
    
    sns.stripplot(feature_barplot.reset_index(),
                  x="feature_type", y="percentage",
                  dodge=True,
                  order=["Intergenic", "Intron", "Exon", "Promoter"],
                  hue="Method", hue_order=["Genome", "Direct hMeDIP", "PCR hMeDIP"],
                  legend=False, color="k",
                  size=3,
                  ax=ax2)
    
    feature_barplot_stats = feature_barplot.groupby(["Method", "feature_type"])
    for feature, hypothesis in zip(["Intergenic", "Intron", "Exon", "Promoter"], ["less", "greater", "greater", "greater"]):
        genome_p = feature_barplot_stats.get_group(("Genome", feature))["proportion"].to_list()[0]

        counts = feature_barplot_stats.get_group(("Direct hMeDIP", feature))["count"].sum()
        nobs = feature_barplot.loc[feature_barplot["Method"] == "Direct hMeDIP", "count"].sum()

        print("Binom:", feature, stats.binomtest(counts, nobs, genome_p, hypothesis))

    ax2.set_ylabel("Percentage of coverage")
    ax2.set_xlabel(None)

    sns.move_legend(ax2, loc="upper right", title=None, frameon=False)
        
    sns.despine()
    sns.despine(ax=ax1, left=True, bottom=True)
    sns.despine(ax=ax4in, right=False, top=False)

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
    parser.add_argument("--fontsize", 
                        action="store", 
                        dest="fontsize",
                        default=5, 
                        required=False) 
    parser.add_argument("--threshold", 
                        action="store", 
                        default=10,
                        dest="threshold", 
                        required=False)

    args = parser.parse_args()

    if args.dryrun:
        dryrun = True
    else: 
        dryrun = False

    fig_main(dryrun, args.fontsize, args.threshold)
    print("Completed.")