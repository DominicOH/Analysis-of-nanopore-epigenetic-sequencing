import argparse
import pandas as pd
import pyranges as pr
import subprocess
from AnalysisTools import common, helpers
import numpy as np
import concurrent.futures
import seaborn as sns
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import stats
import gc

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
    return tiled_df

def calculate_zscore(grouped_df):
    grouped_df = grouped_df.eval("tile_5hmC = (N_5hmC/readCount)")
    grouped_df = grouped_df.eval("mod_density_100bp = (N_5hmC/(readCount * 5))")
    print(f"Average 5hmCpG density per 100bp: {grouped_df['mod_density_100bp'].mean()}")
    grouped_df["asin"] = np.arcsin(grouped_df["tile_5hmC"])
    grouped_df["zscore"] = stats.zscore(grouped_df["asin"])

    enriched_density = grouped_df.loc[grouped_df["zscore"] > 0, ("mod_density_100bp")].mean()
    print(f"Average 5hmCpG density per 100bp in enriched Z>0 regions: {enriched_density}")
    return grouped_df.drop(columns=["N_5hmC", "readCount"])

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

@helpers.timer
def fig_main(dryrun=True):

    blacklist = pr.read_bed("/mnt/data1/doh28/data/reference_genomes/mm39/mm39-blacklist.v2_sorted.bed")
    belongs_on_blacklist = pr.PyRanges(
        pd.DataFrame({"Chromosome" : ["chr13"], 
                      "Start" : [99221000],
                      "End" : [99223000]})
    )
    
    public_hmedip_dirpath = "/mnt/data1/doh28/data/public/PRJNA134133_hMeDIP/ucsc_converted/"
    public_hmedip = pd.concat([table.assign(Replicate = i) for i, table in enumerate(load_hmedip(public_hmedip_dirpath))])

    print("Total length of public hMeDIP peaks: ", public_hmedip.apply(lambda r: r["End"] - r["Start"], axis=1).sum())
    public_hmedip_pr = (pr.PyRanges(public_hmedip).intersect(blacklist, invert=True)
                        .intersect(belongs_on_blacklist, invert=True))
        
    usecols = [["readCount", "N_C", "N_mC", "N_hmC", "percentMeth_5hmC"], ["readCount", "N_mod", "percentMeth_mod"]]
    paths =  ["data/modbases/modbeds/nanopore_base5/", "data/modbases/modbeds/tab/"]

    with concurrent.futures.ThreadPoolExecutor(4) as fetch_executor:
        nano_modbeds_fut, tab_modbeds_fut = [fetch_executor.submit(common.fetch_modbeds, path, sel_cols, dryrun) for path, sel_cols in zip(paths, usecols)]
        nano_modbeds = [future for future in nano_modbeds_fut.result()]
        tab_modbeds = [future.rename(columns={"percentMeth_mod" : "percentMeth_5hmC",
                                              "N_mod" : "N_5hmC"}) for future in tab_modbeds_fut.result()]

    with concurrent.futures.ThreadPoolExecutor(4) as merge_executor:
        nanopore_average, tab_average = merge_executor.map(common.merge_positions, [nano_modbeds, tab_modbeds])
        nanopore_average, tab_average = [df.reset_index() for df in [nanopore_average, tab_average]]
    
    print("Loaded data")
    del nano_modbeds, tab_modbeds
    
    with concurrent.futures.ProcessPoolExecutor(2) as tile_executor:
        print("Tiling Nanopore WGS reference data.")
        tile_dfs = tile_executor.map(tile_wgs, [nanopore_average, tab_average])
        tile_dfs = map(remove_low_count_tiles, tile_dfs)

    del nanopore_average, tab_average
    gc.collect()

    nanopore_wgs_tiles, tab_tiles = [calculate_zscore(df) for df in tile_dfs]
    kde_background = pd.merge(nanopore_wgs_tiles, tab_tiles, on=["Chromosome", "Start", "End"], suffixes=["_Nanopore", "_TAB"])                

    peak_overlay = (pr.PyRanges(kde_background)
                    .join(public_hmedip_pr, report_overlap=True, suffix="_Peak")
                    .as_df()
                    .sort_values("Overlap", ascending=False)
                    .groupby(["Chromosome", "Start_Peak", "End_Peak", "Replicate"], observed=True)
                    .head(1)
                    .eval("Peak_length = End_Peak - Start_Peak")
                    .query("Overlap >= (Peak_length / 2)")) # At least half the peak must sit in the region)

    jg = sns.JointGrid(height=(120/25.4), xlim=(-3, 7), ylim=(-3, 7))
    jg.figure.set_constrained_layout(True)
    
    sns.set_style("whitegrid")
    mpl.rc('font', size=5)

    print("Plotting")
        
    sns.kdeplot(kde_background, x="zscore_TAB", y="zscore_Nanopore", 
                ax=jg.ax_joint, color="#9ecae1", fill=True)
    
    with pd.ExcelWriter('source_data/figs4_peak_overlay.xlsx') as writer:
        peak_overlay.to_excel(writer, 'figs4_peak_overlay')

    sns.scatterplot(peak_overlay.sort_values("fold_enrichment"),
                        x="zscore_TAB",
                        y="zscore_Nanopore",
                        palette="OrRd", 
                        marker=".",
                        size="fold_enrichment", sizes=(.5, 10),
                        hue="fold_enrichment",
                        ax=jg.ax_joint)

    sns.histplot(peak_overlay,
                 x="zscore_TAB",
                 color="#e34a33",
                 ax=jg.ax_marg_x)

    sns.histplot(peak_overlay,
                 y="zscore_Nanopore",
                 color="#e34a33",
                 ax=jg.ax_marg_y)
      
    jg.ax_joint.axhline(kde_background["zscore_TAB"].mean(), ls=":", c="grey", lw=0.8)
    jg.ax_joint.axvline(kde_background["zscore_Nanopore"].mean(), ls=":", c="grey", lw=0.8)
    jg.set_axis_labels(xlabel="Site modification (%) Z-Score (TAB)", 
                       ylabel="Site modification (%) Z-Score (Nanopore)", size=7)
    
    axin = inset_axes(jg.ax_joint, width="20%", height="20%", loc="upper right")
    axin.tick_params("both", left=False, bottom=False, labelleft=False, labelbottom=False)
    axin.set_xlim((-3, 7))
    axin.set_ylim((-3, 7))

    sns.kdeplot(kde_background, x="zscore_TAB", y="zscore_Nanopore", 
                ax=axin, color="#9ecae1", fill=True)
    
    axin.axhline(kde_background["zscore_TAB"].mean(), ls=":", c="grey", lw=0.8)
    axin.axvline(kde_background["zscore_Nanopore"].mean(), ls=":", c="grey", lw=0.8)

    stat = stats.spearmanr(kde_background["zscore_Nanopore"], 
                           kde_background["zscore_TAB"])
    
    print("Spearman (tiles):", stat, "n=", len(kde_background))
        
    if stat.pvalue < 0.0001:
        star = "****"
    elif stat.pvalue < 0.001:
        star = "***"
    elif stat.pvalue < 0.01:
        star = "**"
    elif stat.pvalue < 0.05:
        star = "*"
    
    axin.annotate(f"$\\rho$={round(stat.statistic, 3)}{star}", 
                  xy=(.5, 5.5))
    
    sns.move_legend(jg.ax_joint, "lower right", title="hMeDIP fold\nenrichment")
    print("Peaks over Z0 (Nanopore):", len(peak_overlay.query("zscore_Nanopore > 0"))/len(peak_overlay))
    print("Peaks over Z0 (TAB):", len(peak_overlay.query("zscore_TAB > 0"))/len(peak_overlay))

    print("Spearman (peaks vs. Nanopore):", stats.spearmanr(peak_overlay["zscore_Nanopore"], 
                                                            peak_overlay["fold_enrichment"]), 
                                                            "n=", len(peak_overlay))
    print("Spearman (peaks vs. TAB):", stats.spearmanr(peak_overlay["zscore_TAB"], 
                                                       peak_overlay["fold_enrichment"]), 
                                                       "n=", len(peak_overlay))
    
    print("What proportion of enriched regions in (Nanopore) WGS are detected as hMeDIP peaks?")
    l1 = pr.PyRanges(kde_background.query("zscore_Nanopore > 0")).intersect(public_hmedip_pr, how="first", invert=True) # regions that are enriched AND NOT intersect peaks
    l2 = len(kde_background.query("zscore_Nanopore > 0")) # regions that are enriched
    print(len(l1)/l2)
    l1_bp_len = l1.as_df().apply(lambda r: r["End"] - r["Start"], axis=1).sum()
    print(f"Equivalent to {l1_bp_len} bp")

    jg.savefig("plots/pcr_hmedip_comparison.png", dpi=300)    
    jg.savefig("plots/pcr_hmedip_comparison.svg", dpi=300)

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