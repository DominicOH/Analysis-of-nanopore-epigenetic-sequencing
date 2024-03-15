import argparse
import pandas as pd
import pyranges as pr
import subprocess
from AnalysisTools import common, helpers
import numpy as np
import concurrent.futures
import seaborn as sns
import matplotlib as mpl
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
    return tiled_df.drop(columns="CpG_count")

def assign_enrichment(grouped_df):
    grouped_df = grouped_df.eval("tile_5hmC = (N_5hmC/readCount)*100")
    mean = grouped_df["tile_5hmC"].mean()

    grouped_df = grouped_df.assign(enrichment = lambda r: np.log2((r["tile_5hmC"]+1)/(mean+1)))
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

    blacklist = pr.read_bed("data/hmedip/blacklist/mm39-blacklist.v2_adj.bed")
    belongs_on_blacklist = pr.PyRanges(
        pd.DataFrame({"Chromosome" : ["chr13"], 
                      "Start" : [99221000],
                      "End" : [99223000]})
    )
    
    public_hmedip_dirpath = "/mnt/data1/doh28/data/public/PRJNA134133_hMeDIP/macs2/ucsc_converted/"
    public_hmedip = pd.concat([table.assign(Replicate = i) for i, table in enumerate(load_hmedip(public_hmedip_dirpath))])
    public_hmedip_pr = (pr.PyRanges(public_hmedip).intersect(blacklist, invert=True)
                        .intersect(belongs_on_blacklist, invert=True))

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
    
    del nano_modbeds, tab_modbeds
    with concurrent.futures.ProcessPoolExecutor(2) as tile_executor:
        print("Tiling Nanopore WGS reference data.")
        tile_dfs = tile_executor.map(tile_wgs, [nanopore_average, tab_average])
        tile_dfs = map(remove_low_count_tiles, tile_dfs)

    del nanopore_average, tab_average
    gc.collect()

    nanopore_wgs_tiles, tab_tiles = [assign_enrichment(df) for df in tile_dfs]
    kde_background = pd.merge(nanopore_wgs_tiles, tab_tiles, on=["Chromosome", "Start", "End"], suffixes=["_Nanopore", "_TAB"])
    
    peak_overlay = (pr.PyRanges(kde_background)
                    .join(public_hmedip_pr, report_overlap=True)
                    .as_df()
                    .sort_values("Overlap", ascending=False)
                    .groupby(["Chromosome", "Start", "End"], observed=True)
                    .head(1))

    jg = sns.JointGrid(height=(120/25.4))
    sns.set_style("ticks")

    mpl.rc('font', size=5)

    sns.scatterplot(peak_overlay.sort_values("fold_enrichment"),
                        x="enrichment_TAB",
                        y="enrichment_Nanopore",
                        palette="OrRd", 
                        size="fold_enrichment", sizes=(.5, 10),
                        hue="fold_enrichment",
                        ax=jg.ax_joint)

    sns.kdeplot(peak_overlay.sort_values("fold_enrichment"),
                        x="enrichment_TAB",
                        color="#e34a33",
                        ax=jg.ax_marg_x)

    sns.kdeplot(peak_overlay.sort_values("fold_enrichment"),
                        y="enrichment_Nanopore",
                        color="#e34a33",
                        ax=jg.ax_marg_y)
    
    sns.kdeplot(kde_background, x="enrichment_TAB", y="enrichment_Nanopore", 
                ax=jg.ax_joint, color="#9ecae1")
      
    jg.ax_joint.axhline(1, ls=":", c="grey", lw=0.8)
    jg.ax_joint.axvline(1, ls=":", c="grey", lw=0.8)
    
    jg.set_axis_labels(xlabel="Mean TAB Enrichment", 
                       ylabel="Mean Nanopore Enrichment", size=7)
    
    # jg.set(yticklabels={"size" : 7}, xticklabels={"size" : 7})

    sns.move_legend(jg.ax_joint, "lower right", title="hMeDIP fold\nenrichment")
    print(len(peak_overlay.query("enrichment_Nanopore > 1"))/len(peak_overlay))
    print(len(peak_overlay.query("enrichment_TAB > 1"))/len(peak_overlay))

    print("Spearman:", stats.spearmanr(peak_overlay["enrichment_Nanopore"], 
                                       peak_overlay["fold_enrichment"]))

    if dryrun:
        jg.savefig("plots/tests/pcr_hmedip_comparison.png", dpi=600)
    
    else: 
        jg.savefig("plots/pcr_hmedip_comparison.png", dpi=600)

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