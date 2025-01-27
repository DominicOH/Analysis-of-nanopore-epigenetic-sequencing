import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
import concurrent.futures
import argparse
from AnalysisTools.helpers import timer
import string
import numpy as np
from scipy import stats
import AnalysisTools.ctcf_site_tools as ctcf
import gc
import scikit_posthocs as sp
import pingouin as pg

def produce_violins(file_path, replicate):
    modbed = ctcf.read_merge(file_path, test_run, replicate)
    violin = modbed.site_summit_distances(min_depth, filter_distance)
    return violin

@timer
def main():
    fig = plt.figure(figsize=(89/25.4, 120/25.4), dpi=600, layout="constrained")

    sns.set_style("ticks")
    mpl.rc('font', size=5)

    gs = GridSpec(4, 3, fig)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1:])
    ax3 = fig.add_subplot(gs[1, :])
    ax4 = fig.add_subplot(gs[2:, :])
    
    root_path = "data/duplex_data/"
    files = ["cbm2/CBM_2_rep1.masked.bed", "cbm2/CBM_2_rep2.masked.bed",
             "cbm3/CBM_3_rep1.sorted.bam.bed", "cbm3/CBM_3_rep2.sorted.bam.bed"]

    file_paths = [root_path + file for file in files]

    print("Loaded data. Calculating distances to CTCF summits")
    with concurrent.futures.ProcessPoolExecutor(len(file_paths)) as violin_executor:
        violin_futures = [violin_executor.submit(produce_violins, path, i+1) 
                          for i, path in enumerate(file_paths)]
        violin_all = pd.concat([future.result() for future in violin_futures], ignore_index=True)
    gc.collect()    
    distances = ctcf.DistDF(violin_all.groupby(["Pattern", "Distance", "abs"], observed=True)
                                      .agg({"Count" : sum})
                                      .reset_index())
    
    violin_reads = distances.explode_reads().reset_index(drop=True)
    gc.collect()  

    patterns = ["C:C", "5mC:5mC", "5hmC:5hmC", "C:5mC", "C:5hmC", "5mC:5hmC"]
    groups = [violin_reads['abs'][violin_reads['Pattern'] == pattern] for pattern in patterns]

    print("Kruskal-Wallis test:", stats.kruskal(*groups))
    dunn_result = sp.posthoc_dunn(violin_reads, val_col='abs', group_col='Pattern', p_adjust='holm')
    print("Dunn post hoc analysis", dunn_result)

    violin_reads.to_csv('source_data/fig5b_distances.csv.gz')

    sns.violinplot(violin_reads, 
                   x="Pattern", y="Distance", 
                   color="#99d8c9",
                   bw_method="silverman", bw_adjust=1.5,
                   cut=1,
                   linewidth=.8,
                   order=patterns,    
                   ax=ax3)
    
    print("Correlating patterns with distances")
    for pattern in patterns:
        x = np.where(violin_reads["Pattern"] == pattern, True, False)
        y = violin_reads["abs"]
        test = stats.pointbiserialr(x, y)
        print(pattern, test)    

        ax3.annotate(f"r$_{{pb}}$={round(test.statistic, 2)}", (pattern, 600), ha="center")    
    
    counts = violin_reads.groupby("Pattern").size()
    labels = [pattern + f"\nn={counts[pattern]}" for pattern in patterns]
    del violin_reads
    gc.collect()

    print("Loading data for barplots")
    with concurrent.futures.ProcessPoolExecutor(len(file_paths)) as load_executor:
        all_duplex_modbeds = [load_executor.submit(ctcf.read_merge, path, test_run, replicate+1) for replicate, path in enumerate(file_paths)]
        all_duplex_modbeds = [modbed.result() for modbed in all_duplex_modbeds]

    with concurrent.futures.ThreadPoolExecutor(len(all_duplex_modbeds)) as merge_executor:
        merged_motif_patterns = merge_executor.map(lambda pf: pf.merge_motif_patterns(min_depth), all_duplex_modbeds)

    with concurrent.futures.ThreadPoolExecutor(len(all_duplex_modbeds)) as merge_executor:
        merged_chip_patterns = merge_executor.map(lambda pf: pf.merge_chip_patterns(min_depth), all_duplex_modbeds)

    del all_duplex_modbeds
    gc.collect()

    ax3.axhline(0, ls=":", c="grey")
    ax3.set_xticks(patterns, labels)
    
    ax3.set_ylabel("Distance to summit (bp)")
    ax3.set_xlabel("Density of CpGs")

    with concurrent.futures.ThreadPoolExecutor(4) as pie_executor:
        motif_basecalls = pie_executor.map(lambda pf: pf.piechart_data(), merged_motif_patterns)
        motif_basecalls = pd.concat([ctcf.eval_basecall_proportions(pie.assign(Replicate = i)) for i, pie in enumerate(motif_basecalls)])

    with concurrent.futures.ThreadPoolExecutor(4) as pie_executor:
        chip_patterns = pie_executor.map(lambda pf: pf.piechart_data(), merged_chip_patterns)
        chip_basecalls = pd.concat([ctcf.eval_basecall_proportions(pie.assign(Replicate = i)) for i, pie in enumerate(chip_patterns)])

    gc.collect()
        
    root_path = "data/duplex_data/patterns/"
    files = ["CBM_2_rep1.masked.bed.duplex_patterns.tsv", "CBM_3_rep1.sorted.bam.bed.duplex_patterns.tsv",
             "CBM_2_rep2.masked.bed.duplex_patterns.tsv", "CBM_3_rep2.sorted.bam.bed.duplex_patterns.tsv"]

    file_paths = [root_path + file for file in files]

    genome_basecall_proportions = [(pd.read_table(path)
                        .rename(columns={"N_Pattern" : "Count"})
                        .replace(["-,m,C", "m,-,C", "-,h,C", "h,-,C", "m,h,C", "h,m,C", "-,-,C", "m,m,C", "h,h,C"], 
                                 ["C:5mC", "C:5mC", "C:5hmC", "C:5hmC", "5mC:5hmC", "5mC:5hmC", "C:C", "5mC:5mC", "5hmC:5hmC"])
                        .groupby(["Pattern"])["Count"]
                        .sum()
                        .reset_index()
                        .assign(Replicate = i))
                        for i, path in enumerate(file_paths)]
       
    genome_basecalls = pd.concat(map(ctcf.eval_basecall_proportions, genome_basecall_proportions))
    
    motif_basecalls["Type"], chip_basecalls["Type"], genome_basecalls["Type"] = "CTCF motif", "ChIP summit", "Genome average"
    motif_vs_chip = pd.concat([chip_basecalls, motif_basecalls, genome_basecalls]).reset_index()

    sns.barplot(motif_vs_chip,
                x="Pattern", y="Proportion", 
                hue="Type",
                palette="PuBuGn",
                err_kws={'linewidth': 0.8},
                capsize=.5,
                errorbar=("sd", 1),
                order=["C:C", "5mC:5mC"],
                hue_order=["Genome average", "CTCF motif", "ChIP summit"],
                ax=ax1)
    
    sns.stripplot(motif_vs_chip,
                  x="Pattern", y="Proportion", 
                  hue="Type", 
                  color="k", size=3,
                  legend=False, dodge=True, alpha=.66,
                  order=["C:C", "5mC:5mC"],
                  hue_order=["Genome average", "CTCF motif", "ChIP summit"],
                  ax=ax1)

    sns.move_legend(ax1, "upper center", frameon=False, title=None, ncols=3, bbox_to_anchor=(2, 1.2))
    ax1.get_legend().set_in_layout(False)
    ax1.set_ylabel("Proportion of duplex base-calls")   

    sns.barplot(motif_vs_chip,
            x="Pattern", y="Proportion", 
            hue="Type",
            palette="PuBuGn",
            err_kws={'linewidth': 0.8},
            capsize=.5,
            legend=False,
            errorbar=("sd", 1),
            order=["5hmC:5hmC", "C:5mC", "C:5hmC", "5mC:5hmC"],
            hue_order=["Genome average", "CTCF motif", "ChIP summit"],
            ax=ax2)
    
    sns.stripplot(motif_vs_chip,
                  x="Pattern", y="Proportion", 
                  color="k", size=3, hue="Type",
                  dodge=True, legend=False, alpha=.66,
                  order=["5hmC:5hmC", "C:5mC", "C:5hmC", "5mC:5hmC"],
                  hue_order=["Genome average", "CTCF motif", "ChIP summit"],
                  ax=ax2)

    with pd.ExcelWriter('source_data/fig5a_source_data.xlsx') as writer:
        motif_vs_chip.to_excel(writer, 'fig5a_proportions')
    
    ax2.set_ylabel(None)

    for ax in [ax1, ax2]:
        ax.set_xlabel(None) # type: ignore
    
    ax4.tick_params("both", bottom=False, labelbottom=False, left=False, labelleft=False)

    for i, ax in enumerate([ax1, ax3, ax4]):
        ax.set_title(string.ascii_lowercase[i], fontdict={"fontweight" : "bold"}, loc="left")        

    sns.despine()
    sns.despine(ax=ax4, left=True, bottom=True)
        
    if test_run:
        fig.savefig("plots/tests/ctcf_duplex_analysis.png")

    else:
        fig.savefig("plots/ctcf_duplex_analysis.png")
        fig.savefig("plots/ctcf_duplex_analysis.svg")

    # Stats for ax1 and ax2

    def comparer(exp_set, obs_set, pattern):
        obs = motif_vs_chip.groupby(["Type", "Pattern"]).get_group((obs_set, pattern)).reset_index()
        exp = motif_vs_chip.groupby(["Type", "Pattern"]).get_group((exp_set, pattern)).reset_index()

        exp_sum = exp["Count"].sum()
        test = pg.ttest(obs["Proportion"], exp["Proportion"], True)
        print(f"T-Test {obs_set} vs. {exp_set}:", test, "n1=", len(obs), "n2=", len(exp))       
        return 
    
    for pattern in patterns:    
        # print(pattern, "Genome vs. Genome sanity test :", comparer("Genome average", "Genome average", pattern))
        print(pattern, "Genome vs. Motif :", comparer("Genome average", "CTCF motif", pattern))
        print(pattern, "Genome vs. Summit :", comparer("Genome average", "ChIP summit", pattern))
        print(pattern, "Motif vs. Summit :", comparer("CTCF motif", "ChIP summit", pattern))

##### Main function #####

if __name__=="__main__":
    parser = argparse.ArgumentParser(
                        prog = "ctcf_duplex_analysis",
                        description = "Analyses duplex modified basecalls.")
    parser.add_argument("-t ", "--testrun", 
                        action="store_true", 
                        dest="test_run", 
                        required=False,
                        help="Whether a test output is produced.") 
    parser.add_argument("--min_depth", 
                        dest="min_depth", 
                        default=5,
                        required=False,
                        help="Whether to filter the dataset by depth at CpG.")
    parser.add_argument("--filter_distance", 
                        dest="filter_distance",
                        type=np.int64, 
                        default=500,
                        help="Whether to filter distance CpG dyads when comparing distance to summit.")


    global test_run, min_depth, filter_distance
    args = parser.parse_args()
    test_run, min_depth, filter_distance = args.test_run, args.min_depth, args.filter_distance
    main()    