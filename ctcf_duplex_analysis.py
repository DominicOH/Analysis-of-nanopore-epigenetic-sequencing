import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
import concurrent.futures
import pyranges as pr
import argparse
from AnalysisTools.helpers import timer
import string
import numpy as np

def read_duplex_modbed(path, test_run=True, replicate=None):
    if test_run:
        nrows = 1000000
    else:
        nrows=None
    
    pattern_df = (pd.read_table(path, sep="\t", 
                            names=["Chromosome", "Start", "End", "Pattern", "readCount", "D0", "D1", "D2", "D3", "D4", "percentPattern", "N_Pattern", "N_Canonical", "N_Other", "D5", "D6", "D7", "D8"],
                            usecols=["Chromosome", "Start", "End", "Pattern", "readCount", "N_Pattern"],
                            nrows=nrows)
                            .pivot_table(values="N_Pattern", columns="Pattern", index=["Chromosome", "Start", "End"], fill_value=0)
                            .reset_index())
    print(f"Found {len(pattern_df)} sites in {path}")

    if replicate:
        pattern_df["Replicate"] = replicate

    return pattern_df

def merge_pattern(df):
        df = (df.rename(columns={"-,-,C" : "CC",
                                 "m,m,C" : "MM",
                                 "h,h,C" : "HH",
                                 "-,m,C" : "CM",
                                 "m,-,C" : "MC",
                                 "h,-,C" : "HC",
                                 "-,h,C" : "CH",
                                 "m,h,C" : "MH",
                                 "h,m,C" : "HM"})
                .eval("total_homometh = CC + MM + HH")
                .eval("total_heterometh = CM + MC + CH + HC + MH + HM")
                .eval("readCount = total_homometh + total_heterometh"))
        return PatternFrame(df)

def read_merge(path, test_run=True, replicate=None):
    return merge_pattern(read_duplex_modbed(path, test_run, replicate))

class PatternFrame:
    def __init__(self, pattern_df):
        self.pattern_df = pattern_df
        self.ctcf_patterns = ctcf_intersect(pattern_df)
    
    def quantify_homomodification_type(self):
        """
        Quantifies modification type for homomodified CTCF-overlapping CpG sites. 
        """
        counts_df = self.ctcf_patterns
        counts_df = counts_df.loc[:,("CC", "HH", "MM")].sum().reset_index(name="Count")
        counts_df["Percentage"] = (counts_df["Count"]/counts_df["Count"].sum()).multiply(100)
        counts_df = counts_df.replace(["CC", "MM", "HH"], ["C", "5mC", "5hmC"])        
        return counts_df
        
    def quantify_homomodified_vs_heteromodified(self):
        """
        Finds the proportion of CTCF-bind overlapping CpG basecalls that are homo/heteromodified. 
        """
        df = self.ctcf_patterns
        summary = ((df[["total_homometh", "total_heterometh"]].sum()
                   .div(df["readCount"].sum()).multiply(100))
                   .reset_index(name="Percentage"))
        
        return summary
    
    def piechart_data(self):
        """
        Outputs a dataframe counting all homo-modification states. Hetero- and hemi-modification states are grouped. 
        """
        df = self.ctcf_patterns
        
        pie_data = pd.DataFrame({
            "Pattern" : ["Homo-C", "Homo-5mC", "Homo-5hmC", "Hemi-", "Hetero-"],
            "Count" : [df["CC"].sum(), df["MM"].sum(), df["HH"].sum(), 
                    df["MC"].sum() + df["CM"].sum() + df["HC"].sum() + df["CH"].sum(),
                    df["MH"].sum() + df["HM"].sum()]
            })

        return pie_data
    
    def extract_heteromodified_site(self, min_depth=1):
        """
        Finds the proportion of CTCF-bind overlapping CpG basecalls that are homo/heteromodified. 
        """
        df = self.ctcf_patterns
        df = df.query(f"readCount >= {min_depth}")
        hetero = df.query("total_heterometh > total_homometh")
        
        return hetero
    
    def quantify_homo_hetero_site_proportion(self, min_depth=1):
        """
        Finds the proportion of CTCF-bind overlapping CpG sites that are homo/heteromodified. 
        """
        df = self.ctcf_patterns
        df = df.query(f"readCount >= {min_depth}")
        hetero = len(df.query("total_heterometh > total_homometh"))
        result = pd.DataFrame({"Pattern" : ["total_heterometh", "total_homometh"],
                               "Percentage" : [(hetero/len(df))*100, ((len(df) - hetero)/len(df))*100]})
        
        return result

    def hemihetero_sites(self, min_depth):
        """
        Restricts the dataset to just sites that are hetero or hemi modified in the majority of reads. 
        """
        df = self.ctcf_patterns.query(f"readCount >= {min_depth}")
        return df.loc[df.eval("total_heterometh > total_homometh")]

    def quantify_hemihetero_sites(self, merge_sites, min_depth):
        if merge_sites:
            df = self.hemihetero_sites(min_depth)
        else:
            df = self.ctcf_patterns
            df = df.loc[df.eval("(MC + CM + HC + CH + HM + MH) > 0")]
        df = df.loc[:, ("Chromosome", "Start", "End", "MC", "CM", "HC", "CH", "HM", "MH")]

        # need strand information from CTCF sites
        ctcf_join = (pr.PyRanges(df)
                          .join(ctcf_bs, apply_strand_suffix=True, suffix="_CTCF")
                          .as_df())
        print(f"Identified {len(ctcf_join)} overlaps with CTCF binding sites.")           

        # remove those with no hemimethlyation
        ctcf_join = ctcf_join.melt(id_vars=["Chromosome", "Start", "End", "Strand_CTCF"], 
                                             var_name="Pattern", 
                                             value_vars=["MC", "CM", "HC", "CH", "HM", "MH"], 
                                             value_name="Count")
        ctcf_join = ctcf_join.loc[ctcf_join.eval("Count > 0")]

        summary = (ctcf_join.groupby(["Strand_CTCF", "Pattern"])["Count"]
                   .sum()
                   .reset_index())
        
        summary["Percentage"] = (summary["Count"]
                                 .div(summary["Count"].sum())
                                 .multiply(100))

        summary = pd.concat([summary, (summary.Pattern.str.split("", n=2, expand=True)
                                .drop(columns=0)
                                .rename(columns={1 : "pos_mod", 
                                                 2 : "neg_mod"}))], axis=1)

        summary["motif_mod"] = summary["pos_mod"].where(cond=(summary["Strand_CTCF"] == "+"), other=summary["neg_mod"])
        summary["opp_mod"] = summary["neg_mod"].where(cond=(summary["Strand_CTCF"] == "+"), other=summary["pos_mod"])

        summary = summary.groupby(["motif_mod", "opp_mod"]).sum(numeric_only=True).reset_index()
        summary = summary.replace(["C", "M", "H"], ["C", "5mC", "5hmC"])

        return  HemiHeteroPatterns(summary)
    
class HemiHeteroPatterns:
    def __init__(self, df):
        self.df = df

    def quantify_motif_strand(self):
        summary = self.df
        return summary.groupby("motif_mod").sum(numeric_only=True)
    
    def quantify_opp_strand(self):
        summary = self.df
        return summary.groupby("opp_mod").sum(numeric_only=True)
    
    def quantify_mod_combinations(self):
        summary = self.df.drop(columns="Percentage")
        motif_total = summary.groupby("motif_mod")["Count"].transform("sum")

        summary["Percentage"] = (summary["Count"]/motif_total)*100

        return  summary
    
ctcf_bs = pr.PyRanges(pd.read_table("feature_references/CTCF_mm39_jaspar_sorted.tsv")
                          .query(f"qvalue < 0.05")
                          .rename(columns={"seqnames" : "Chromosome", "start" : "Start", "end" : "End", "strand" : "Strand"}))

def ctcf_intersect(df, **intersect_kwargs):
        # Loading JASPAR CTCF binding sites 
        return pr.PyRanges(df).intersect(ctcf_bs, strandedness=False, **intersect_kwargs).as_df()

@timer
def main(test_run=True, min_depth=1, merge_sites=False):
    root_path = "data/duplex_data/"
    files = ["cbm2/CBM_2_rep1.masked.bed", "cbm2/CBM_2_rep2.masked.bed",
             "cbm3/CBM_3_rep1.sorted.bam.bed", "cbm3/CBM_3_rep2.sorted.bam.bed"]

    file_paths = [root_path + file for file in files]

    
    with concurrent.futures.ProcessPoolExecutor(len(file_paths)) as load_executor:
        all_duplex_modbeds = [load_executor.submit(read_merge, path, test_run, replicate+1) for replicate, path in enumerate(file_paths)]
        all_duplex_modbeds = [modbed.result() for modbed in all_duplex_modbeds]
                              
    fig = plt.figure(figsize=(89/25.4, 100/25.4), dpi=600, layout="constrained")

    sns.set_style("ticks")
    mpl.rc('font', size=5)

    gs = GridSpec(3, 4, fig)

    ax1 = fig.add_subplot(gs[0, :2])
    ax2 = fig.add_subplot(gs[0, 2:]) 
    ax3 = fig.add_subplot(gs[1, :2])
    ax4 = fig.add_subplot(gs[1, 2:])

    sgs = gs[2, :2].subgridspec(1, 3)

    ax7 = fig.add_subplot(sgs[0, 0])
    ax8 = fig.add_subplot(sgs[0, 1])
    ax9 = fig.add_subplot(sgs[0, 2])

    palette = {
    "C" : "#e0f3db", 
    "5mC" : "#a8ddb5", 
    "5hmC" : "#43a2ca"
           }
    
    with concurrent.futures.ThreadPoolExecutor(4) as stat_executor:
        homohemi_proportions = [stat_executor.submit(pf.piechart_data) for pf in all_duplex_modbeds]
        homohemi_proportions = pd.concat([(proportion.result()) for proportion in homohemi_proportions]).reset_index()

    pie_data = homohemi_proportions.groupby("Pattern")["Count"].sum().reset_index()
    homohemi_proportions = homohemi_proportions.replace(["total_homometh", "total_heterometh"], 
                                                        ["Homo-\nmodified", "Hetero-\nmodified"])
    
    ax2.pie(pie_data["Count"], 
            labels=pie_data["Pattern"], 
            colors=sns.color_palette("YlGn_r"),
            startangle=45, counterclock=False,
            autopct='%.1f',
            radius=0.75, rotatelabels=False,
            labeldistance=1.25, pctdistance=.75, 
            explode=(0, 0, 0, 0, .15))

    with concurrent.futures.ThreadPoolExecutor(4) as stat_executor:
        hetero_mods_at_ctcf = [stat_executor.submit(pf.quantify_hemihetero_sites, merge_sites, min_depth) for pf in all_duplex_modbeds]
        hetero_mods_at_ctcf = [summary.result() for summary in hetero_mods_at_ctcf]

        motifs_only = [summary.quantify_motif_strand() for summary in hetero_mods_at_ctcf]
        motif_mod_summary = pd.concat([summary.assign(Replicate = i, 
                                                      Kind = "Motif") for i, summary in enumerate(motifs_only)])

        opp_only = [summary.quantify_opp_strand() for summary in hetero_mods_at_ctcf]
        opp_mod_summary = pd.concat([summary.assign(Replicate = i, 
                                                    Kind = "Opposite") for i, summary in enumerate(opp_only)])
        
        mod_combinations = pd.concat([summary.quantify_mod_combinations()
                                        .assign(Replicate = i) for i, summary in enumerate(hetero_mods_at_ctcf)]).reset_index()

    hetero_mod_sum = pd.concat([motif_mod_summary, opp_mod_summary]).reset_index(names="Mod")

    sns.barplot(hetero_mod_sum, 
        x="Mod", y="Percentage",
        order=["C", "5mC", "5hmC"],
        hue="Kind", palette="BuGn_r", dodge=True,
        errorbar=("sd", 1), err_kws={"lw" : .8}, capsize=.5,
        width=.8,
        ax=ax4)
        
    ax4.set_ylabel("Percent of strand basecalls")
    ax4.set_ylim(0, 100)
    sns.move_legend(ax4, "upper left", title="Relation to motif strand", frameon=True, ncol=2)
    c, m, h = [mod_combinations.groupby("motif_mod").get_group(mod) for mod in ["C", "5mC", "5hmC"]]

    sns.barplot(c, 
            x="motif_mod", y="Percentage",
            hue_order=["5mC", "5hmC"],
            hue="opp_mod", palette=palette,
            errorbar=("sd", 1), err_kws={"lw" : .8}, capsize=.5,
            width=0.8, 
            legend=False,
            ax=ax7)
    
    ax7.sharey(ax8)
    ax7.set_ylabel("Percent modification\nopposite motif")

    sns.barplot(m, 
                x="motif_mod", y="Percentage",
                hue="opp_mod", palette=palette,
                errorbar=("sd", 1), err_kws={"lw" : .8}, capsize=.5,
                width=0.8, 
                legend=False,
                ax=ax8)

    sns.barplot(h, 
                x="motif_mod", y="Percentage",
                hue="opp_mod", palette=palette,
                errorbar=("sd", 1), err_kws={"lw" : .8}, capsize=.5,
                width=0.8, 
                legend=False,
                ax=ax9)
    ax9.sharey(ax8)
    
    sns.despine()

    [ax.set_xlabel(None) for ax in [ax2, ax4, ax7, ax8, ax9]]
    for ax in [ax8, ax9]:
        ax.set_ylabel(None)
        ax.tick_params("both", left=False, labelleft=False)
        sns.despine(ax=ax, left=True)

    for i, ax in enumerate(fig.axes[:(len(fig.axes)-2)]):
        ax.set_title(string.ascii_lowercase[i], fontdict={"fontweight" : "bold"}, loc="left")

    for ax in [ax1, ax3]:
        sns.despine(ax=ax, bottom=True, left=True)
        ax.tick_params("both", bottom=False, left=False, 
                        labelbottom=False, labelleft=False)
    
    if test_run:
        fig.savefig("plots/tests/ctcf_duplex_analysis.png")
    else:
        fig.savefig("plots/ctcf_duplex_analysis.png")
        fig.savefig("plots/ctcf_duplex_analysis.svg")

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
                        default=1,
                        required=False,
                        help="Whether to filter the dataset by depth at CpG.")
    parser.add_argument("--merge_sites", 
                        default=False,
                        action="store_true",
                        dest="merge_sites", 
                        help="Whether to reads from overlapping CpG sites.")

    args = parser.parse_args()

    if args.test_run:
        test_run = True
    else: 
        test_run = False

    args = parser.parse_args()
    main(test_run, args.min_depth, args.merge_sites)    