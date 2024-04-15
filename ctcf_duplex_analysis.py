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
import upsetplot
import warnings

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
                .eval("readCount = CC + MM + HH + CM + MC + HC + CH + MH + HM"))
        return PatternFrame(df)

def read_merge(path, test_run=True, replicate=None):
    return merge_pattern(read_duplex_modbed(path, test_run, replicate))

class PatternFrame:
    def __init__(self, pattern_df):
        self.pattern_df = pattern_df
        self.ctcf_patterns = ctcf_intersect(pattern_df)

    def merge_patterns(self, min_depth):
        df = (self.ctcf_patterns
              .query(f"readCount > {min_depth}")
              .eval("Hetero = MH + HM")
              .eval("Hemi = CM + MC + CH + HC"))
        
        df = df.assign(Majority = lambda r: r.loc[:, ("CC", "MM", "HH", "Hetero", "Hemi")].idxmax(axis=1))

        return MergedSites(df, min_depth)
           
    def piechart_data(self):
        """
        Outputs a dataframe counting all homo-modification states. Hetero- and hemi-modification states are grouped. 
        """
        df = self.ctcf_patterns
        
        pie_data = pd.DataFrame({
            "Pattern" : ["Homo-C", "Homo-5mC", "Homo-5hmC", "Hemi-dyad", "Hetero-dyad"],
            "Count" : [df["CC"].sum(), df["MM"].sum(), df["HH"].sum(), 
                    df["MC"].sum() + df["CM"].sum() + df["HC"].sum() + df["CH"].sum(),
                    df["MH"].sum() + df["HM"].sum()]
            })

        return pie_data

    def quantify_hemihetero_sites(self, min_depth):
        df = self.ctcf_patterns

        # need strand information from CTCF sites
        ctcf_join = (pr.PyRanges(df, int64=True)
                          .join(chip_merge, apply_strand_suffix=True, suffix="_CTCF")
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
    
class MergedSites(PatternFrame):
    def __init__(self, pattern_df, min_depth):
        super().__init__(pattern_df)
        self.min_depth = min_depth

    def piechart_data(self):
        """
        Outputs a dataframe counting all homo-modification states. Hetero- and hemi-modification states are grouped. 
        """
        df = self.pattern_df
        pie_data = df["Majority"].value_counts(normalize=True)

        return pie_data
    
    def extract_hemi_states(self):
        hemi = self.pattern_df.query("Majority == 'Hemi'")
        return MergedSites(hemi, self.min_depth)

    def extract_hetero_states(self):
        hetero = self.pattern_df.query("Majority == 'Hetero'")
        return MergedSites(hetero, self.min_depth)
    
class HemiHeteroPatterns:
    def __init__(self, df):
        self.df = df

    def hemi(self):
        summary = self.df.query("motif_mod == 'C' | opp_mod == 'C'")
        total_count = summary["Count"].sum()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            summary["Percentage"] = (summary["Count"].div(total_count))*100
        return HemiHeteroPatterns(summary)

    def hetero(self):
        summary = self.df.query("(motif_mod == '5mC' & opp_mod == '5hmC') | (motif_mod == '5hmC' & opp_mod == '5mC')")
        total_count = summary["Count"].sum()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            summary["Percentage"] = (summary["Count"].div(total_count))*100
        return HemiHeteroPatterns(summary)

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
        
ctcf_chip = pr.PyRanges(pd.read_table("data/ctcf/ChIP2MACS2/MACS2/ENCSR000CBN_peaks.narrowPeak", 
                                      names=["Chromosome", "Start", "End", "Name", "Pileup", 
                                             "Strand", "FoldDifference", "pValue", "qValue", 
                                             "Peak"],
                                      usecols=["Chromosome", "Start", "End", "FoldDifference", 
                                               "pValue", "qValue"]), 
                                      int64=True)
ctcf_motif = pr.PyRanges(pd.read_table("feature_references/CTCF_mm39_jaspar_sorted.tsv",
                                       names=["Chromosome", "Start", "End", "Width", "Strand",
                                              "Name", "Score", "pValue", "qValue", "Sequence"],
                                       usecols=["Chromosome", "Start", "End", "Strand"],
                                       skiprows=1), 
                                       int64=True)
with warnings.catch_warnings():
    warnings.simplefilter("ignore", FutureWarning)
    chip_merge = ctcf_chip.join(ctcf_motif, how="left", 
                                apply_strand_suffix=False, suffix="_Motif")
    
def ctcf_intersect(df, **intersect_kwargs):
        # Loading JASPAR CTCF binding sites 
        return pr.PyRanges(df, int64=True).intersect(chip_merge, strandedness=False, **intersect_kwargs).as_df()

@timer
def main(test_run=True, min_depth=5, merge_sites=False):
    root_path = "data/duplex_data/"
    files = ["cbm2/CBM_2_rep1.masked.bed", "cbm2/CBM_2_rep2.masked.bed",
             "cbm3/CBM_3_rep1.sorted.bam.bed", "cbm3/CBM_3_rep2.sorted.bam.bed"]

    file_paths = [root_path + file for file in files]

    with concurrent.futures.ProcessPoolExecutor(len(file_paths)) as load_executor:
        all_duplex_modbeds = [load_executor.submit(read_merge, path, test_run, replicate+1) for replicate, path in enumerate(file_paths)]
        all_duplex_modbeds = [modbed.result() for modbed in all_duplex_modbeds]
        merged_patterns = [pf.merge_patterns(5) for pf in all_duplex_modbeds]
                              
    fig = plt.figure(figsize=(120/25.4, 89/25.4), dpi=600, layout="constrained")

    sns.set_style("ticks")
    mpl.rc('font', size=5)

    gs = GridSpec(2, 3, fig)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])

    sgs2 = gs[1, 0].subgridspec(3, 1)
    c_ax = fig.add_subplot(sgs2[0, 0])
    mc_ax = fig.add_subplot(sgs2[1, 0])
    hmc_ax = fig.add_subplot(sgs2[2, 0])

    ax4 = fig.add_subplot(gs[1, 1:])

    palette = {
    "C" : "#e0f3db", 
    "5mC" : "#a8ddb5", 
    "5hmC" : "#43a2ca"
           }
        
    pie_data = pd.concat([merged_pattern.piechart_data().reset_index(name="Proportion") for merged_pattern in merged_patterns])
    pie_data = (pie_data.groupby("index")
                .mean(numeric_only=True)
                .reset_index()
                .replace(["CC", "MM", "HH", "Hetero", "Hemi"], 
                         ["Homo-C", "Homo-M", "Homo-H", "Hetero-dyad", "Hemi-dyad"]))
    
    ax1.pie(pie_data["Proportion"], 
            labels=pie_data["index"], 
            colors=sns.color_palette("GnBu"),
            startangle=45, counterclock=False,
            autopct='%.1f',
            radius=0.75, rotatelabels=False,
            labeldistance=1.25, pctdistance=.75, 
            explode=(0, .15, .15, .15, .15))

    all_duplex_modbeds[0].merge_patterns(min_depth).extract_hemi_states().quantify_hemihetero_sites(min_depth)
    exit()

    with concurrent.futures.ThreadPoolExecutor(4) as stat_executor:
        hemihetero_mods_at_ctcf = [stat_executor.submit(pf.quantify_hemihetero_sites, merge_sites, min_depth) for pf in all_duplex_modbeds]
        hemihetero_mods_at_ctcf = [summary.result() for summary in hemihetero_mods_at_ctcf]

        # Hemi 

        motifs_only = [summary.hemi().quantify_motif_strand() for summary in hemihetero_mods_at_ctcf]
        motif_mod_summary = pd.concat([summary.assign(Replicate = i, 
                                                      Kind = "Motif") for i, summary in enumerate(motifs_only)])
        
        opp_only = [summary.hemi().quantify_opp_strand() for summary in hemihetero_mods_at_ctcf]
        opp_mod_summary = pd.concat([summary.assign(Replicate = i, 
                                                    Kind = "Opposite") for i, summary in enumerate(opp_only)])
        
        mod_combinations = pd.concat([summary.quantify_mod_combinations()
                                        .assign(Replicate = i) for i, summary in enumerate(hemihetero_mods_at_ctcf)]).reset_index()

    hemi_mod_sum = pd.concat([motif_mod_summary, opp_mod_summary]).reset_index(names="Mod")

    sns.barplot(hemi_mod_sum, 
        x="Mod", y="Percentage",
        order=["C", "5mC", "5hmC"],
        hue="Kind", palette="Paired", dodge=True,
        errorbar=("sd", 1), err_kws={"lw" : .8}, capsize=.5,
        width=.8,
        ax=ax2)
        
    ax2.set_ylabel("Percent of strand basecalls")
    ax2.set_ylim(0, 70)
    sns.move_legend(ax2, "upper left", title="Strand", frameon=False, ncol=2)

        # Hetero
    
    motifs_only = [summary.hetero().quantify_motif_strand() for summary in hemihetero_mods_at_ctcf]
    motif_mod_summary = pd.concat([summary.assign(Replicate = i, 
                                                    Kind = "Motif") for i, summary in enumerate(motifs_only)])
    
    opp_only = [summary.hetero().quantify_opp_strand() for summary in hemihetero_mods_at_ctcf]
    opp_mod_summary = pd.concat([summary.assign(Replicate = i, 
                                                Kind = "Opposite") for i, summary in enumerate(opp_only)])
    
    mod_combinations = pd.concat([summary.quantify_mod_combinations()
                                    .assign(Replicate = i) for i, summary in enumerate(hemihetero_mods_at_ctcf)]).reset_index()

    hetero_mod_sum = pd.concat([motif_mod_summary, opp_mod_summary]).reset_index(names="Mod")

    sns.barplot(hetero_mod_sum, 
        x="Mod", y="Percentage",
        order=["5mC", "5hmC"],
        hue="Kind", palette="Paired", dodge=True,
        errorbar=("sd", 1), err_kws={"lw" : .8}, capsize=.5,
        width=.8,
        ax=ax3)
        
    ax3.set_ylabel("Percent of strand basecalls")
    ax3.set_ylim(0, 70)
    sns.move_legend(ax3, "upper left", title="Strand", frameon=False, ncol=2)

    c, m, h = [mod_combinations.groupby("motif_mod").get_group(mod) for mod in ["C", "5mC", "5hmC"]]

    sns.barplot(c, 
            y="motif_mod", x="Percentage",
            hue_order=["5mC", "5hmC"],
            hue="opp_mod", palette=palette,
            errorbar=("sd", 1), err_kws={"lw" : .8}, capsize=.5,
            width=0.8, 
            legend=False,
            ax=c_ax)
    
    c_ax.sharex(mc_ax)
    # ax7.set_ylabel("Percent modification\nopposite motif")

    sns.barplot(m, 
                y="motif_mod", x="Percentage",
                hue="opp_mod", palette=palette,
                errorbar=("sd", 1), err_kws={"lw" : .8}, capsize=.5,
                width=0.8, 
                legend=False,
                ax=mc_ax)

    sns.barplot(h, 
                y="motif_mod", x="Percentage",
                hue="opp_mod", palette=palette,
                errorbar=("sd", 1), err_kws={"lw" : .8}, capsize=.5,
                width=0.8, 
                legend=False,
                ax=hmc_ax)
    hmc_ax.sharex(mc_ax)
    hmc_ax.set_xlabel("Percent opposite strand basecall")
    
    sns.despine()

    [ax.set_ylabel(None) for ax in [c_ax, mc_ax, hmc_ax]]

    for ax in [c_ax, mc_ax]:
        ax.set_xlabel(None)
        ax.tick_params("both", bottom=False, labelbottom=False)
        sns.despine(ax=ax, bottom=True)

    for i, ax in enumerate([ax1, ax2, ax3, c_ax, ax4]):
        ax.set_title(string.ascii_lowercase[i], fontdict={"fontweight" : "bold"}, loc="left")        
        
    # bbox = ax3.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    # uw, uh = bbox.width, bbox.height

    ctcf_concat = pd.concat([patternframe.ctcf_patterns for patternframe in all_duplex_modbeds])
    ctcf_upset_mshp = upsetplot.from_memberships(memberships=[
        ["C"], 
        ["5mC"],
        ["C", "5mC"],
        ["5hmC"],
        ["5mC", "5hmC"],
        ["C", "5hmC"]],
        data = [
            ctcf_concat["CC"].sum(), 
            ctcf_concat["MM"].sum(), 
            ctcf_concat["CM"].sum() + ctcf_concat["MC"].sum(), 
            ctcf_concat["HH"].sum(), 
            ctcf_concat["HM"].sum() + ctcf_concat["MH"].sum(), 
            ctcf_concat["CH"].sum() + ctcf_concat["HC"].sum()])
    
    upset_fig = plt.figure(dpi=600)

    ctcf_upset = upsetplot.UpSet(ctcf_upset_mshp, sort_by="cardinality", sort_categories_by="input", show_percentages=True, 
                                 facecolor="grey", 
                                 element_size=None, 
                                 intersection_plot_elements=3)

    ctcf_upset.plot(upset_fig)
    
    if test_run:
        fig.savefig("plots/tests/ctcf_duplex_analysis.png")
        upset_fig.savefig("plots/tests/ctcf_upset.png")

    else:
        fig.savefig("plots/ctcf_duplex_analysis.png")
        fig.savefig("plots/ctcf_duplex_analysis.svg")
        upset_fig.savefig("plots/ctcf_upset.svg")
        upset_fig.savefig("plots/ctcf_upset.png")

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