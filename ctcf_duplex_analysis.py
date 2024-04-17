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
from scipy import stats

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

    def merge_ctcf_patterns(self, min_depth):
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
    
    def summarise_ctcf_duplex_states(self):
        ctcf_join = (pr.PyRanges(self.pattern_df, int64=True)
                          .join(chip_merge, apply_strand_suffix=True, suffix="_CTCF")
                          .as_df())
        
        print(f"Identified {len(ctcf_join)} overlaps with CTCF binding sites.")           

        ctcf_join = ctcf_join.melt(id_vars=["Chromosome", "Start", "End", "Strand_CTCF"], 
                                             var_name="Pattern", 
                                             value_vars=["CC", "MM", "HH", "MC", "CM", "HC", "CH", "HM", "MH"], 
                                             value_name="Count")
        
        assert ctcf_join.loc[:, ("Chromosome", "Start", "End")].duplicated().all() == False

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
        
        # identifies motif vs. opposite modifications
        summary["motif_mod"] = summary["pos_mod"].where(cond=(summary["Strand_CTCF"] == "+"), other=summary["neg_mod"])
        summary["opp_mod"] = summary["neg_mod"].where(cond=(summary["Strand_CTCF"] == "+"), other=summary["pos_mod"])

        summary = summary.groupby(["motif_mod", "opp_mod"]).sum(numeric_only=True).reset_index()
        summary = summary.replace(["C", "M", "H"], ["C", "5mC", "5hmC"])

        return  CTCF_Patterns(summary)
    
class CTCF_Patterns:
    def __init__(self, df):
        self.df = df

    def filter_hemi_reads(self):
        summary = self.df.query("motif_mod == 'C' | opp_mod == 'C'")
        total_count = summary["Count"].sum()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            summary["Percentage"] = (summary["Count"].div(total_count))*100
        return CTCF_Patterns(summary)

    def filter_hetero_reads(self):
        summary = self.df.query("(motif_mod == '5mC' & opp_mod == '5hmC') | (motif_mod == '5hmC' & opp_mod == '5mC')")
        total_count = summary["Count"].sum()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            summary["Percentage"] = (summary["Count"].div(total_count))*100
        return CTCF_Patterns(summary)

    def quantify_motif_strand(self):
        summary = self.df
        return summary.groupby("motif_mod").sum(numeric_only=True)
    
    def quantify_opp_strand(self):
        summary = self.df
        return summary.groupby("opp_mod").sum(numeric_only=True)
    
    def quantify_strands(self):
        motif = self.quantify_motif_strand().reset_index(names="Mod").assign(Strand = "Motif")
        opp = self.quantify_opp_strand().reset_index(names="Mod").assign(Strand = "Opposite")
        return pd.concat([motif, opp])
    
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
                                       int64=True).merge()

# removal of overlapping binding sites to prevent duplication downstream 
# binding site with least qValue chosen
with warnings.catch_warnings():
    warnings.simplefilter("ignore", FutureWarning)
    chip_merge = pr.PyRanges(ctcf_motif.join(ctcf_chip, apply_strand_suffix=False, suffix="_Peak")
                             .as_df()
                             .sort_values("qValue", ascending=True)
                             .groupby(["Chromosome", "Start", "End"])
                             .head(1)
                             .reset_index(drop=True),
                             int64=True)
    
def ctcf_intersect(df, **intersect_kwargs):
        # Loading JASPAR CTCF binding sites 
        intersect = pr.PyRanges(df, int64=True).intersect(chip_merge, strandedness=False, **intersect_kwargs).as_df()

        assert intersect.loc[:, ("Chromosome", "Start", "End")].duplicated().all() == False
        
        return intersect

@timer
def main(test_run=True, min_depth=5):
    root_path = "data/duplex_data/"
    files = ["cbm2/CBM_2_rep1.masked.bed", "cbm2/CBM_2_rep2.masked.bed",
             "cbm3/CBM_3_rep1.sorted.bam.bed", "cbm3/CBM_3_rep2.sorted.bam.bed"]

    file_paths = [root_path + file for file in files]

    with concurrent.futures.ProcessPoolExecutor(len(file_paths)) as load_executor:
        all_duplex_modbeds = [load_executor.submit(read_merge, path, test_run, replicate+1) for replicate, path in enumerate(file_paths)]
        all_duplex_modbeds = [modbed.result() for modbed in all_duplex_modbeds]
        merged_patterns = [pf.merge_ctcf_patterns(min_depth) for pf in all_duplex_modbeds]

    fig = plt.figure(figsize=(120/25.4, 89/25.4), dpi=600, layout="constrained")

    sns.set_style("ticks")
    mpl.rc('font', size=5)

    gs = GridSpec(2, 3, fig)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 2])

    sgs2 = gs[0, 1].subgridspec(3, 1)
    c_ax = fig.add_subplot(sgs2[0, 0])
    mc_ax = fig.add_subplot(sgs2[1, 0])
    hmc_ax = fig.add_subplot(sgs2[2, 0])

    ax3 = fig.add_subplot(gs[1, 0])
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
            radius=1, rotatelabels=False,
            labeldistance=1.25, pctdistance=.75,
            wedgeprops = {'linewidth': 0.1})
    
    with concurrent.futures.ThreadPoolExecutor(4) as stat_executor:
        mod_combinations = stat_executor.map(lambda pf: (pf.summarise_ctcf_duplex_states()
                                                           .quantify_mod_combinations()),
                                             merged_patterns)
    mod_combinations_summary = pd.concat([mod_combination.assign(Replicate = i) for i, mod_combination in enumerate(mod_combinations)])
    
    c, m, h = [mod_combinations_summary.groupby("motif_mod").get_group(mod) for mod in ["C", "5mC", "5hmC"]]

    sns.barplot(c, 
            y="motif_mod", x="Percentage",
            hue_order=["C", "5mC", "5hmC"],
            hue="opp_mod", palette=palette,
            errorbar=("sd", 1), err_kws={"lw" : .8}, capsize=.5,
            width=0.8, 
            legend=True,
            ax=c_ax)

    sns.move_legend(c_ax, loc="upper right", title=None, frameon=False)
    c_ax.get_legend().set_in_layout(False)
    # ax7.set_ylabel("Percent modification\nopposite motif")

    sns.barplot(m, 
                y="motif_mod", x="Percentage",
                hue="opp_mod", palette=palette,
                hue_order=["C", "5mC", "5hmC"],
                errorbar=("sd", 1), err_kws={"lw" : .8}, capsize=.5,
                width=0.8, 
                legend=False,
                ax=mc_ax)

    sns.barplot(h, 
                y="motif_mod", x="Percentage",
                hue="opp_mod", palette=palette,
                hue_order=["C", "5mC", "5hmC"],
                errorbar=("sd", 1), err_kws={"lw" : .8}, capsize=.5,
                width=0.8, 
                legend=False,
                ax=hmc_ax)
    
    for ax in [mc_ax, hmc_ax]:
        ax.sharex(c_ax)
    hmc_ax.set_xlabel("Percent opposite strand basecall")
    
    sns.despine()

    [ax.set_ylabel(None) for ax in [c_ax, mc_ax, hmc_ax]]

    for ax in [c_ax, mc_ax]:
        ax.set_xlabel(None)
        ax.tick_params("both", bottom=False, labelbottom=False)
        sns.despine(ax=ax, bottom=True)
    mc_ax.set_ylabel("Modification at motif")

    # Hemi # 

    with concurrent.futures.ThreadPoolExecutor(4) as stat_executor:
        hemi_mods_at_ctcf = stat_executor.map(lambda pf: (pf.extract_hemi_states()
                                                            .summarise_ctcf_duplex_states()
                                                            .filter_hemi_reads()
                                                            .quantify_strands()), 
                                              merged_patterns)

    hemi_mod_summary = pd.concat([hemi.assign(Replicate = i) for i, hemi in enumerate(hemi_mods_at_ctcf)])

    sns.barplot(hemi_mod_summary, 
        x="Mod", y="Percentage",
        order=["C", "5mC", "5hmC"],
        hue="Strand", palette=sns.color_palette("Paired", 4)[:2], dodge=True,
        errorbar=("sd", 1), err_kws={"lw" : .8}, capsize=.5,
        width=.8,
        ax=ax2)
    
    for mod in ["C", "5mC", "5hmC"]:
        mod_group = hemi_mod_summary.groupby("Mod").get_group(mod)
        motif = np.array(mod_group.query("Strand == 'Motif'")["Percentage"])
        opp = np.array(mod_group.query("Strand == 'Opposite'")["Percentage"])
        print(mod, stats.ttest_ind(motif, opp))
          
    ax2.set_ylabel("Percent of strand modification")
    ax2.set_ylim(0)
    sns.move_legend(ax2, "upper right", title="Strand", frameon=False)
    ax2.set_title(f"Hemi-dyads\nn={hemi_mod_summary['Count'].sum()} reads", loc="center")

    # Hetero # 

    with concurrent.futures.ThreadPoolExecutor(4) as stat_executor:
        hetero_mods_at_ctcf = stat_executor.map(lambda pf: (pf.extract_hetero_states()
                                                              .summarise_ctcf_duplex_states()
                                                              .filter_hetero_reads()
                                                              .quantify_strands()), 
                                                merged_patterns)

    hetero_mod_summary = pd.concat([hemi.assign(Replicate = i) for i, hemi in enumerate(hetero_mods_at_ctcf)]).reset_index(drop=True)

    sns.barplot(hetero_mod_summary, 
        x="Mod", y="Percentage",
        order=["5mC", "5hmC"],
        hue="Strand", palette=sns.color_palette("Paired", 4)[2:], dodge=True,
        errorbar=("sd", 1), err_kws={"lw" : .8}, capsize=.5,
        width=.8,
        ax=ax3)
    
    for mod in ["5mC", "5hmC"]:
        mod_group = hetero_mod_summary.groupby("Mod").get_group(mod)
        motif = np.array(mod_group.query("Strand == 'Motif'")["Percentage"])
        opp = np.array(mod_group.query("Strand == 'Opposite'")["Percentage"])
        print(mod, stats.ttest_ind(motif, opp))
    
    ax3.set_ylabel("Percent of strand basecalls")
    ax3.set_ylim(0, 70)
    sns.move_legend(ax3, "upper left", title="Strand", frameon=False, ncol=2)
    ax3.set_title(f"Hetero-dyads\nn={hetero_mod_summary['Count'].sum()} reads", loc="center")

    for ax in [ax2, ax3]:
        ax.set_xlabel(None)

    for i, ax in enumerate([ax1, c_ax, ax2, ax3, ax4]):
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

    ax4.tick_params("both", bottom=False, labelbottom=False, left=False, labelleft=False)
    sns.despine(ax=ax4, left=True, bottom=True)

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

    with pd.ExcelWriter("data/ctcf/ctcf_intersects.xlsx", mode="w") as writer:
        for i, pattern in enumerate(merged_patterns):
            ((pr.PyRanges(pattern.pattern_df, int64=True)
                .join(chip_merge, apply_strand_suffix=True, suffix="_CTCF"))
                .as_df()
                .sort_values("qValue", ascending=False)
                .groupby(["Chromosome", "Start", "End"], observed=True)
                .head(1)
                .to_excel(writer, sheet_name=f"Sheet {i+1}", index=False))

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

    args = parser.parse_args()

    if args.test_run:
        test_run = True
    else: 
        test_run = False

    args = parser.parse_args()
    main(test_run, args.min_depth)    