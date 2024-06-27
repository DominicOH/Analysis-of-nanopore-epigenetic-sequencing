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
from scipy.stats import contingency
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
    def __init__(self, pattern_df: pd.DataFrame):
        self.pattern_df = pattern_df
        self.motif_patterns = motif_join(pattern_df)
        self.chip_patterns = chip_join(pattern_df)
            
    def merge_chip_patterns(self, min_depth):
        df = (self.chip_patterns
              .query(f"readCount > {min_depth}"))
        
        df = df.melt(id_vars=["Chromosome", "Start", "End", "Strand_ChIP"], 
                     value_vars=["CC", "MM", "HH", "CM" , "MC" , "HC" , "CH" , "MH" , "HM"],
                     var_name="Pattern", value_name="Count")

        return MeltDF(df, min_depth)
    
    def merge_motif_patterns(self, min_depth):
        df = (self.motif_patterns
              .query(f"readCount > {min_depth}"))
        
        df = df.melt(id_vars=["Chromosome", "Start", "End", "Strand_Motif"], 
                     value_vars=["CC", "MM", "HH", "CM" , "MC" , "HC" , "CH" , "MH" , "HM"],
                     var_name="Pattern", value_name="Count")

        return MeltDF(df, min_depth)
           
class MeltDF(PatternFrame):
    def __init__(self, pattern_df, min_depth):
        super().__init__(pattern_df)
        self.min_depth = min_depth

    def piechart_data(self):
        """
        Outputs a dataframe counting all constitutive-modification states. Hetero-modification and hemi-methylation states are grouped. 
        """
        df = self.pattern_df
        df = df.replace(["CM", "MC", "CH", "HC", "MH", "HM", "CC", "MM", "HH"], 
                        ["C:5mC", "C:5mC", "C:5hmC", "C:5hmC", "5mC:5hmC", "5mC:5hmC", "C:C", "5mC:5mC", "5hmC:5hmC"])
        pie_data = df.groupby("Pattern")["Count"].sum().reset_index()

        return pie_data
                
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

print(len(ctcf_motif), "CTCF motifs")

ctcf_summits = pr.PyRanges(pd.read_table("data/ctcf/ChIP2MACS2/MACS2/ENCSR000CBN_summits.bed", 
                                         names=["Chromosome", "Start", "End", "Name", "Score"]),
                           int64=True).merge()

print(len(ctcf_summits), "CTCF summits")

# removal of overlapping binding sites to prevent duplication downstream 
# binding site with least qValue chosen
with warnings.catch_warnings():
    warnings.simplefilter("ignore", FutureWarning)
    # will find the intersect of sites that are present on CTCF motifs, 
    # that overlap with summits from confident peaks
    chip_merge = pr.PyRanges(ctcf_motif
                             .join(ctcf_summits, apply_strand_suffix=False, suffix="_Summit")
                             .join(ctcf_chip, apply_strand_suffix=False, suffix="_Peak")
                             .as_df()
                             .sort_values("qValue", ascending=True)
                             .groupby(["Chromosome", "Start", "End"])
                             .head(1)
                             .reset_index(drop=True),
                             int64=True)
    
print(len(chip_merge), "CTCF motifs that overlap summits and peaks")

def motif_join(df: pd.DataFrame) -> pd.DataFrame:
    """
    Intersects with CTCF motifs - not necessarily those within ChIP summit peaks. 
    Overlapping motifs are merged. 
    """
    # Loading JASPAR CTCF binding sites 
    intersect = pr.PyRanges(df, int64=True).join(ctcf_motif, strandedness=False, apply_strand_suffix=True, suffix="_Motif").as_df()

    assert intersect.loc[:, ("Chromosome", "Start", "End")].duplicated().all() == False # type: ignore
        
    return intersect
    
def chip_join(df: pd.DataFrame) -> pd.DataFrame:
    """
    Intersects with CTCF motifs present in ChIP summit peaks. 
    Note that this joins with CTCF summits to remove adjacent motifs that may not be bound.
    """
    # Loading JASPAR CTCF binding sites 
    intersect = pr.PyRanges(df, int64=True).join(chip_merge, strandedness=False, apply_strand_suffix=True, suffix="_ChIP").as_df()

    assert intersect.loc[:, ("Chromosome", "Start", "End")].duplicated().all() == False # type: ignore
        
    return intersect

@timer
def main(test_run=True, min_depth=10, upset_plot=False):

    root_path = "data/duplex_data/"
    files = ["cbm2/CBM_2_rep1.masked.bed", "cbm2/CBM_2_rep2.masked.bed",
             "cbm3/CBM_3_rep1.sorted.bam.bed", "cbm3/CBM_3_rep2.sorted.bam.bed"]

    file_paths = [root_path + file for file in files]

    with concurrent.futures.ProcessPoolExecutor(len(file_paths)) as load_executor:
        all_duplex_modbeds = [load_executor.submit(read_merge, path, test_run, replicate+1) for replicate, path in enumerate(file_paths)]
        all_duplex_modbeds = [modbed.result() for modbed in all_duplex_modbeds]
        
    with concurrent.futures.ThreadPoolExecutor(len(all_duplex_modbeds)) as merge_executor:
        merged_motif_patterns = merge_executor.map(lambda pf: pf.merge_motif_patterns(min_depth), all_duplex_modbeds)

    with concurrent.futures.ThreadPoolExecutor(len(all_duplex_modbeds)) as merge_executor:
        merged_chip_patterns = merge_executor.map(lambda pf: pf.merge_chip_patterns(min_depth), all_duplex_modbeds)
 
    fig = plt.figure(figsize=(89/25.4, 120/25.4), dpi=600, layout="constrained")

    sns.set_style("ticks")
    mpl.rc('font', size=5)

    gs = GridSpec(4, 2, fig)

    ax1 = fig.add_subplot(gs[0, :])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1])
    ax4 = fig.add_subplot(gs[2:, :])

    def eval_basecall_proportions(df: pd.DataFrame):
        count_sum = df["Count"].sum()
        df.eval("Proportion = Count / @count_sum", inplace=True)
        
        return df

    def merge_pie_replicates(df: pd.DataFrame):
        df = (df.groupby("Pattern")
              .sum(numeric_only=True)
              .reset_index())
        df = eval_basecall_proportions(df)

        return df

    with concurrent.futures.ThreadPoolExecutor(4) as pie_executor:
        motif_basecalls = pie_executor.map(lambda pf: pf.piechart_data(), merged_motif_patterns)
        motif_basecalls = pd.concat([eval_basecall_proportions(pie.assign(Replicate = i)) for i, pie in enumerate(motif_basecalls)])

    with concurrent.futures.ThreadPoolExecutor(4) as pie_executor:
        chip_patterns = pie_executor.map(lambda pf: pf.piechart_data(), merged_chip_patterns)
        chip_basecalls = pd.concat([eval_basecall_proportions(pie.assign(Replicate = i)) for i, pie in enumerate(chip_patterns)])
        
        chip_piechart = merge_pie_replicates(chip_basecalls)

    chip_piechart["Pattern"] = pd.Categorical(chip_piechart["Pattern"], 
                                              categories=["C:C", "5mC:5mC", "5hmC:5hmC", "C:5mC", "C:5hmC", "5mC:5hmC"],
                                              ordered=True)

    chip_piechart = chip_piechart.sort_values("Pattern")



    # ax2.pie(chip_piechart["Proportion"], 
    #         labels=chip_piechart["Pattern"], 
    #         colors=sns.color_palette("viridis", 6),
    #         startangle=90,
    #         counterclock=False,
    #         center=(0, 0),
    #         explode=(0.15, 0, 0, 0, 0, 0),
    #         autopct='%.1f',
    #         radius=.85,
    #         labeldistance=1.45, pctdistance=1.13,
    #         wedgeprops=dict(linewidth = 0.1))
    
    # ax2.set_title("Duplex patterns\nat CTCF peaks", loc="center", y=.9)
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
       
    genome_basecalls = pd.concat(map(eval_basecall_proportions, genome_basecall_proportions))
    
    motif_basecalls["Type"], chip_basecalls["Type"], genome_basecalls["Type"] = "CTCF motif", "ChIP summit", "Genome average"
    motif_vs_chip = pd.concat([chip_basecalls, motif_basecalls, genome_basecalls])

    sns.barplot(motif_vs_chip,
                x="Pattern", y="Proportion", 
                hue="Type",
                palette="PuBuGn",
                err_kws={'linewidth': 0.8},
                capsize=.5,
                errorbar=("sd", 1),
                order=["C:C", "5mC:5mC", "5hmC:5hmC", "C:5mC", "C:5hmC", "5mC:5hmC"],
                hue_order=["Genome average", "CTCF motif", "ChIP summit"],
                ax=ax1)
    
    ax1.set_xlabel(None) # type: ignore
    ax1.set_ylabel("Proportion of duplex base-calls")
    sns.move_legend(ax1, "upper right", frameon=False, title=None, ncols=3)

    test_groups = motif_vs_chip.groupby(["Pattern", "Type"])
    def ttester(pattern: str, types: list[str]):
        t1 = test_groups.get_group((pattern, types[0]))["Proportion"]
        t2 = test_groups.get_group((pattern, types[1]))["Proportion"]

        test = stats.ttest_rel(t1, t2)
        return print("T-Test:", pattern, types, round(test.pvalue, 4))
    
    [ttester(pattern, ["Genome average", "CTCF motif"]) for pattern in motif_vs_chip["Pattern"].unique()]
    [ttester(pattern, ["CTCF motif", "ChIP summit"]) for pattern in motif_vs_chip["Pattern"].unique()]

    # Locate heterogenous CpG sites at ChIP summits

    all_summits = (pd.concat([modbed.chip_patterns.assign(Replicate =  i+1) for i, modbed in enumerate(all_duplex_modbeds)])
                   .groupby(["Chromosome", "Start", "End", "Strand_ChIP"], observed=True)
                   .sum(numeric_only=True)
                   .reset_index())

    all_summits["Hemi_Motif_C"] = all_summits.eval("CM + CH").where((all_summits["Strand_ChIP"] == "+"), 
                                                                     all_summits.eval("MC + HC"))
    all_summits["All_Motif_5mC"] = all_summits.eval("MM + MH + MC").where((all_summits["Strand_ChIP"] == "+"), 
                                                                         all_summits.eval("MM + HM + CM"))
    all_summits["All_Motif_5hmC"] = all_summits.eval("HM + HH + HC").where((all_summits["Strand_ChIP"] == "+"), 
                                                                            all_summits.eval("MH + HH + CH"))
    
    pie_all_data = pd.DataFrame({"mod" : ["C:C", "Hemi-C", "5mC" , "5hmC"],
                                 "count" : [all_summits["CC"].sum(),
                                            all_summits["Hemi_Motif_C"].sum(),
                                            all_summits["All_Motif_5mC"].sum(), 
                                            all_summits["All_Motif_5hmC"].sum()]})
    
    all_summits["Hemi_Motif_5mC"] = all_summits["CM"].where((all_summits["Strand_ChIP"] == "+"), 
                                                                         all_summits["MC"])

    all_summits["Hemi_Motif_5hmC"] = all_summits["CH"].where((all_summits["Strand_ChIP"] == "+"), 
                                                                           all_summits["HC"])
    mini_pie_data = pd.DataFrame({"mod" : ["C:5mC", "C:5hmC"],
                                  "count" : [all_summits["Hemi_Motif_5mC"].sum(), 
                                             all_summits["Hemi_Motif_5hmC"].sum()]})
    
    def autopcter(pct, allvals):
        """
        From matplotlib: "Labeling a pie and a donut". 
        """
        absolute = int(np.round(pct/100.*np.sum(allvals)))
        return f"{pct:.1f}%\n({absolute:.1e})"

    ax2.pie(pie_all_data["count"], 
            labels=pie_all_data["mod"],
            explode=(0, 0.25, 0, 0),
            colors=sns.color_palette("Paired", 4),
            startangle=40, 
            autopct=lambda pct: autopcter(pct, pie_all_data["count"].to_list()),
            radius=1.25)

    ax2.set_title("Base-call at motif", fontsize=5)
    
    ax3.pie(mini_pie_data["count"],
            labels=mini_pie_data["mod"],
            colors=sns.color_palette("Paired", 4)[2:],
            radius=.75,
            autopct=lambda pct: autopcter(pct, mini_pie_data["count"].to_list()),
            startangle=90)
    
    ax3.set_title("Base-call opposite C", fontsize=5)

    for col in ["CC",  "CH",  "CM",  "HC",  "HH",  "HM",  "MC",  "MH",  "MM"]:
        all_summits[col] = all_summits.eval(f"({col} / readCount)*100")

    multiple_sites = (all_summits.groupby(["Chromosome", "Start_Summit", "End_Summit"], as_index=True)
                      .sum(numeric_only=True))
    multiple_sites["CpGs_in_summit"] = all_summits.groupby(["Chromosome", "Start_Summit", "End_Summit"], 
                                                           as_index=True)["Start"].count()

    with pd.ExcelWriter("data/ctcf/ctcf_intersects.xlsx", mode="w") as writer:
        all_summits.to_excel(writer, columns=
                             ["Chromosome", "Start", "End", "CC", "MM", "HH", "Hemi_Motif_C", "All_Motif_5mC", "All_Motif_5hmC", "Hemi_Motif_5mC", "Hemi_Motif_5hmC", "readCount"], 
                             sheet_name="CpG_Sites", index=False)
        multiple_sites.reset_index().to_excel(writer, 
                                              columns=["Chromosome", "Start_Summit", "End_Summit", "CC",  "CH",  "CM",  "HC",  "HH",  "HM",  "MC",  "MH",  "MM", "CpGs_in_summit"],
                                              sheet_name="Motifs", index=False)
        
    ax4.tick_params("both", bottom=False, labelbottom=False, left=False, labelleft=False)

    for i, ax in enumerate([ax1, ax2, ax4]):
        ax.set_title(string.ascii_lowercase[i], fontdict={"fontweight" : "bold"}, loc="left")        

    sns.despine()
    sns.despine(ax=ax4, left=True, bottom=True)
        
    if test_run:
        fig.savefig("plots/tests/ctcf_duplex_analysis.png")

    else:
        fig.savefig("plots/ctcf_duplex_analysis.png")
        fig.savefig("plots/ctcf_duplex_analysis.svg")

    # Stats for ax1 and ax2

    def comparer(exp, obs):
        merged_replicates = motif_vs_chip.groupby(["Type", "Pattern"]).sum(numeric_only=True)
        obs = merged_replicates.groupby("Type").get_group(obs).reset_index()
        exp = merged_replicates.groupby("Type").get_group(exp).reset_index()

        exp_sum = exp["Count"].sum()
        exp["Proportion"] = exp.eval("Count / @exp_sum")

        obs_sum = obs["Count"].sum()
        exp_freq = exp["Proportion"].mul(obs_sum)

        # print(obs_sum, exp_freq.sum())

        assert round(exp_freq.sum()) == obs_sum

        data = np.array([obs["Count"].to_numpy(), exp_freq.to_numpy()])
        data = data.astype(int)

        # print(data)

        v = contingency.association(data)
        stat, p, dof, _ = stats.chi2_contingency(data, lambda_="log-likelihood")        
        return v, p
    
    print("Genome vs. Genome sanity test :", comparer("Genome average", "Genome average"))
    print("Genome vs. Motif :", comparer("Genome average", "CTCF motif"))
    print("Genome vs. Summit :", comparer("Genome average", "ChIP summit"))
    print("Motif vs. Summit :", comparer("CTCF motif", "ChIP summit"))

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
                        default=10,
                        required=False,
                        help="Whether to filter the dataset by depth at CpG.")
    parser.add_argument("--upset_plot", 
                        dest="upset_plot", 
                        action="store_true",
                        default=False,
                        help="Whether to filter the dataset by depth at CpG.")


    args = parser.parse_args()

    if args.test_run:
        test_run = True
    else: 
        test_run = False

    args = parser.parse_args()
    main(test_run, args.min_depth, args.upset_plot)    