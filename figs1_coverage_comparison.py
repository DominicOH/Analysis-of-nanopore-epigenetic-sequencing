"""
Will read the table files to compare genomic coverage depth of multiple datasets. 

Takes a while to run as it needs to generate new intermediate data files for each dataset. 
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
from AnalysisTools.helpers import timer
import pandas as pd
import string
import subprocess
import concurrent.futures
import os
from itertools import chain
from scipy import stats

CYTOSINES_IN_CPG = 18518250 # Count of all genomic CpG sites. Ignores masked positions. Produced using: 
# awk '{IGNORECASE=0}/CG/{ ++count }END{ print count*2 }' mm39.fa 

class GenomeCov:
    """
    File produced by bedtools genomecov with -ibam <sorted_alignments>.
    """
    def __init__(self, path, replicate) -> None:
        self.path = path
        self.replicate = replicate
        self.__df = self.__read_genomecov() 
        self.__median = self.median_genomecov()

    def __read_genomecov(self) -> pd.DataFrame:
        names = ["SeqName", "Depth", "NBases", "ChromSize", "FractionCoverage"]
        df = (pd.read_table(self.path, names=names)
              .query("SeqName == 'genome'")
              .assign(Replicate = self.replicate)
              .reset_index(drop=True))

        return df
    
    @property
    def df(self):
        return self.__df
    
    def median_genomecov(self):
        """
        Calculate median coverage from the genomecov dataset.
        """
        df = self.df

        cov_fract = 0
        iterrows = df.iterrows()
        next(iterrows) # the first line is always depth=0, so skipping this to ensure it's the median depth of covered sites
        for _, row in iterrows:
            cov_fract += row["FractionCoverage"]
            if cov_fract >= .5:
                median_depth = row["Depth"]                    
                return median_depth     
    
    @property
    def median(self):
        return self.__median
    
    def covered_proportion(self):
        df = self.df
        return (1 - df.loc[df["Depth"] < 5, "FractionCoverage"].sum())
   
class FeatureCoverage:
    """
    Files produced by AnalysisTools/median_feature_depth, which calculates the median depth of CpG sites of a given type of feature. 
    """
    def __init__(self, path) -> None:
        self.path = path
        self.__df = self.__read_coverage()
        pass
    
    def __read_coverage(self):        
        return pd.read_table(self.path, names=["Context", "MedianDepth"]).assign(FPath = os.path.basename(self.path))
    
    @property
    def df(self):
        return self.__df
    
class Bed4:
    """
    Bed files where the 4th field is comprised of the depth at a position. Also accepts Nanopore bedMethyl files.
    """
    def __init__(self, path) -> None:
        self.path = path
        self.__kind = self.__detect_kind()
        self.__df = self.__read_file()
        self.__median = self.__calculate_median()

        if print_stats:
            print(self.median)

        pass

    def __detect_kind(self):
        with open(self.path) as file:
            fline = file.readline()
        if len(fline.split()) < 6:
            return "Bismark bed4 converted"
        else: 
            return "Nanopore ModKit"
        
    def __read_file(self):
        # print(f"Reading {self.path}...")
        if self.__kind == "Bismark bed4 converted":
                    usecols = 3
                    skiprows = None
        else: 
            usecols = 4
            skiprows = lambda x: x % 2 == 0

        if test_run:
            nrows = 100000
        else:
            nrows = None
        names = ["Count", "Depth"]
                
        return pd.read_table(self.path, sep="\t", usecols=[0, usecols], nrows=nrows, skiprows=skiprows, names=names) 

    @property
    def df(self):
        return self.__df
    
    def __calculate_median(self):
        return self.df["Depth"].median()
    
    @property
    def median(self):
        return self.__median

    def group_depth(self):
        df = self.df
        
        gb = df.groupby("Depth", as_index=False).count()
        return gb
    
    def compare_to_features(self, path):
        feature_depth = FeatureCoverage(path).df
        features = feature_depth.assign(DiffToGenome = lambda r: ((r["MedianDepth"] -  self.median) / self.median)*100,
                                        GPath = os.path.basename(self.path),
                                        GMedian = self.median)
        return features
    
    def compare_to_gc(self, path):
        gc = GCBed.build(path)

        gc = gc.assign(DiffToGenome = lambda r: ((r["MedianDepth"] -  self.median) / self.median)*100)
        return gc

class GCBed:
    """
    Takes CpG sites from a file, such as a bedMethyl or bismark BED file, and bins sites according to %GC using a %GC reference

    Instructions for producing this are as follows: 
    1. bedtools makewindows is used to produce non-overlapping 100bp windows from the mm39 genome file (UCSC mm39.chrom.sizes)
    2. UCSC bigWigAverageOverBed is run using the UCSC gc5BaseBw file for mm39 and the windowed bed file above with the --bedOut option.
    3. bedtools intersect is run with: -a <bedMethyl | bismark output> and -b <the GC% reference bed in 2.> -wa -wb 
    4. cut is used to select relevant columns; here: BED3 + Depth + GC%
    """
    def __init__(self, path, df) -> None:
        self.path = path 
        self.df = df
        pass

    @classmethod
    def build(cls, path) -> pd.DataFrame:     
        """
        Returns the median depth for CpG sites of a given GC% bin.
        """
        names=["Chromosome", "Start", "End", "MedianDepth", "GC%"]
        if test_run:
            nrows = 100000
        else:
            nrows = None

        gc_df = pd.read_table(path, names=names, nrows=nrows)  
        gc_df["GCBin"] = pd.cut(gc_df["GC%"], bins=np.arange(0, 105, 5), labels=np.arange(5, 105, 5))
        gb = gc_df.groupby("GCBin", as_index=False)["MedianDepth"].median()

        return cls(path, gb).df

def files_from_dir(dirpath, **kwargs) -> list[str]:
    """
    Fetches filepaths from directory. 
    """    
    directory_ls = [dirpath + path for path in subprocess.check_output(["ls", dirpath], **kwargs).decode("utf-8").split("\n")]
    directory_ls.pop(-1)
    return directory_ls

def open_genomecovs(dirpath):
    return [GenomeCov(path, i) for i, path in enumerate(files_from_dir(dirpath))]

def open_bed4(dirpath, **kwargs):
    """
    **kwargs: Keywords passed to subprocess.check_output().
    """
    return [Bed4(path) for path in files_from_dir(dirpath, **kwargs)]

def make_depth_hist(data: list[Bed4], ax: plt.Axes, replicates: list[str], **kwargs):
    with concurrent.futures.ProcessPoolExecutor(len(data)) as tpe:
        grouped_depths = tpe.map(Bed4.group_depth, data)
        depths = pd.concat([df.assign(Replicate = i) for df, i in zip(grouped_depths, replicates)]).reset_index(drop=True)
    
    depths = depths.assign(PercentBases = lambda r: (r["Count"]/depths["Count"].sum())*100)
    plot = sns.lineplot(depths, x="Depth", y="PercentBases", 
                        hue="Replicate", 
                        ax=ax,
                        **kwargs) 
    ax.set_ylabel("Percentage of CpG sites (%)")
    ax.set_xlim(0, 30) 
    sns.move_legend(ax, "best", frameon=False)
        
    return plot

def make_coverage_comparison(data: list[GenomeCov], ax: plt.Axes):
    all_proportions = pd.DataFrame({
        "Tech" : [*["Nanopore"]*4, *["oxBS"]*2, *["TAB"]*3],
        "Replicate" : [genomecov.replicate for genomecov in data], 
        "Proportion covered" : [genomecov.covered_proportion() for genomecov in data]})

    sns.barplot(all_proportions, 
                y="Tech", x="Proportion covered",
                hue="Tech", palette="Accent",
                errorbar="sd", err_kws={"lw": .6}, capsize=.3, 
                order=["Nanopore", "oxBS", "TAB"],
                hue_order=["Nanopore", "TAB", "oxBS"],
                width=.6,
                ax=ax)
    
    sns.stripplot(all_proportions, 
                  y="Tech", x="Proportion covered",
                  color="k", size=3,
                  ax=ax)
    
    ax.set_ylabel(None)
    ax.set_xlabel(r"Proportion covered at depth$\geq$5x")
    return

def tech_context_depth(bed4s: list[Bed4], path, tech):
    feature_files = files_from_dir(path)
    all_reps = pd.concat([bed4.compare_to_features(path).assign(Replicate = i, Tech = tech)
                          for i, bed4, path in zip(range(len(bed4s)), bed4s, feature_files)])
    
    return all_reps

@timer
def main():
    sns.set_style("ticks")
    mpl.rc('font', size=5)

    fig = plt.figure(figsize=(180/25.4, 120/25.4), dpi=600, layout="constrained")
    gs = GridSpec(3, 6, fig)

    nano_hist_ax = fig.add_subplot(gs[0, :2])
    oxbs_hist_ax = fig.add_subplot(gs[0, 2:4])
    tab_hist_ax = fig.add_subplot(gs[0, 4:])
    median_barplot_ax = fig.add_subplot(gs[1, 0])
    median_cpg_ax = fig.add_subplot(gs[1, 1])
    gc_depth_ax = fig.add_subplot(gs[1, 2:4])
    gc_vmedian_ax = fig.add_subplot(gs[1, 4:])
    feature_coverage_ax = fig.add_subplot(gs[2, :4])
    coverage_comparison_ax = fig.add_subplot(gs[2, 4:])

    print("Loading bedtools genomecov datasets...")    
    genomecov_paths = ["/mnt/data1/doh28/data/nanopore_hmc_validation/nanopore_wgs/coverage_stats/genome/", 
                       "/mnt/data1/doh28/data/public/CRR008807_TAB/coverage_stats/genome/", 
                       "/mnt/data1/doh28/data/public/CRR008808_oxBS/coverage_stats/genome/"]

    nanopore_genomecovs, tab_genomecovs, oxbs_genomecovs = map(open_genomecovs, genomecov_paths)

    all_medians = pd.DataFrame({
        "Tech" : [*["Nanopore"]*4, *["TAB"]*3, *["oxBS"]*2],
        "Replicate" : [genome_cov.replicate for genome_cov in [*nanopore_genomecovs, *tab_genomecovs, *oxbs_genomecovs]],
        "Median" : [genome_cov.median for genome_cov in [*nanopore_genomecovs, *tab_genomecovs, *oxbs_genomecovs]]
    })

    print("Plotting median depth")
    sns.barplot(all_medians, 
                x="Tech", y="Median",
                hue="Tech", palette="Accent",
                order=["TAB", "oxBS", "Nanopore"],
                errorbar="sd", err_kws={"lw": .6}, capsize=.3, width=.6,
                ax=median_barplot_ax)
    
    sns.stripplot(all_medians, 
                  x="Tech", y="Median",
                  order=["TAB", "oxBS", "Nanopore"],
                  color="k", size=3, 
                  ax=median_barplot_ax)
    
    median_barplot_ax.set_ylabel("Median depth (genomic)")
    median_barplot_ax.set_xlabel(None)

    if print_stats:
        print("Median depths by dataset:")
        print(all_medians)

    ### Line plots of coverage depth per c ###
    bed_paths = ["data/modbases/nanopore/both_reps/", 
                 "data/modbases/public/CRR008808_oxBS/bed_convert/", 
                 "data/modbases/public/CRR008807_TAB/bed_convert/"]
    all_beds = [open_bed4(path) for path in bed_paths]

    print("Plotting C-depth lineplots")
    nano_hists, ox_hists, tab_hists  = *all_beds[:4], *all_beds[4:7], *all_beds[7:] 

    make_depth_hist(nano_hists, nano_hist_ax, ["CBM2_1", "CBM2_2", "CBM3_1","CBM3_2"], palette="Greens")
    make_depth_hist(ox_hists, oxbs_hist_ax, ["1", "2"], palette="Oranges")
    make_depth_hist(tab_hists, tab_hist_ax, ["1", "2", "3"], palette="Purples")
    
    for ax, title in zip([nano_hist_ax, oxbs_hist_ax, tab_hist_ax], ["Nanopore", "oxBS-seq", "TAB-seq"]):
        ax.set_xlabel("Depth at CpG")
        ax.set_title(title, ha="center")

    print("Plotting median CpG depth")

    cpg_median_depth_df = pd.DataFrame({
        "Tech" : [*["Nanopore"]*4, *["TAB"]*3, *["oxBS"]*2],
        "Replicate" : [i for i in range(9)], 
        "CpGDepth" : [bed.median for bed in chain.from_iterable(all_beds)]  
    })

    sns.barplot(cpg_median_depth_df,
                x="Tech", y="CpGDepth",
                hue="Tech", palette="Accent",
                order=["TAB", "oxBS", "Nanopore"],
                errorbar="sd", err_kws={"lw": .6}, capsize=.3, width=.6,
                ax=median_cpg_ax)
    
    sns.stripplot(cpg_median_depth_df,
                  x="Tech", y="CpGDepth",
                  hue="Tech", color="k", size=3,
                  order=["TAB", "oxBS", "Nanopore"],
                  ax=median_cpg_ax)
    
    median_cpg_ax.set_xlabel(None)
    median_cpg_ax.set_ylabel("Median depth at CpG (stranded)")

    for ax in [median_barplot_ax, median_cpg_ax]:
        ax.tick_params("x", labelrotation=15)

    print("Plotting coverage comparison")
    make_coverage_comparison([*nanopore_genomecovs, *oxbs_genomecovs, *tab_genomecovs], coverage_comparison_ax)

    sns.despine()   
    [ax.set_title(letter, loc="left", fontweight="bold") for ax, letter in zip(fig.axes, string.ascii_lowercase)]

    ### Feature context coverage ###
    print("Plotting context/depth comparison")
    feature_comparison_paths = ["data/depth_analysis/nanopore/",
                                "data/depth_analysis/oxbs/",
                                "data/depth_analysis/tab/"]
    
    all_feature_comparisons = pd.concat([tech_context_depth(bed4s, path, tech)
                                         for bed4s, path, tech in zip(all_beds,
                                                                      feature_comparison_paths, 
                                                                      ["Nanopore", "oxBS", "TAB"])])
    sns.stripplot(all_feature_comparisons, 
                  x="Context", y="DiffToGenome",
                  hue="Tech", palette="Accent", 
                  order=["Intergenic", "Gene", "Exons", "Intron", "Promoter", "CGI"], 
                  hue_order=["Nanopore", "TAB", "oxBS"],
                  dodge=True, size=5, jitter=True,
                  ax=feature_coverage_ax)
        
    feature_coverage_ax.set_ylabel("Difference to median\nCpG depth (%)")
    feature_coverage_ax.set_xlabel("CpG context")
    feature_coverage_ax.axhline(0, lw=.8, c="grey", ls=":")
    feature_coverage_ax.axvline(1.5, lw=.8, c="k")
    feature_coverage_ax.axvline(3.5, lw=.8, c="k")
    
    sns.move_legend(feature_coverage_ax, "lower left", title=None)

    ### GC% depth comparison ### 
    # mouse mean gc % content: 
    
    mm39 = pd.read_table("/mnt/data1/doh28/data/reference_genomes/mm39/mm39.fa.nucBed")
    total_bases = mm39[["6_num_A", "7_num_C", "8_num_G", "9_num_T"]].sum().sum()
    total_gc = mm39[["7_num_C", "8_num_G"]].sum().sum()

    mean_gc = (total_gc/total_bases)*100


    print("Plotting GC% context vs depth")
    # Follow documentation under GCBed for instructions on producing these files.
    
    gc_paths = files_from_dir("data/depth_analysis/gc/")
    with concurrent.futures.ProcessPoolExecutor(3) as ppe:
        gc_futures = [ppe.submit(bed.compare_to_gc, path)
                      for bed, path in zip(chain.from_iterable(all_beds), gc_paths)]
        gc_data = pd.concat([future.result().assign(Tech = tech,
                                                    Replicate = i) 
                             for future, tech, i in zip(gc_futures, 
                                                        [*["Nanopore"]*4, *["TAB"]*3, *["oxBS"]*2],
                                                        range(len(gc_futures)))]).reset_index()

    sns.lineplot(gc_data,
                 x="GCBin", y="MedianDepth", 
                 hue="Tech", palette="Accent", 
                 lw=1, style="Tech",
                 ax=gc_depth_ax)
    
    gc_depth_ax.axvline(mean_gc, c="grey", ls="dashdot", lw=0.8)
    
    handles_median, leg_labels_median = gc_depth_ax.get_legend_handles_labels()
    gc_data_grouped = gc_data.groupby("Tech")

    new_labels_median = []
    new_labels_percentage = []
    for group, label in zip(gc_data_grouped.groups, leg_labels_median):
        data = gc_data.groupby("Tech").get_group(group).dropna()
        median_stat = stats.spearmanr(data["MedianDepth"], data["GCBin"]).statistic
        median_p = stats.spearmanr(data["MedianDepth"], data["GCBin"]).pvalue

        print(f"Median depth vs. GC% for {group}: $\\rho$={median_stat}, p={median_p}, n={len(data)}")

        percentage_stat = stats.spearmanr(data["DiffToGenome"], data["GCBin"]).statistic
        percentage_p = stats.spearmanr(data["DiffToGenome"], data["GCBin"]).pvalue

        print(f"Percentage from median depth vs. GC% for {group}: $\\rho$={percentage_stat}, p={percentage_p}, n={len(data)}")

        new_labels_median.append(label + f", $\\rho$={round(median_stat, 2)}")
        new_labels_percentage.append(label + f", $\\rho$={round(percentage_stat, 2)}")

    gc_depth_ax.legend(handles_median, new_labels_median)
    gc_depth_ax.set_ylabel("Depth at CpG")
    
    sns.lineplot(gc_data,
                 x="GCBin", y="DiffToGenome", 
                 hue="Tech", palette="Accent", 
                 lw=1, style="Tech",
                 ax=gc_vmedian_ax)
    
    handles_percent, _ = gc_vmedian_ax.get_legend_handles_labels()
    gc_vmedian_ax.axvline(mean_gc, c="grey", ls="dashdot", lw=0.8)
    gc_vmedian_ax.set_ylabel("Difference to median\nCpG depth (%)")
    gc_vmedian_ax.legend(handles_percent, new_labels_percentage)

    for ax in [gc_depth_ax, gc_vmedian_ax]:
        ax.set_xlim(5, 100)
        ax.set_xlabel("GC in 100bp")
        sns.move_legend(ax, "lower left", 
                    title=None, frameon=True,
                    bbox_to_anchor=(.5, .8)
                    )

    fig.savefig("plots/compare_coverage.png")
    if not test_run:
        fig.savefig("plots/compare_coverage.svg")
    return 

### main function ##### 
if __name__=="__main__":
    parser = argparse.ArgumentParser(
        prog = "compare_coverage",
        description = "Compares the coverage of the different datasets.")
    parser.add_argument("-t", "--test-run", dest="test_run",  action="store_true", default=False)
    parser.add_argument("--print-stats", dest="print_stats",  action="store_true", default=False)

args = parser.parse_args()

global test_run
test_run = args.test_run

global print_stats
print_stats = args.print_stats

main()
print("Completed.")