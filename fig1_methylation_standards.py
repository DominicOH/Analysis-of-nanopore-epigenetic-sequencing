import seaborn as sns
from sklearn import metrics
from sklearn.utils import resample
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import patches, lines
import logomaker
import argparse
import numpy as np
import gc
from AnalysisTools.common import fetch_controls
import concurrent.futures
from string import ascii_lowercase
from AnalysisTools.helpers import timer
import pandas as pd
import pyranges as pr
from scipy import stats

class Control:
    def __init__(self, df: pd.DataFrame) -> None:
        self.df = df

    def calculate_fpr_nonspecific(self):
        fpr = ((self.df["N_5mC"].sum() + self.df["N_5hmC"].sum())/self.df["readCount"].sum())*100
        return fpr

    def calculate_fpr_modspecific(self, mod):
        fpr = (self.df[mod].sum()/self.df["readCount"].sum())*100

        return fpr
    
    def calculate_feature_fpr(self, ref_pr: pr.PyRanges, mod=None) -> pd.DataFrame:
        annotated_error = pr.PyRanges(self.df).join(ref_pr, suffix="_Feature").as_df()
        target_repeat_types = ["Alu", "MIR", "L1", "L2", "Low_complexity", "Simple_repeat", "Satellite", "CGI"]

        if test_run:
            target_repeat_types.remove("Satellite") # can't find any satellites in the first 1 million basecalls

        if not mod:
            grouped_by_feature = (annotated_error
                    .groupby("Score", as_index=True)[["readCount", "N_5hmC", "N_5mC"]]
                    .sum()
                    .loc[target_repeat_types]
                    .assign(FPR_5mC = lambda r: r["N_5mC"]/r["readCount"], 
                            FPR_5hmC = lambda r: r["N_5hmC"]/r["readCount"])
                            .reset_index()
                            .melt(id_vars=["Score", "readCount"], value_vars=["FPR_5mC", "FPR_5hmC"]))
        else:
            grouped_by_feature = (annotated_error
                    .groupby("Score", as_index=True)[["readCount", mod]]
                    .sum()
                    .loc[target_repeat_types]
                    .assign(FPR = lambda r: r[mod]/r["readCount"])
                            .reset_index())
            
        grouped_by_feature = (grouped_by_feature.replace(["DNA", "Simple_repeat", "Low_complexity"], ["Repeat", "Simple repeat", "Low complexity"]))
        return grouped_by_feature

def concat_controls(controls: list[Control]):
    concat = (pd.concat([control.df for control in controls])
              .groupby(["Chromosome", "Start", "End"], observed=True, as_index=False, sort=False)
              .sum(numeric_only=True))
    return concat.drop(columns=["Chromosome", "Start", "End"])

def precision_recall(modified_control: pd.DataFrame, 
                     unmodified_control: pd.DataFrame, 
                     ax):
    
    print("TPR:", (modified_control["N_5mC"].sum()/modified_control["readCount"].sum())*100) # TPR
    print("FPR:", ((unmodified_control["N_5mC"].sum() + unmodified_control["N_5hmC"].sum())/unmodified_control["readCount"].sum())*100) # FPR

    # Note: 5hmC basecalls are considered false negatives in the methylated dataset (as they are off-target) but are false positives in the unmodified control. 

    pred = np.concatenate([np.ones(int(modified_control["N_5mC"].sum())), # True positive calls in the methylated control
                           np.zeros(int(modified_control["N_C"].sum() + modified_control["N_5hmC"].sum())), # False negative calls in the methylated control
                           np.ones(int(unmodified_control["N_5mC"].sum() + unmodified_control["N_5hmC"].sum())), # False positive (5mC or 5hmC) calls in the unmodified control
                           np.zeros(int(unmodified_control["N_C"].sum()))])# True negative calls (C) calls in the unmodified data
    
    gt = np.concatenate([np.ones(int(modified_control["readCount"].sum())), # All CpG basecalls in the methylated control (assumed ground truth positive)
                         np.zeros(int(unmodified_control["readCount"].sum()))]) # All CpG basecalls in the unmodified control (assumed ground truth negative)
    
    if not test_run:
        n_samples=100000000
    else:
        n_samples=1000000
    print(f"Resampling to {n_samples} basecalls and plot PR Curve")

    gt_r, pred_r = resample(gt, pred, n_samples=n_samples, random_state=42)
    
    print("Plot statistics:")
    print("AUC:", round(metrics.roc_auc_score(gt_r, pred_r), 3))
    print("F1:", metrics.f1_score(gt_r, pred_r))
    print("Precision:", metrics.precision_score(gt_r, pred_r))
    print("Recall:", metrics.recall_score(gt_r, pred_r))
    
    return metrics.PrecisionRecallDisplay.from_predictions(gt_r, pred_r, ax=ax)

def make_confusion_matrix(modified_control: pd.DataFrame, 
                         unmodified_control: pd.DataFrame, 
                         ax):
    
    truth = np.concatenate([np.full(int(unmodified_control["readCount"].sum()), "C"),  
                        np.full(int(modified_control["readCount"].sum()), "5mC")])

    pred = np.concatenate([np.full(int(unmodified_control["N_C"].sum()), "C"),
                           np.full(int(unmodified_control["N_5mC"].sum()), "5mC"),
                           np.full(int(unmodified_control["N_5hmC"].sum()), "5hmC"),
                           np.full(int(modified_control["N_C"].sum()), "C"),
                           np.full(int(modified_control["N_5mC"].sum()), "5mC"),
                           np.full(int(modified_control["N_5hmC"].sum()), "5hmC"),
                       ])
    
    if not test_run:
        n_samples=100000000
    else:
        n_samples=1000000

    print(f"Resampling to {n_samples} basecalls and plotting p Curve")
    truth_resample, pred_resample = resample(truth, pred, n_samples=n_samples, random_state=42)

    assert len(truth) == len(pred)
    cm = metrics.confusion_matrix(truth_resample, pred_resample, normalize="true")
    mat = np.delete(cm, (0), axis=0)

    for i, a in enumerate(mat):
        mat[i]=(a[::-1])

    return sns.heatmap(mat,
                xticklabels=["C", "5mC", "5hmC"], 
                yticklabels=["5mC", "C"], 
                cmap="Blues",
                annot=True, linewidths=1, 
                vmax=.1, vmin=0, cbar=False,
                square=False,
                ax=ax)

def repeat_barplot(controls: list[Control], ax=plt.Axes, mod=None):
    """
    ::param mod:: must be N_5mC or N_5hmC
    """
    feature_ref = pr.read_bed("feature_references/hg38/repeats/hg38_rmsk_cgi.tsv")

    if not mod:   
        with concurrent.futures.ProcessPoolExecutor(2) as ppe:
            futures = [ppe.submit(control.calculate_feature_fpr, feature_ref) for control in controls]
            results = pd.concat([future.result().assign(Replicate = i) for i, future in enumerate(futures)])

        results = results.replace(["FPR_5mC", "FPR_5hmC"], ["5mC", "5hmC"])

        sns.barplot(results, 
                    x="Score", y="value", 
                    hue="variable",
                    err_kws={"lw" : .6}, errorbar="sd", 
                    width=.6, capsize=.3,
                    palette="BuGn",
                    ax=ax)

    else: 
        with concurrent.futures.ProcessPoolExecutor(2) as ppe:
            futures = [ppe.submit(control.calculate_feature_fpr, feature_ref, mod) for control in controls]
            results = pd.concat([future.result().assign(Replicate = i) for i, future in enumerate(futures)])

        sns.barplot(results, 
                    x="Score", y="FPR", 
                    err_kws={"lw" : .6}, errorbar="sd", 
                    width=.3, capsize=.15, legend=True,
                    color=sns.color_palette("BuGn", 2)[1],
                    ax=ax)
        
    ax.set_xlabel(None)
    ax.set_ylabel("FPR")
    return 

def kmer_count(df: pd.DataFrame):
    return df.groupby("query_kmer", as_index=False).agg(count=pd.NamedAgg(column="query_kmer", aggfunc="count"))

def gc_plots(ax: plt.Axes):
    binned_gc_paths = ["zymo_wga_unmodified_rep1.sorted.bedMethyl.percentGC.bed.gcBinned", 
                       "zymo_wga_unmodified_rep2.sorted.bedMethyl.percentGC.bed.gcBinned"]
    base_path = "data/gc/controls/"
    paths = [base_path + path for path in binned_gc_paths]
    names = ["GCBin", "readCount", "N_5mC", "N_5hmC", "FPR"]

    dataset = pd.concat([pd.read_table(path, names=names).assign(Replicate = i) for i, path in enumerate(paths)])
    
    dataset_melt = (dataset.melt(id_vars=["GCBin", "readCount"], value_vars=["N_5mC", "N_5hmC"], value_name="FP", var_name="modType")
                           .assign(FPR = lambda r: r["FP"]/r["readCount"])
                           .replace(["N_5mC", "N_5hmC"], ["5mC", "5hmC"])
                           )

    dataset = dataset.assign(Proportion_of_calls = lambda r: r["readCount"] / dataset["readCount"].sum())

    ax.axhline(y=(dataset["N_5mC"].sum()+dataset["N_5hmC"].sum())/dataset["readCount"].sum(), 
                    c="grey", ls=":", label="Genome mean FPR (Total)")    
    sns.lineplot(dataset_melt, 
                 x="GCBin", y="FPR", 
                 hue="modType", palette="BuGn", lw=.8,
                 ax=ax)
    
    sns.lineplot(dataset, 
                 x="GCBin", y="FPR", 
                 label="Total (5mC+5hmC)", 
                 c="black", lw=.8,
                 ax=ax)
    print("Pearson r: FPR vs. GCBin", stats.pearsonr(dataset["FPR"], dataset["GCBin"]))
    
    sns.move_legend(ax, "upper left", frameon=False, title="False positive call")

    ax.set_xlim(5, 100)
    ax.set_ylim(0)

    ax.set_ylabel("False positive rate")
    ax.set_xlabel("GC in 100bp")    

    return 
    
@timer
def main(test_run=True):
    print("Fetching data")
    controls = fetch_controls(["readCount", "N_C", "N_mC", "N_hmC"], test_run=test_run)
    controls = [control for control in map(lambda df: Control(df), controls)]
    print("Data fetched")

    fig = plt.figure(figsize=(180/25.4, 120/25.4), dpi=600, layout="constrained")
    mpl.rc('font', size=5)
    gs = fig.add_gridspec(5, 6)

    precision_recall_ax = fig.add_subplot(gs[:2, :2])   
    gc_bar_plot_ax = fig.add_subplot(gs[:2, 2:])
    logos_m_ax = fig.add_subplot(gs[2, :3])
    logos_h_ax = fig.add_subplot(gs[2, 3:])
    repeat_unmodified_ax = fig.add_subplot(gs[3, :4])
    repeat_modified_ax = fig.add_subplot(gs[4, :4])
    confusion_matrix_ax = fig.add_subplot(gs[3:, 4:])

    # FPR barplot by context # 

    print("Assessing context specific FPR")
    pos_ctrls, neg_ctrls = controls[:2], controls[2:]   
    del controls
    
    # logomaker presentation of sequence motifs # 
    extract_paths = ["/mnt/data1/doh28/data/nanopore_hmc_validation/nanopore_wgs/zymo_unmodified/extract/" + name 
                     for name in ["zymo_wga_unmodified_rep1.sorted.extract.tsv", "zymo_wga_unmodified_rep2.sorted.extract.tsv"]]
    names = ["readID", "forward_read_position", "ref_position", "chrom", "mod_strand", "ref_strand", "ref_mod_strand", "fw_soft_clipped_start", "fw_soft_clipped_end", "read_length",
             "call_prob", "call_code", "base_qual", "ref_kmer", "query_kmer", "canonical_base", "modified_primary_base", "fail", "inferred", "within_alignment", "flag"]
    
    if test_run:
        nrows=1000
    else:
        nrows=None

    print("Merging control datasets for PR curve")
    m_merge, u_merge = [concat_controls(controls) for controls in [pos_ctrls, neg_ctrls]]

    unmod_fpr = (u_merge["N_5hmC"].sum() + u_merge["N_5mC"].sum())/u_merge["readCount"].sum()
    unmod_hmc_fpr = u_merge["N_5hmC"].sum()/u_merge["readCount"].sum()
    mod_hmc_fpr = m_merge["N_5hmC"].sum()/m_merge["readCount"].sum()
        
    repeat_unmodified_ax.axhline(unmod_fpr, c="grey", ls=":", lw=.8, label="Mean FPR (Total)")

    print("Writing repeat-context specific FPR barplots")
    repeat_barplot(neg_ctrls, repeat_unmodified_ax)
     
    sns.move_legend(repeat_unmodified_ax, "lower left", frameon=False, 
                    ncols=3, title=None, bbox_to_anchor=(.25, 1))
    
    repeat_barplot(pos_ctrls, repeat_modified_ax, "N_5hmC")

    patch_legend = [patches.Patch(color=sns.color_palette("BuGn", 2)[1], label="5hmC"),
                    lines.Line2D([0], [0], c="grey", ls=":", lw=.8, label="Mean FPR (5hmC)")]

    repeat_modified_ax.legend(handles=patch_legend, loc="center")
    repeat_modified_ax.axhline(mod_hmc_fpr, c="grey", ls=":", lw=.8)

    sns.move_legend(repeat_modified_ax, "lower left", frameon=False, 
                    ncols=2, title=None, reverse=True, bbox_to_anchor=(.25, 1))

    repeat_unmodified_ax.set_title("Modification-negative control", y=1.25)
    repeat_modified_ax.set_title("5mC-positive control", y=1.25)

    [ax.axvline(6.5, c="k", lw=.8) for ax in [repeat_unmodified_ax, repeat_modified_ax]]

    print("Writing GC plot")
    gc_plots(gc_bar_plot_ax)

    print("Making motif logo")
    extract_dfs = pd.concat([pd.read_table(path, names=names, usecols=["readID", "call_code", "query_kmer"], nrows=nrows) 
                             for path in extract_paths])
    
    m_only = extract_dfs.query("call_code == 'm'")
    h_only = extract_dfs.query("call_code == 'h'")
    m_only, h_only = map(kmer_count, [m_only, h_only])

    m_only_mat, h_only_mat = [logomaker.alignment_to_matrix(frame["query_kmer"], counts=frame["count"]) 
                              for frame in [m_only, h_only]]

    for mat, ax in zip([m_only_mat, h_only_mat], [logos_m_ax, logos_h_ax]):
        mat = logomaker.transform_matrix(mat, normalize_values=True)
        logomaker.Logo(mat, fade_probabilities=True, color_scheme="base_pairing", ax=ax)

        print("GC%:", mat["G"].mean() + mat["C"].mean())

    del extract_dfs, m_only_mat, h_only_mat
    gc.collect()

    for ax in [logos_m_ax, logos_h_ax]:
        ax.set_ylabel("Base probability")
        ax.set_xlabel("K-mer position")

    logos_m_ax.set_title("5mC", loc="center")
    logos_h_ax.set_title("5hmC", loc="center")

    # making a precision-recall curve
    
    print(f"Total 5hmC FPR: \n5mC-positive standard {mod_hmc_fpr}\nModification negative standard: {unmod_hmc_fpr}")   
    print("Defining Precision-Recall and confusion matrix")
    precision_recall(m_merge, u_merge, precision_recall_ax)
    precision_recall_ax.get_legend().remove()
    make_confusion_matrix(m_merge, u_merge, confusion_matrix_ax)

    print("Finished plotting. Calculating stats using full datasets")
    precision_recall_ax.set_ylim(0, 1)
    precision_recall_ax.set_xlim(0, 1)
    precision_recall_ax.set_ylabel("Precision")
    precision_recall_ax.set_xlabel("Recall")
    # precision_recall_ax.set_aspect("equal")

    confusion_matrix_ax.set_ylabel("True label")
    confusion_matrix_ax.set_xlabel("Predicted label")
    
    for i, ax in enumerate(fig.get_axes()):
        ax.set_title(ascii_lowercase[i], loc="left", fontweight="bold")

    sns.despine()

    fig.savefig("plots/control_roc.png")
    if not test_run:
        fig.savefig("plots/control_roc.svg")

##### Main function #####

if __name__=="__main__":
    parser = argparse.ArgumentParser(
                        prog = "cpg_methylation_compare",
                        description = "Compares the coverage of the different datasets.")
    parser.add_argument("-t ", "--test-run", 
                        action="store_true", 
                        dest="test_run", 
                        required=False,
                        help="Produce a test output using the first 100k base positions from the pileup data.") 

    args = parser.parse_args()

    global test_run
    if args.test_run:
        test_run = True
    else: 
        test_run = False

    main(test_run)    