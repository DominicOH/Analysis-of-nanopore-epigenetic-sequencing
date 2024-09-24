import seaborn as sns
from sklearn import metrics
from sklearn.utils import resample
import matplotlib.pyplot as plt
import matplotlib as mpl
import logomaker
import argparse
import numpy as np
import gc
from AnalysisTools.common import merge_positions, fetch_controls
import concurrent.futures
from string import ascii_lowercase
from AnalysisTools.helpers import timer
import pandas as pd
import pyranges as pr

class Control:
    def __init__(self, df: pd.DataFrame) -> None:
        self.df = df

    def calculate_fpr_nonspecific(self):
        fpr = ((self.df["N_5mC"].sum() + self.df["N_5hmC"].sum())/self.df["readCount"].sum())*100
        return fpr

    def calculate_fpr_modspecific(self, mod):
        fpr = (self.df[mod].sum()/self.df["readCount"].sum())*100

        return fpr

    def calculate_annotated_fpr(self, ref_path):
        feature_annotations = pr.read_bed(ref_path)
        annotated_ctrl = pr.PyRanges(self.df).join(feature_annotations, strandedness=False, suffix="_Feature", apply_strand_suffix=False).as_df()
        fpr = ((annotated_ctrl["N_5mC"].sum() + annotated_ctrl["N_5hmC"].sum())/annotated_ctrl["readCount"].sum())*100

        return fpr
    
    def calculate_annotated_fpr_multiple(self, ref_list):
        with concurrent.futures.ThreadPoolExecutor(len(ref_list)) as tpe: 
            fpr_futures = [tpe.submit(self.calculate_annotated_fpr, ref) for ref in ref_list]
            fpr_results = [future.result() for future in fpr_futures]

        return fpr_results
    
def concat_controls(controls: list[Control]):
    concat = (pd.concat([control.df for control in controls])
              .groupby(["Chromosome", "Start", "End"], observed=True, as_index=False, sort=False)
              .sum(numeric_only=True))
    return concat.drop(columns=["Chromosome", "Start", "End"])

def pr_dataframe(modified_control: pd.DataFrame, 
                 unmodified_control: pd.DataFrame):
    
    print("TPR:", (modified_control["N_5mC"].sum()/modified_control["readCount"].sum())*100) # TPR
    print("FPR:", ((unmodified_control["N_5mC"].sum() + unmodified_control["N_5hmC"].sum())/unmodified_control["readCount"].sum())*100) # FPR

    # Note: 5hmC basecalls are considered false negatives in the methylated dataset (as they are off-target) but are false positives in the unmodified control. 

    pred = np.concatenate([np.ones(int(modified_control["N_5mC"].sum())), # True positive calls in the methylated control
                           np.zeros(int(modified_control["N_C"].sum() + modified_control["N_5hmC"].sum())), # False negative calls in the methylated control
                           np.ones(int(unmodified_control["N_5mC"].sum() + unmodified_control["N_5hmC"].sum())), # False positive (5mC or 5hmC) calls in the unmodified control
                           np.zeros(int(unmodified_control["N_C"].sum()))])# True negative calls (C) calls in the unmodified data
    
    gt = np.concatenate([np.ones(int(modified_control["readCount"].sum())), # All CpG basecalls in the methylated control (assumed ground truth positive)
                         np.zeros(int(unmodified_control["readCount"].sum()))]) # All CpG basecalls in the unmodified control (assumed ground truth negative)
    
    return pred, gt
    
@timer
def main(test_run=True):
    print("Fetching data")
    controls = fetch_controls(["readCount", "N_C", "N_mC", "N_hmC"], test_run=test_run)
    controls = [control for control in map(lambda df: Control(df), controls)]
    print("Data fetched")

    fig = plt.figure(figsize=(120/25.4, 89/25.4), dpi=600, layout="constrained")
    mpl.rc('font', size=5)
    gs = fig.add_gridspec(4, 3)

    pr_plot = fig.add_subplot(gs[:2, 0])   
    neg_context_bplot = fig.add_subplot(gs[:2, 1:])
    logos_m = fig.add_subplot(gs[2, :2])
    logos_h = fig.add_subplot(gs[3, :2])
    hmc_fpr_bplot = fig.add_subplot(gs[2:, 2])

    # 5hmC FPR in unmodified and methylated controls # 

    all_hmc_fprs = map(lambda control: control.calculate_fpr_modspecific("N_5hmC"), controls)
    hmc_fpr_bplot_df = pd.DataFrame({"Control" : ["Methylated", "Methylated", "Unmodified", "Unmodified"],
                                     "Replicate" : [1, 2]*2,
                                     "FPR%" : all_hmc_fprs})
    
    sns.barplot(hmc_fpr_bplot_df,
                x="Control", y="FPR%",
                hue="Control", palette="BuGn", 
                order=["Unmodified", "Methylated"],
                capsize=.25, err_kws={"linewidth" : .8}, errorbar=("sd", 1),
                width=.6,
                ax=hmc_fpr_bplot)
    
    hmc_fpr_bplot.set_ylabel("5hmC false positive rate (%)")
    
    # FPR barplot by context # 

    print("Assessing context specific FPR")
    feature_dir = "feature_references/hg38/"
    refs = [feature_dir + name for name in ["hg38_genes.bed", "hg38_intergenic.bed", "hg38_cpg.bed", "hg38_rmsk.tsv"]]
    pos_ctrls, neg_ctrls = controls[:2], controls[2:]   
    del controls
    gc.collect()

    fprs = []
    for ctrl in neg_ctrls:
        fprs.extend(ctrl.calculate_annotated_fpr_multiple(refs))
        fprs.append(ctrl.calculate_fpr_nonspecific())

    mod_c_fprs = pd.DataFrame({
        "Replicate" : [*[1]*5, *[2]*5],
        "Context" : ["Genic", "Intergenic", "CGI", "Repeat", "Average"]*2,
        "FPR%" : fprs
        })
    
    sns.barplot(mod_c_fprs, 
                x="Context", y="FPR%", 
                hue="Context", palette="Set2",
                capsize=.25, err_kws={"linewidth" : .8}, 
                order=["Genic", "Intergenic", "Repeat", "CGI", "Average"],
                width=.6,
                ax=neg_context_bplot)
    
    neg_context_bplot.set_ylabel("5modC false positive rate (%)")
    
    # logomaker presentation of sequence motifs # 
    extract_paths = ["/mnt/data1/doh28/data/nanopore_hmc_validation/nanopore_wgs/zymo_unmodified/extract/" + name 
                     for name in ["zymo_wga_unmodified_rep1.sorted.extract.tsv", "zymo_wga_unmodified_rep2.sorted.extract.tsv"]]
    names = ["readID", "forward_read_position", "ref_position", "chrom", "mod_strand", "ref_strand", "ref_mod_strand", "fw_soft_clipped_start", "fw_soft_clipped_end", "read_length",
             "call_prob", "call_code", "base_qual", "ref_kmer", "query_kmer", "canonical_base", "modified_primary_base", "fail", "inferred", "within_alignment", "flag"]
    
    if test_run:
        nrows=1000
    else:
        nrows=None

    print("Making motif logo")
    extract_dfs = pd.concat([pd.read_table(path, names=names, usecols=["readID", "call_code", "query_kmer"], nrows=nrows) 
                             for path in extract_paths])
    # extract_dfs = [df.loc]
    m_only = extract_dfs.query("call_code == 'm'")
    h_only = extract_dfs.query("call_code == 'h'")

    m_only_mat, h_only_mat = [logomaker.alignment_to_matrix(frame["query_kmer"]) for frame in [m_only, h_only]]

    for mat, ax in zip([m_only_mat, h_only_mat], [logos_m, logos_h]):
        mat = logomaker.transform_matrix(mat, normalize_values=True)
        logomaker.Logo(mat, fade_probabilities=True, ax=ax)

        print("GC%:", mat["G"].mean() + mat["C"].mean())

    fig.supylabel("Base probability", y=.3, x=.01, fontsize=5)

    del extract_dfs, m_only_mat, h_only_mat
    gc.collect()

    logos_m.set_xticks([])
    logos_h.set_xlabel("K-mer position")

    # making a precision-recall curve
    print("Merging control datasets for PR curve")
    m_merge, u_merge = [concat_controls(controls) for controls in [pos_ctrls, neg_ctrls]]
    del pos_ctrls, neg_ctrls
    gc.collect()

    fpr_hmc = pd.DataFrame({"Control" : ["Unmethylated", "Methylated"], 
                            "FPR%" : [(u_merge["N_5hmC"].sum()/u_merge["readCount"].sum())*100, # how many 5hmC basecalls in the unmethylated control? 
                                      (m_merge["N_5hmC"].sum()/m_merge["readCount"].sum())*100]}) # how many 5hmC basecalls in the methylated control? 
    
    print(f"Total 5hmC FPR% {fpr_hmc}")   
    print("Defining Precision-Recall curve predictions and ground-truth")
    gt, pred = pr_dataframe(m_merge, u_merge)
    del m_merge, u_merge
    gc.collect()

    print("Resampling to 1,000,000 basecalls and plotting PR Curve")
    gt_r, pred_r = resample(gt, pred, n_samples=1000000, random_state=42)

    metrics.PrecisionRecallDisplay.from_predictions(gt_r, pred_r, ax=pr_plot)
    sns.move_legend(pr_plot, "lower left", frameon=False)

    print("Finished plotting. Calculating stats using full datasets")
    pr_plot.set_ylim(0, 1)
    pr_plot.set_xlim(0, 1)
    
    for i, ax in enumerate(fig.get_axes()):
        ax.set_title(ascii_lowercase[i], loc="left", fontweight="bold")

    print("Plot statistics:")
    print("AUC:", round(metrics.roc_auc_score(gt_r, pred_r), 3))
    print("F1:", metrics.f1_score(gt_r, pred_r))

    sns.despine()

    fig.savefig("plots/control_roc.png")
    # fig.savefig("plots/control_roc.svg", dpi=600)

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

    if args.test_run:
        test_run = True
    else: 
        test_run = False

    main(test_run)    