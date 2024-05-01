import seaborn as sns
from sklearn import metrics
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import numpy as np
import gc
from AnalysisTools.common import merge_positions, fetch_controls
import concurrent.futures
from string import ascii_lowercase
from AnalysisTools.helpers import timer

@timer
def main(dryrun=True):
    print("Fetching data")
    pos_ctrl, neg_ctrl = fetch_controls(["readCount", "N_C", "N_mC", "N_hmC"], dryrun=dryrun)
    
    print("Data fetched")

    with concurrent.futures.ProcessPoolExecutor(2) as ppe: 
        m_average, u_average = ppe.map(merge_positions, [pos_ctrl, neg_ctrl])

    del pos_ctrl, neg_ctrl
    gc.collect()

    print("TPR:", (m_average["N_5mC"].sum()/m_average["readCount"].sum())*100) # TPR
    print("FPR:", (u_average["N_5mC"].sum()/u_average["readCount"].sum())*100) # FPR

    t_positives = np.ones(int(m_average["N_5mC"].sum())) # All methylated basecalls in the methylated control = TP
    f_negatives = np.zeros(int(m_average["N_C"].sum() + m_average["N_5hmC"].sum())) # All non 5mC basecalls in the methylated control = FN

    f_positives = np.ones(int(u_average["N_5mC"].sum() + u_average["N_5hmC"].sum())) # All methlyated basecalls in the unmethylated control = FP 
    t_negatives = np.zeros(int(u_average["N_C"].sum())) # All unmethlyated basecalls in the unmethylated control = TN 

    gt_positives = np.ones(int(m_average["readCount"].sum())) # All methylated bases in the methylated control = T
    gt_negatives = np.zeros(int(u_average["readCount"].sum())) # All unmethylated bases in the unmethylated control = F

    del m_average, u_average
    gc.collect()

    pred = np.concatenate([t_positives, f_negatives, f_positives, t_negatives])
    gt = np.concatenate([gt_positives, gt_negatives]) 

    fig, axes = plt.subplots(1, 2, figsize=(120/25.4, 60/25.4), layout="constrained")
    ax1, ax2 = axes
    mpl.rc('font', size=5)
    
    print("Plotting ROCs")

    metrics.RocCurveDisplay.from_predictions(gt, pred, ax=ax1)
    metrics.PrecisionRecallDisplay.from_predictions(gt, pred, ax=ax2)

    ax1.set_ylabel("True Positive Rate")
    ax1.set_xlabel("False Positive Rate")
    sns.move_legend(ax1, "lower right", frameon=False)
    
    ax2.set_ylabel("Precision")
    ax2.set_xlabel("Recall")
    sns.move_legend(ax2, "lower left", frameon=False)

    for i, ax in enumerate(axes):
        ax.set_title(ascii_lowercase[i], loc="left", fontweight="bold")
        ax.set_ylim(0, 1)
        ax.set_xlim(0, 1)

    print("AUC:", round(metrics.roc_auc_score(gt, pred), 3))
    print("F1:", metrics.f1_score(gt, pred))

    sns.despine()

    fig.savefig("plots/control_roc.png", dpi=600)
    fig.savefig("plots/control_roc.svg", dpi=600)

##### Main function #####

if __name__=="__main__":
    parser = argparse.ArgumentParser(
                        prog = "cpg_methylation_compare",
                        description = "Compares the coverage of the different datasets.")
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

    main(dryrun)    