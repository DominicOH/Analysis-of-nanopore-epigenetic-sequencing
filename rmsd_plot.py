import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import string

def main():

    mpl.rc('font', size=5)

    fig, axes = plt.subplots(4, 1, figsize=(120/25.4, 120/25.4), dpi=600, layout="constrained")
    ax1, ax2, ax3, ax4 = axes

    nanopore_ia_mc = pd.read_table("data/rmsd/nanopore_intrassay_5mC.tsv").assign(Modification = "5mC", Method="Nanopore")
    nanopore_ia_hmc = pd.read_table("data/rmsd/nanopore_intrassay_5hmC.tsv").assign(Modification = "5hmC", Method="Nanopore")

    tab_ia = pd.read_table("data/rmsd/tab_intrassay.tsv").assign(Modification = "5hmC", Method="TAB")
    ox_ia = pd.read_table("data/rmsd/oxbs_intrassay.tsv").assign(Modification = "5mC", Method="oxBS")

    all_iarmsd = (pd.concat([nanopore_ia_mc, nanopore_ia_hmc, ox_ia, tab_ia])
                  .drop_duplicates(['Depth', 'Size', 'RMSD', 'Modification', 'Method']))

    sns.lineplot(all_iarmsd, 
                x="Depth", y="RMSD", 
                errorbar=None,
                style="Modification", 
                hue="Method", palette="Accent", hue_order=["Nanopore", "TAB", "oxBS"],
                ax=ax1)
    
    ax1.set_title("Intra-assay deviation")

    sns.lineplot(all_iarmsd, 
                x="Depth", y="Mean_X", 
                style="Modification",
                hue="Method", palette="Accent", hue_order=["Nanopore", "TAB", "oxBS"],
                ax=ax2)
    
    ax2.set_ylabel("Mean CpG\nmodification (%)")
    ax2.set_ylim(0)

    ox_rmsds = pd.read_table("data/rmsd/nanopore_vs_ox_rmsd.tsv").rename(columns={"Other_rep" : "oxBS-seq rep.",
                                                                                  "Comparison" : "Nanopore rep."})
    tab_rmsds = pd.read_table("data/rmsd/nanopore_vs_tab_rmsd.tsv").rename(columns={"Other_rep" : "TAB-seq rep.",
                                                                                    "Comparison" : "Nanopore rep."})
    sns.lineplot(ox_rmsds, 
                x="Depth", y="RMSD",
                errorbar=None,
                style="oxBS-seq rep.", hue="Nanopore rep.", 
                palette="BuGn", 
                ax=ax3)

    sns.lineplot(tab_rmsds, 
                x="Depth", y="RMSD",
                errorbar=None,
                style="TAB-seq rep.", hue="Nanopore rep.",
                palette="PuBuGn",
                ax=ax4)

    ax3.set_title("5mC")
    ax4.set_title("5hmC")   

    for i, ax in enumerate(axes):
        ax.set_xticks(range(5, 16, 1))
        sns.move_legend(ax, "lower left", bbox_to_anchor=(-.31, -.25), frameon=False) 
        ax.set_xlim(5, 15)
        ax.set_xlabel("Minimum depth at CpG site")
        ax.set_title(string.ascii_lowercase[i], fontweight="bold", loc="left", x=-.3)
        ax.get_legend().set_in_layout(False)
    
    for ax in [ax1, ax3, ax4]:
        ax.set_ylim(0, 25) 
        ax.set_yticks(range(0, 26, 5))
    sns.despine()
    fig.savefig("plots/rmsd_plot.png")
    fig.savefig("plots/rmsd_plot.svg")

    return print("Done")

if __name__ == "__main__":
    main()