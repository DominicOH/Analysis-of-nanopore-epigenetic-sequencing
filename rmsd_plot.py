import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import string

def main():
    ox_rmsds = pd.read_table("nanopore_vs_ox_rmsd.tsv").rename(columns={"Bisulphite_rep" : "oxBS-seq rep.",
                                                                        "Comparison" : "Nanopore rep."})
    tab_rmsds = pd.read_table("nanopore_vs_tab_rmsd.tsv").rename(columns={"Bisulphite_rep" : "TAB-seq rep.",
                                                                        "Comparison" : "Nanopore rep."})

    mpl.rc('font', size=5)

    fig, axes = plt.subplots(2, 1, figsize=(120/25.4, 89/25.4), dpi=600, layout="constrained")
    ax1, ax2 = axes

    sns.scatterplot(ox_rmsds, 
                x="Depth", y="RMSD",
                style="oxBS-seq rep.", hue="Nanopore rep.", 
                palette="BuGn",
                ax=ax1)

    sns.scatterplot(tab_rmsds, 
                x="Depth", y="RMSD",
                style="TAB-seq rep.", hue="Nanopore rep.",
                palette="PuBuGn", 
                ax=ax2)

    ax1.set_title("5mC")
    ax2.set_title("5hmC")   

    for i, ax in enumerate(axes):
        sns.move_legend(ax, "lower left", bbox_to_anchor=(-.33, .15), frameon=False)
        ax.set_xticks(range(5, 26, 5))
        ax.set_ylim(0, 25)
        ax.set_xlabel("Minimum depth at CpG site")
        ax.set_title(string.ascii_lowercase[i], fontweight="bold", loc="left")

    sns.despine()
    fig.savefig("plots/rmsd_plot.png")
    fig.savefig("plots/rmsd_plot.svg")

    return print("Done")

if __name__ == "__main__":
    main()