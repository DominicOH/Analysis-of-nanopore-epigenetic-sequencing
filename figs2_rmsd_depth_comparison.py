import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import string


def mean_of_means(df: pd.DataFrame):
    df["Mean"] = df[["Mean_X", "Mean_Y"]].mean(axis=1)
    return df

def main():

    mpl.rc('font', size=5)
    sns.set_style("whitegrid")

    fig, axes = plt.subplots(3, 1, figsize=(89/25.4, 120/25.4), dpi=600, layout="constrained")
    ax1, ax2, ax3 = axes

    nanopore_ia_mc = pd.read_table("data/rmsd/nanopore_intrassay_5mC.tsv").assign(Method="Nanopore 5mC")
    nanopore_ia_hmc = pd.read_table("data/rmsd/nanopore_intrassay_5hmC.tsv").assign(Method="Nanopore 5hmC")

    tab_ia = pd.read_table("data/rmsd/tab_intrassay.tsv").assign(Method="TAB")
    ox_ia = pd.read_table("data/rmsd/oxbs_intrassay.tsv").assign(Method="oxBS")
    
    nanopore_ia_mc, nanopore_ia_hmc, tab_ia, ox_ia = map(mean_of_means, [nanopore_ia_mc, nanopore_ia_hmc, tab_ia, ox_ia])

    all_iarmsd = (pd.concat([nanopore_ia_mc, nanopore_ia_hmc, ox_ia, tab_ia])
                  .drop_duplicates(['Depth', 'Size', 'RMSD']))

    sns.lineplot(all_iarmsd, 
                x="Depth", y="RMSD", 
                errorbar=None,
                hue="Method", style="Method",
                palette="Paired", hue_order=["Nanopore 5mC", "Nanopore 5hmC", "TAB", "oxBS"],
                ax=ax1)
    
    ax1.set_title("Intra-assay deviation")
    
    sns.lineplot(all_iarmsd, 
                x="Depth", y="Mean", 
                errorbar=None,
                hue="Method", style="Method",
                palette="Paired", hue_order=["Nanopore 5mC", "Nanopore 5hmC", "TAB", "oxBS"],
                ax=ax2)
    
    ax2.set_ylabel("Mean CpG\nmodification (%)")
    ax2.set_ylim(0)
    ax2.set_title("Mean modification rate")

    ox_rmsds = pd.read_table("data/rmsd/nanopore_vs_ox_rmsd.tsv").assign(Modification = "5mC (vs. oxBS-seq)")
    tab_rmsds = pd.read_table("data/rmsd/nanopore_vs_tab_rmsd.tsv").assign(Modification = "5hmC (vs. TAB-seq)")

    rmsds = pd.concat([ox_rmsds, tab_rmsds])
    sns.lineplot(rmsds, 
                x="Depth", y="RMSD",
                errorbar=None,
                hue="Modification", 
                palette="GnBu", 
                ax=ax3)  
    
    with pd.ExcelWriter('source_data/figs2_rmsd_source_data.xlsx') as writer:
        all_iarmsd.to_excel(writer, 'figs2a-b_intraassay')
        rmsds.to_excel(writer, 'figs2c_interassay')
    
    ax3.set_title("Inter-assay deviation")

    for i, ax in enumerate(axes):
        ax.set_xticks(range(5, 16, 1))
        sns.move_legend(ax, "lower left", frameon=True, title=None) 
        ax.set_xlim(5, 15)
        ax.set_xlabel("Minimum depth at CpG site")
        ax.set_title(string.ascii_lowercase[i], fontweight="bold", loc="left")
        ax.get_legend().set_in_layout(False)

    sns.move_legend(ax2, "center left")
    
    for ax in [ax1, ax3]:
        ax.set_ylim(0, 25) 
        ax.set_yticks(range(0, 26, 5))
        ax.set_ylabel("Root Mean Square Deviation (%)")

    sns.despine()
    fig.savefig("plots/rmsd_plot.png")
    fig.savefig("plots/rmsd_plot.svg")

    return print("Done")

if __name__ == "__main__":
    main()