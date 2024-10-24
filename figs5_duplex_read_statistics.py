import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from string import ascii_lowercase
import pandas as pd

def replacer(df):
    return df.replace(["-,m,C", "m,-,C", "-,h,C", "h,-,C", "m,h,C", "h,m,C", "-,-,C", "m,m,C", "h,h,C"], 
                      ["C:5mC", "C:5mC", "C:5hmC", "C:5hmC", "5mC:5hmC", "5mC:5hmC", "C", "5mC", "5hmC"])

def make_pie(patterns: pd.DataFrame, mod: str, ax: plt.Axes, palette):
    other_mods = ["-", "m", "h"]
    other_mods.remove(mod)

    mod_df = (patterns.loc[patterns["Pattern"].str.contains(mod)]
              .groupby("Pattern")
              .sum()
              .reset_index())
    
    s1, s2 = [mod_df.loc[mod_df["Pattern"].str.contains(mod_type), "N_Pattern"].sum() for mod_type in other_mods]
    constitutive = (mod_df["N_Pattern"].sum() - s1) - s2

    summary = pd.DataFrame({"Pattern" : [*other_mods, mod], 
                            "Count" : [s1, s2, constitutive]})

    summary = summary.replace(["-", "m", "h"], ["C", "5mC", "5hmC"])
    summary["Pattern"] = pd.Categorical(summary["Pattern"], ["C", "5mC", "5hmC"], ordered=True)
    summary = summary.sort_values("Pattern")

    return ax.pie(summary["Count"], labels=summary["Pattern"], 
                   colors=palette,
                   autopct='%.1f')

def main():
    dx_count_data = pd.read_table("data/duplex_data/dx_counts.tsv", index_col=0)

    bar_df = dx_count_data.melt(id_vars="Bam", value_vars=["Duplex", "Orphan simplex"], var_name="Read type", value_name="Read count")

    bar_df["Proportion"] = bar_df["Read count"] / bar_df.groupby("Bam")["Read count"].transform(sum)
    bar_df = bar_df.replace(["CBM_2_rep1", "CBM_2_rep2", "CBM_3_rep1", "CBM_3_rep2", "zymo_wga_methylated_rep1", "zymo_wga_methylated_rep2", "zymo_wga_unmodified_rep1", "zymo_wga_unmodified_rep2", "Orphan simplex"],
                            ["CBM2_1", "CBM2_2", "CBM3_1", "CBM3_2", "ZM1", "ZM2", "ZU1", "ZU2", "Simplex only"])

    root_path = "data/duplex_data/patterns/"
    files = ["CBM_2_rep1.masked.bed.duplex_patterns.tsv", "CBM_3_rep1.sorted.bam.bed.duplex_patterns.tsv",
             "CBM_2_rep2.masked.bed.duplex_patterns.tsv", "CBM_3_rep2.sorted.bam.bed.duplex_patterns.tsv"]

    file_paths = [root_path + file for file in files]

    patterns = pd.concat([pd.read_table(path) for path in file_paths])

    mpl.rc('font', size=5)

    fig = plt.figure(dpi=600, figsize=(89/25.4, 50/25.4), layout="constrained")
    gs = fig.add_gridspec(2, 3)

    ax1 = fig.add_subplot(gs[0, :])

    sns.barplot(bar_df, 
                x="Bam", y="Read count",
                hue="Read type", 
                palette="BuPu",
                ax=ax1)

    sns.move_legend(ax1, "upper right", ncols=2, title=None, frameon=False)
    ax1.set_xlabel(None)

    palette = {"C" : "#edf8b1",
                "5mC" : "#7fcdbb",
                "5hmC" : "#2c7fb8"}

    ax2, ax3, ax4 = [fig.add_subplot(gs[1, i]) for i in range(3)]

    for ax, mod in zip([ax2, ax3, ax4], ["-", "m", "h"]):
        make_pie(patterns, mod, ax, palette.values())

    [ax.set_title(ascii_lowercase[i], fontweight="bold", loc="left") for i, ax in enumerate(fig.axes)]
    [ax.set_title(title) for ax, title in zip([ax2, ax3, ax4], ["C", "5mC", "5hmC"])]

    sns.despine()
    fig.savefig("plots/duplex_sfig.png")

if __name__=="__main__":
    main()