
import pandas as pd
from pathlib import Path
from AnalysisTools import fetch_reads_from_modkit
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

def main():
    path_base = "data/duplex_data/dmr/readcalls"
    files = [file for file in Path(path_base).glob("*")]
    include_bed = "feature_references/dmr/mm39_dmr_coordinates_modified.bed"
    all_replicates = []

    for i, file in enumerate(files):
        dup_read_table = fetch_reads_from_modkit.ModkitExtract(file).cpg_table
        dup_read_table.set_include_bed(include_bed)

        all_genes = []
        m_reads_all = []
        u_reads_all = []
        for gene in dup_read_table.include_bed["Name"]:
            try:
                m_reads = len(dup_read_table.select_gene(gene).get_methylated_readIDs(minimum_read_proportion=.25, min_cpg_proportion=.25))
                u_reads = len(dup_read_table.select_gene(gene).get_unmethylated_readIDs(minimum_read_proportion=.25, min_cpg_proportion=.25))
                gene_df = pd.DataFrame(dict(gene_name = gene, read_type = ["Methylated", "Unmethylated"], count=[m_reads, u_reads]))
                all_genes.append(gene_df)
            except:
                print(f"Failed at {gene} in replicate {i}")
                pass

        all_genes_df = pd.concat(all_genes).assign(Replicate = i)
        all_replicates.append(all_genes_df)

    for read_ls, allele in zip([m_reads_all, u_reads_all], ["methylated", "unmethylated"]):
        reads = pd.concat(read_ls)
        reads.to_csv(f"data/duplex_data/duplex_icrs/{allele}_reads.txt", index=False, header=False)
    all_replicates_df = pd.concat(all_replicates).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(180/25.4, 2), dpi=600, layout="constrained")

    mpl.rc(("size", 5))

    sns.barplot(all_replicates_df, 
                x="gene_name",
                y="count", 
                palette="GnBu",
                hue="read_type",
                hue_order=["Unmethylated", "Methylated"],
                err_kws=dict(linewidth=0.8),
                capsize=.5,
                ax=ax)

    ax.set_ylabel("Count of reads")
    ax.set_xlabel(None)
    ax.tick_params("x", labelrotation=30)
    sns.move_legend(ax, "upper center", ncols=2, frameon=False, title="Allele", bbox_to_anchor=(.5, 1.5))
    sns.despine()
    
    return fig.savefig("plots/all_icrs_readcount.svg")

if __name__ == "__main__":
    main()