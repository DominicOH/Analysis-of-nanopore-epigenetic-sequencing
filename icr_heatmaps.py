import argparse
from AnalysisTools import fetch_reads_from_modkit
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from matplotlib import cm
from AnalysisTools.helpers import timer

@timer
def fig_main(extract_path, include_bed, outpath, width, height):
    merged_read_table = fetch_reads_from_modkit.ModkitExtract(extract_path).cpg_table
    print("Loading modkit extract data")
    merged_read_table.set_include_bed(include_bed)

    fig, axes = plt.subplots(5, 3,
                            dpi=600, 
                            figsize=(width/25.4, height/25.4), 
                            layout="constrained")

    mpl.rc(("fontsize", 5))

    axes = axes.flatten()

    cmap = (mpl.colors.ListedColormap(sns.color_palette("YlGnBu", 3)))
    norm = mpl.colors.BoundaryNorm([0, 1, 2, 3], cmap.N)
    cbar_ax = fig.add_axes([1.01, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cbar_ax)
    cbar.minorticks_off()
    cbar.set_ticks([.5, 1.5, 2.5], labels=["C", "5mC", "5hmC"])
    cbar_ax.set_in_layout(True)

    gene_table = pd.read_table("feature_references/dmr/mm39_dmr_coordinates_modified.bed", 
                            names=["Chromosome", "Start", "End", "Name", "Allele"])["Name"].unique()

    for i, name in enumerate(gene_table):
        gene_table = merged_read_table.select_gene(name)
        gene_table.heatmap(min_cpg_proportion=.25, minimum_read_proportion=.25, ax=axes[i])
        print(f"Completed plotting {name}")

    fig.savefig(outpath)

##### main function ##### 

if __name__=="__main__":
    parser = argparse.ArgumentParser(
                        prog = "icr_heatmaps",
                        description = "Produce a figure of heatmaps, across multiple DMRs, with reads split by allele.")
    parser.add_argument("extract_path", 
                        help="Path to modkit extract output.")
    parser.add_argument("include_bed", 
                        help="Path to DMR bed file.")
    parser.add_argument("-w ", "--width", 
                        action="store", 
                        dest="width", 
                        default=160,
                        help="Width of output image (mm).") 
    parser.add_argument("-h ", "--height", 
                        action="store", 
                        dest="height", 
                        default=160,
                        help="Height of output image (mm).")
    parser.add_argument("-o ", "--outpath", 
                        action="store", 
                        dest="outpath", 
                        required=True,
                        help="Filename for output.")
    args = parser.parse_args()

    fig_main(args.extract_path, args.include_bed, args.outpath, args.width, args.height)