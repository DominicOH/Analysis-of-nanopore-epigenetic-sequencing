import argparse
import pandas as pd
import upsetplot
import matplotlib.pyplot as plt
from AnalysisTools.helpers import timer

@timer
def make_upset(path, w, h, outpath):
    """
    Plots an upsetplot of Nanopore duplex modiifcation states. Requires inputs in the format produced by count_duplex_pattern.py"
    """

    print("Loading data")
    if type(path) == str:
        patterns = pd.read_table(path)
    elif type(path) == list:
        patterns = pd.concat([pd.read_table(p) for p in path])

    print("Loaded data")
    patterns_grouped = patterns.groupby("Pattern")["N_Pattern"].sum()

    upsetplot_reads = upsetplot.from_memberships(
        memberships=[
    ["C"], 
    ["5mC"],
    ["C", "5mC"],
    ["5hmC"],
    ["5mC", "5hmC"],
    ["C", "5hmC"]],
        data = [
        patterns_grouped["-,-,C"], 
        patterns_grouped["m,m,C"], 
        patterns_grouped["-,m,C"] + patterns_grouped["m,-,C"] , 
        patterns_grouped["h,h,C"], 
        patterns_grouped["h,m,C"] + patterns_grouped["m,h,C"], 
        patterns_grouped["-,h,C"] + patterns_grouped["h,-,C"]] 
    )
    upset_reads = upsetplot.UpSet(upsetplot_reads, sort_by="input", sort_categories_by="input", show_percentages=True,
                                    element_size=None,
                                    facecolor="grey",
                                    intersection_plot_elements=3)

    upset_fig = plt.figure(figsize=(w/25.4, h/25.4), dpi=600)
    print("Plotting upset")
    upset_reads.plot(fig=upset_fig)

    return upset_fig.savefig(outpath, dpi=600)

if __name__=="__main__":
    parser = argparse.ArgumentParser(
                        prog = "compare_hmedip",
                        description = "Produces an upsetplot of nanopore duplex modification states. Requires output of count_duplex_pattern.py.")
    parser.add_argument("filenames", nargs="+", help="Filepaths for tables produced by count_duplex_pattern.py")
    parser.add_argument("-w ", "--width", 
                        action="store", 
                        dest="width", 
                        default=89,
                        help="Width of output image (mm).") 
    parser.add_argument("-h ", "--height", 
                        action="store", 
                        dest="height", 
                        default=89,
                        help="Height of output image (mm).")
    parser.add_argument("-o ", "--outpath", 
                        action="store", 
                        dest="outpath", 
                        required=True,
                        help="Filename for output.")
    
    args = parser.parse_args()
    make_upset(args.filenames, args.width, args.height, args.outpath)