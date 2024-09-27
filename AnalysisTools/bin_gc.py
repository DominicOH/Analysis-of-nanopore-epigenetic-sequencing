import pandas as pd
import argparse
import numpy as np

def main(file: str):
    names = ["chromosome", "start", "end", "modType", "readCount", "modCalls", "GC%"]
    load_file = pd.read_table(file, sep="\t", names=names)

    df_pivot = load_file.pivot(index=["chromosome", "start", "end", "readCount", "GC%"], columns="modType", values="modCalls").reset_index()
    del load_file
    df_pivot["GCBin"] = pd.cut(df_pivot["GC%"], bins=np.arange(0, 105, 5), labels=np.arange(5, 105, 5))
    gb = (df_pivot.groupby("GCBin", as_index=False)[["readCount", "m", "h"]].sum(numeric_only=True)
          .assign(FPR = lambda r: (r["m"] + r["h"])/r["readCount"]))

    return gb.to_csv(f"{file}.gcBinned", "\t", header=False, index=False)

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        prog = "bin_GC",
        description = "Reads a bedfile of modified base information and GC%. Bins modified base calls by GC%")
    parser.add_argument("target_file")

    args = parser.parse_args()
    main(args.target_file)