import argparse
import pandas as pd
import concurrent.futures
import os
from helpers import timer

def count_patterns(path):
    print(f"Reading patterns in {path}...")
    pattern_df = pd.read_table(path, sep="\t", 
                            names=["Chromosome", "Start", "End", "Pattern", "readCount", "D0", "D1", "D2", "D3", "D4", "percentPattern", "N_Pattern", "N_Canonical", "N_Other", "D5", "D6", "D7", "D8"],
                            usecols=["Chromosome", "Start", "End", "Pattern", "readCount", "N_Pattern"])
    
    pattern_count = pattern_df.groupby("Pattern")["N_Pattern"].sum().reset_index()

    outpath = outdir + os.path.basename(path) + ".duplex_patterns.tsv"

    print(f"Finished counting patterns at {len(pattern_df)} positions. Saving to {outpath}.")
    return pattern_count.to_csv(outpath, sep="\t", index=False)

@timer
def main(paths):
    if type(paths) == list:
        with concurrent.futures.ThreadPoolExecutor(len(paths)) as tpe:
            multiple_hemi_pileups = tpe.map(count_patterns, paths)
            return multiple_hemi_pileups
    else:
        return count_patterns(paths)

##### Main function #####

if __name__=="__main__":
    parser = argparse.ArgumentParser(
                        prog = "count_reads_dx",
                        description = "Reads one or more hemi-pileup files and tallies CpG modification states.")
    parser.add_argument("filenames", nargs="+")
    parser.add_argument("-o ", "--outdir", type=str, required=False) 

    args = parser.parse_args()

    if args.outdir:
        global outdir
        outdir = args.outdir
    
    main(args.filenames)