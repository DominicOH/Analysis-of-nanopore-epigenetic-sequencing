import argparse
import pandas as pd
from AnalysisTools.helpers import timer
import concurrent.futures

def read_file(file):
    print(f"Reading file {file}")
    
    with open(file, "r") as open_file:
        asymm_states = ["-,m,C", "m,-,C", "m,h,C", "h,m,C", "-,h,C", "h,-,C"]
        chromnames, starts, ends, asym_reads = [], [], [], []

        for line in open_file.readlines():
            readline = line.split("\t")
            modstate = readline[3].strip()
            if modstate in asymm_states:
                if len(starts) > 0:
                    if readline[1] == starts[-1]:
                        asym_reads[-1] += int(readline[11])
                    else:
                        chromnames.append(readline[0])
                        starts.append(readline[1])
                        ends.append(readline[2])
                        asym_reads.append(int(readline[11])) 
                else:
                    chromnames.append(readline[0])
                    starts.append(readline[1])
                    ends.append(readline[2])
                    asym_reads.append(int(readline[11])) 
            else:
                continue
        
        file_df = pd.DataFrame({"Chromosome" : chromnames,
                                "Start" : starts,
                                "End" : ends,
                                "Asym_reads" : asym_reads})
        
        print(f"Found {sum(file_df['Asym_reads'])} asymmetrical sites across {len(file_df)} positions")

    return file_df
        
@timer
def main(filenames, outpath):
    with concurrent.futures.ThreadPoolExecutor(len(filenames)) as read_executor:
        all_files = read_executor.map(read_file, filenames)
        all_files = [file for file in all_files]

    print("Merging files")
    all_files_df = (pd.concat(all_files)
                    .groupby(["Chromosome", "Start", "End"], observed=True, as_index=True)
                    .sum(numeric_only=True)
                    .reset_index())
    
    print(f"Saving {len(all_files_df)} lines to {outpath}")
    return all_files_df.to_csv(outpath, sep="\t", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                        prog = "filter_hemi",
                        description = "Compares the overlap of duplex modified bases, oxBS-seq and TAB-seq.")
    parser.add_argument("filenames", action="store", nargs="+")
    parser.add_argument("-o", "--outpath", action="store", type=str, required=True)

    args = parser.parse_args()
    main(args.filenames, args.outpath)  
