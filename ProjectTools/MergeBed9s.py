import pandas as pd
import argparse

def load_files(paths):
    
    headers=["chromosome", "chromStart", "chromEnd", "name", "score", "strand",
         "thickStart", "thickEnd", "reserved", "readCount", "percentMeth", "refContext",
         "calledContext", "genotypeQual"]
    
    bed_dict = {}

    for path in paths:
        bed = pd.read_csv(path, sep="\t", names=headers, header=None, index_col=None)
        bed.drop(["name", "score", "thickStart", "thickEnd", "reserved", "refContext", "calledContext", "genotypeQual"], 
                 axis=1, inplace=True)
        bed_dict.update({path : bed})
    
    return bed_dict

def merge_bed_dict_to_df(bed_dict):
    concatenated_beds = pd.concat(bed_dict.values()).reset_index(drop=True)
    
    read_counts = concatenated_beds.groupby(["chromosome", "chromStart", "chromEnd", "strand"], as_index=False)["readCount"].sum()
    percent_meth = concatenated_beds.groupby(["chromosome", "chromStart", "chromEnd", "strand"], as_index=False)["percentMeth"].mean()

    merged_df = pd.merge(read_counts, percent_meth)
    
    return merged_df

if __name__=="__main__":
    parser = argparse.ArgumentParser(
                        prog = "MergeBed9s",
                        description = "Merges two or more modified base-containing bed files (bed9) into a 6-column tab-separated output file.")
    parser.add_argument('paths', action="store", type=str, nargs='+', help='Filepaths for each bed9 file')
    parser.add_argument("-o ", "--output_path", action="store", dest="out_path", default="./merged_bed.tsv") 

    args = parser.parse_args()
    paths = args.paths
    out_path = args.out_path
    
    print(f"Loading {paths}")
    bed_dict = load_files(paths)
    
    print(f"Merging together")
    df = merge_bed_dict_to_df(bed_df)
    
    print(f"Saving to {args.paths}")
    df.to_csv(out_path, sep="\t", index=False)