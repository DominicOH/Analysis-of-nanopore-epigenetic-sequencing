import pandas as pd
import argparse

def make_ob_table(path):  
    header = ["read", "sense", "chromosome", "start", "modification_state"]
    ob_df = pd.read_table(path, names=header, skiprows=1)

    ob_df = ob_df.groupby(["chromosome", "start", "modification_state"]).size().reset_index(name="count")
    ob_df["bed_start"] = ob_df["start"].sub(2)
    ob_df["bed_end"] = ob_df["start"].sub(1)

    ob_df = ob_df.pivot(index=["chromosome", "bed_start", "bed_end"], columns="modification_state", values="count").fillna(0)
    ob_df = ob_df.rename(columns={"Z":"modified_reads", "z":"unmodified_reads"})

    return ob_df

def make_ot_table(path): 
    header = ["read", "sense", "chromosome", "start", "modification_state"]
    ot_df = pd.read_table(path, names=header, skiprows=1)

    ot_df = ot_df.groupby(["chromosome", "start", "modification_state"]).size().reset_index(name="count")
    ot_df["bed_start"] = ot_df["start"].sub(1)
    ot_df["bed_end"] = ot_df["start"]

    ot_df = ot_df.pivot(index=["chromosome", "bed_start", "bed_end"], columns="modification_state", values="count").fillna(0)
    ot_df = ot_df.rename(columns={"Z":"modified_reads", "z":"unmodified_reads"})

    return ot_df

def merge_strands(ot_path, ob_path):
    ot_table = make_ot_table(ot_path)
    ob_table = make_ot_table(ob_path)

    return pd.concat([ot_table, ob_table])

def output_bed(table, path):
    return table.to_csv(path, sep="\t", header=False)

if __name__=="__main__":
    parser = argparse.ArgumentParser(
                        prog = "BismarkToBed",
                        description = "Converts Bismark output files to BED file format and merges separated strand files.")
    parser.add_argument("-t", "--original_top", required=True, action="store", dest="ot_path")
    parser.add_argument("-b", "--original_bottom", required=True, action="store", dest="ob_path")
    parser.add_argument("-o ", "--output_path", action="store", dest="out_path", default="./BismarkBedOut.bed") 

    args = parser.parse_args()
    print(args.ot_path, args.ob_path, args.out_path)
    
    ot_path = args.ot_path
    ob_path = args.ob_path
    out_path = args.out_path
    
    bed_output = merge_strands(ot_path, ob_path)
    output_bed(bed_output, out_path)
    
    
    