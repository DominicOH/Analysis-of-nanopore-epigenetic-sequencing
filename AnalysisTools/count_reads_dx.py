import argparse
from helpers import timer
import pysam as ps
import concurrent.futures
import pandas as pd
import os

@timer
def count_dx(bam):
    read_counter = 0
    duplex_reads = 0
    orphan_simplex = 0

    with ps.AlignmentFile(bam, "rb") as af:
        for read in af.fetch():
            dx = read.get_tag("dx")
            read_counter+=1

            if dx == 1:
                duplex_reads+=1
            elif dx == -1:
                orphan_simplex+=1

        print(f"{os.path.basename(bam).split('.')[0]}. Duplex: {duplex_reads}. Orphan simplex: {orphan_simplex}.")

        result_dict = {
            "Bam" : os.path.basename(bam.split(".")[0]), 
            "Duplex" : duplex_reads,
            "Orphan simplex" : orphan_simplex,
            "Total reads" : read_counter
        }

    return result_dict

def pool_jobs(paths):
    with concurrent.futures.ThreadPoolExecutor(len(paths)) as executor:
        count_futures = executor.map(count_dx, paths)

    results = [result for result in count_futures]
    return results

@timer
def main(paths, outpath):
    results = pool_jobs(paths)

    result_df = pd.DataFrame().from_dict(results, "columns")
    return result_df.to_csv(outpath, sep="\t")

##### Main function #####

if __name__=="__main__":
    parser = argparse.ArgumentParser(
                        prog = "count_reads_dx",
                        description = "Reads one or more bamfiles and counts the values associated with the DX tag.")
    parser.add_argument("filenames", nargs="+")
    parser.add_argument("-o ", "--outpath", type=str, required=True) 

    args = parser.parse_args()

    main(args.filenames, args.outpath)
    print("Done")