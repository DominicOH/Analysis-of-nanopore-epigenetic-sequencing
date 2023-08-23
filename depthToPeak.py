import csv
import argparse

def readSamtoolsDepth(inpath, outpath):
    with open(inpath, newline="", mode="rt") as depth_file, open(outpath, "wt") as new_file:
        reader = csv.reader(depth_file, delimiter="\t")
        writer = csv.writer(new_file, delimiter="\t")

        row = next(reader)

        chrom = row[0]
        start = row[1]
        end = int(start) + 1
        depth = row[2]
        peak_counter = 0

        for row in reader:
            if int(row[1]) in range(int(end) - 10, int(end) + 10):
                end = row[1]
                if row[2] > depth:
                    depth = row[2]
                current_peak = [chrom, start, end, depth]
            else: 
                writer.writerow(current_peak)
                current_peak += 1
                chrom = row[0]
                start = row[1]
                end = int(start) + 1
                depth = row[2]

    return print(f"Finished writing {peak_counter} peaks.")

def readBedtoolsCoverage(inpath, outpath):
    with open(inpath, newline="", mode="rt") as depth_file, open(outpath, "wt") as new_file:
        reader = csv.reader(depth_file, delimiter="\t")
        writer = csv.writer(new_file, delimiter="\t")
        row = next(reader)

        chrom = row[0]
        start = row[1]
        end = row[2]
        depth = row[3]
        current_peak = [chrom, start, end, depth]
        peak_counter = 0

        for row in reader:        
            if int(row[1]) in range(int(end) - 10, int(end) + 10): # if the start of the next line is approximately equal to the end of the previous
                end = row[2] # the end of the peak is updated
                start = start
                if row[3] > depth:
                    depth = row[3] # the highest depth is taken
                current_peak = [chrom, start, end, depth]
            else:
                writer.writerow(current_peak)
                peak_counter += 1
                chrom = row[0]
                start = row[1]
                end = row[2]
                depth = row[3]
                current_peak = [chrom, start, end, depth]

    return print(f"Finished writing {peak_counter} peaks.")

if __name__=="__main__":
    parser = argparse.ArgumentParser(
                        prog = "depthToPeak",
                        description = "Reads the output of a depth/coverage file and outputs joined peaks.")
    parser.add_argument("input_tool", choices=["samtools", "bedtools"], action="store", dest="choice", required=True) 
    parser.add_argument("-i ", "--input_path", action="store", dest="outpath", required=True) 
    parser.add_argument("-o ", "--output_path", action="store", dest="inpath", required=True) 

    args = parser.parse_args()

    inpath = args.inpath
    outpath = args.outpath

    if args.choice == "samtools": 
        readSamtoolsDepth(inpath, outpath)
    elif args.choice == "bedtools":
        readBedtoolsCoverage(inpath, outpath)
    
    

