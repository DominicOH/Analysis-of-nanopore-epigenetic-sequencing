from modbampy import ModBam
import pandas as pd
import ImprintPhasing

def CpGsFromModbam(path : str, 
                   regions_of_interest: dict
                   ):
    """
    Extracts a pandas DataFrame table of modified basecalls from the target Nanopore modbam file using ONT's modbampy package. 
    Coordinates are selected using a dictionary object containing the name of the target region and its chromosomal position.

    :param str path: Path of the target Nanopore modbam file. File must be sorted and indexed using samtools index. 
    :param dict regions_of_interest: Dictionary with keys/values like so:
                                                    {0: {'Chromosome': 'chr2', # str
                                                    'Start': 174136155, # int 
                                                    'End': 174142856, # int
                                                    'Name': 'Nespas-Gnasxl', # any
                                                    } ... 
                                                    } 
    """

    with ModBam(path) as bam:
        all_regions = []
        for region in regions_of_interest.values():
            for read in bam.reads(region["Chromosome"], region["Start"], region["End"]):
                [all_regions.append([region["Name"], region["Chromosome"], region["Start"], region["End"], read.reference_start, read.reference_end, *pos_mod]) for pos_mod in read.mod_sites]

    region_df = pd.DataFrame(all_regions).rename(columns={
        0 : "Name",
        1 : "Chromosome",
        2 : "regionStart",
        3 : "regionEnd",
        4 : "readStart",
        5 : "readEnd", 
        6 : "readID",
        7 : "refPos", 
        8 : "qPos",                                                            
        9 : "strand", 
        10 : ".",
        11 : "cBase",
        12 : "modBase", 
        13 : "modScore"}).drop(columns=[".", "cBase"])
    return ImprintPhasing.ModBaseAssemblage(region_df)
    
def selectCpGsFromModBam(path: str, 
                         chrom: str, 
                         start: int, 
                         end: int
                         ):
    """
    Extracts a pandas DataFrame table of modified basecalls from the target Nanopore modbam file using ONT's modbampy package. Target coordinates are supplied as arguments. 
    
    :param str path: Path of the target Nanopore modbam file. File must be sorted and indexed using samtools index. 
    :param str chrom: Chromosome of target region in standard format; eg.: "chr2"
    :param int start: Start position of targe region; eg. 174136155
    :param int end: End position of targe region; eg. 174142856
    """

    with ModBam(path) as bam:
        read_list = []
        for read in bam.reads(chrom, start, end):
            [read_list.append([chrom, start, end, read.reference_start, read.reference_end, *pos_mod]) for pos_mod in read.mod_sites]

    read_df = pd.DataFrame(read_list).rename(columns={
        0 : "Chromosome",
        1 : "regionStart",
        2 : "regionEnd",
        3 : "readStart",
        4 : "readEnd", 
        5 : "readID",
        6 : "refPos", 
        7 : "qPos",                                                            
        8 : "strand", 
        9 : ".",
        10 : "cBase",
        11 : "modBase", 
        12 : "modScore"}).drop(columns=[".", "cBase"])
    return ImprintPhasing.ModBaseAssemblage(read_df)

