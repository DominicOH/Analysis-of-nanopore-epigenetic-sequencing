import subprocess
import pandas as pd
import pyranges as pr
import concurrent.futures
import numpy as np

class Reference:
    """
    Objects used as genomic or feature references extracted from files. Input files should be in BED4, BED6, BED8, and BED12 format files are supported.
    """
    def __init__(self, 
                 path: str=None):
        self._path = path
        self.pr = pr.read_bed(self._path)
        self._df = None

    @property
    def path(self): # getter function for read-only path value
        return self._path
    
    def __get_feature_type(self):
        """
        Uses the file extension to determine the type of feature. Feature type must be stored prior to the file extension within the file name. 
        """
        filename = self.path.split("/").pop(-1)
        feature_type = filename.split("_").pop(-1).split(".").pop(0)
        
        return feature_type
    
    def __merge_overlaps(self):
        pr = self.pr
        merged_pr = pr.merge()
        return merged_pr 
        
    @property
    def df(self):
        if self._df == None:
            pr = self.__merge_overlaps()
            df = pr.as_df()
            
            df["feature_type"] = self.__get_feature_type()
            self._df = df
            return df
        else: 
            return self._df
        
def fetch_feature_PyRange(dir_path: str, p_threads=1):
    """
    Takes a directory of BED4 or BED6 files containing lists of features and feature coordinates to be used for annotation purposes. 
    """
    gene_feature_list = subprocess.check_output(["ls", dir_path]).decode("utf-8").split("\n") 
    gene_feature_list.pop(-1) # removes the current directory dot node 

    def add_reference(filepath):
        path = dir_path + filepath
        feature_tsv = Reference(path)
        return feature_tsv.df
    
    if p_threads > 1:
        with concurrent.futures.ThreadPoolExecutor(6) as tpe:
            df_futures = tpe.map(add_reference, gene_feature_list)
            feature_reference_df = pd.concat([df for df in df_futures]).drop(columns=["Score", "ThickStart", "ThickEnd"], errors="ignore")
    else:
        df_generator = map(add_reference, gene_feature_list)
        feature_reference_df = pd.concat(df_generator).drop(columns=["Score", "ThickStart", "ThickEnd"], errors="ignore")

    return pr.PyRanges(feature_reference_df)

def annotate(df, dir_path):
    feature_pr = fetch_feature_PyRange(dir_path)

    annotated_df = pr.PyRanges(df).join(feature_pr, strandedness=False, suffix="_Feature", apply_strand_suffix=False).as_df()
    return annotated_df
        