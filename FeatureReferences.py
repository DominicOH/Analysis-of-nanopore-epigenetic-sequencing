import subprocess
import pandas as pd
import pyranges as pr

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
    
    def __check_num_columns(self):
        """
        Checks and returns the number of columns present in the TSV file.
        """
        first_line = subprocess.check_output(["head", "-n 1", f"{self.path}"]).decode("utf-8")
        num_columns = len(first_line.split("\t"))

        return num_columns
    
    def __get_column_names(self):
        """
        Uses the number of columns to predict column name labels. 
        """
        num_columns = self.__check_num_columns()
        names = ["Chromosome", "Start", "End", "Name"]

        if num_columns == 6:
            names.extend(["Score", "Strand"])
        elif num_columns == 8:
            names.extend(["Score", "Strand", "ThickStart", "ThickEnd"])
        elif num_columns == 12:
            names.extend(["Score", "Strand", "ThickStart", "ThickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"])
        return names
    
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
            
            # Currently not working with names
            # names = self.__get_column_names()
            # df.columns = names
            
            df["feature_type"] = self.__get_feature_type()
            self._df = df
            return df
        else: 
            return self._df
        
def featureRefPyRange(dir_path: str):
    """
    Takes a directory of BED4 or BED6 files containing lists of features and feature coordinates to be used for annotation purposes. 
    """
    gene_feature_list = subprocess.check_output(["ls", dir_path]).decode("utf-8").split("\n") 
    gene_feature_list.pop(-1) # removes the current directory dot node 

    df_list = []
    for file in gene_feature_list:
        path = dir_path + file
        feature_tsv = Reference(path)
        feature_df = feature_tsv.df
        df_list.append(feature_df)

    feature_reference_df = pd.concat(df_list).drop(columns=["Score", "ThickStart", "ThickEnd"], errors="ignore")
    return pr.PyRanges(feature_reference_df)
        