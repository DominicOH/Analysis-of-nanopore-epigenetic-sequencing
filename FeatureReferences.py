import subprocess
import pandas as pd

class Reference:
    """
    Objects used as genomic or feature references extracted from files. Input files are preferably in BED4 or BED6 format.
    """
    def __init__(self, 
                 path: str=None):
        self._path = path

    @property
    def path(self): # getter function for read-only path value
        return self._path
    
    def __check_number_of_columns(self):
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
        num_columns = self.__check_number_of_columns()
        names = ["Chromosome", "Start", "End", "Name"]

        if num_columns == 6:
            names.extend(["Score", "Strand"])
        elif num_columns == 8:
            names.extend(["Score", "Strand", "ThickStart", "ThickEnd"])
        return names
    
    def __get_feature_type(self):
        """
        Uses the file extension to determine the type of feature. Feature type must be stored prior to the file extension within the file name. 
        """
        filename = self.path.split("/").pop(-1)
        feature_type = filename.split("_").pop(-1)

        try:
            return feature_type.strip(".bed")
        except:
            return feature_type
        
    def as_dataframe(self):
        names = self.__get_column_names()
        dataframe = pd.read_csv(self.path, sep="\t", names=names)
        dataframe["feature_type"] = self.__get_feature_type()
        return dataframe