import subprocess
import pandas as pd

class Reference:
    """
    Objects used as genomic or feature references. May be in the form of tab-separated variables (TSV) files.
    """
    def __init__(self, path=None):
        self.path = path
    
    def checkNumColumns(self):
        """
        Checks and returns the number of columns present in the TSV file.
        """
        first_line = subprocess.check_output(["head", "-n 1", f"{self.path}"]).decode("utf-8")
        num_columns = len(first_line.split("\t"))

        return num_columns
    
    def getColNames(self):
        """
        Uses the number of columns to predict column name labels. 
        """
        names = ["Chromosome", "Start", "End"]
        return names
    
    def getFirstLine(self):
        first_line = subprocess.check_output(["head", "-n 1", f"{self.path}"]).decode("utf-8").split("\t")
        return first_line

class Features(Reference):
    """
    Regularly used object type containing information about genomic features. May be a tab-separated variables (TSV) file. Columns must contain follow BED 'Chromosome', 'Start', 'End', 'Name' format. 
    """

    def __init__(self, path, dataframe=None):
        super().__init__(path)
        self.dataframe = dataframe
   
    def getColNames(self):
        """
        Uses the number of columns to predict column name labels. 
        """
        col_length = super().checkNumColumns()
        names =  super().getColNames()
        
        if col_length == 4:
            names.extend(["Name"])
        elif col_length == 5:
            names.extend(["Name", "Strand"])
        elif col_length == 6:
            names.extend(["Name", "Score", "Strand"])
        else:
            raise ValueError("Please check number of columns in file. Must equal 4, 5, or 6.")
        return names
    
    def retrieveFeatureType(self):
        """
        Uses the file extension to determine the type of feature. File must be saved in a './1_2_3_feature-type_4.bed' fashion. 
        """
        filename = self.path.split("/").pop(-1)
        feature_type = filename.split("_").pop(3)

        try:
            return feature_type.strip(".bed")
        except:
            return feature_type
    
    def toDF(self):
        """
        Shows the feature file as a pandas DataFrame.
        """
        if not self.dataframe:
            df = pd.read_csv(self.path, sep="\t", names=self.getColNames())
            df["feature_type"] = self.retrieveFeatureType()
            self.dataframe = df
        return df

class CGIs(Reference):
    """
    Subclass of reference tab-separated variable file containing information about CpG island positions. file
    """
    def __init__(self, path):
        self.path = path

    def retrieveFeatureType(self):
        first_line = super().getFirstLine()
        
        return first_line[-1].strip("\n")
    
    def getColNames(self):
        """
        Uses the number of columns to predict column name labels. 
        """
        col_length = super().checkNumColumns()
        names =  super().getColNames()
        first_line = super().getFirstLine()

        if first_line[3].__contains__("CpG"):
            names.extend(["NCpGs", "feature_type"])
        else:
            names.extend(["Name", "feature_type"])
        
        return names
    
    def toDF(self):
        """
        Shows the feature file as a pandas DataFrame.
        """
        df = pd.read_csv(self.path, sep="\t", names=self.getColNames())
        
        return df
    
