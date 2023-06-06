import subprocess
import pandas as pd

class Checkpoint:
    """
    Class to support saving and loading dataframes after intensive data processing steps.  
    """
    def __init__(self, dataframe=None, name=None, path=None):
        self.dataframe = dataframe
        self.name = name
        self.path = path

    def saveCheckpoint(self, name=None, path=None):
        """
        Saves the Checkpoint object as a TSV file in a folder for intermediate data files.
        """
        if self.dataframe.empty:
            raise ValueError("Missing dataframe")
        
        if self.name:
            name = self.name
        else: 
            name = name
            if not name:
                raise ValueError("No name")

        if not path:
            path = self.path
            if not path:
                path = f"./intermediates/{name}.tsv"
                self.path = path

        print(f"Checkpointing {self.name}.")
        self.dataframe.to_csv(path, sep="\t", index=False, header=False)
        print(f"Saved to {self.path}.")
        return 
    
    def getHeaders(self):
        """
        Uses number of columns to estimate column names - based on standards in current analysis.
        """
        colnames = ["chromosome", "chromStart", "chromEnd"]
        first_line = subprocess.check_output(["head", "-n 1", f"{self.path}"]).decode("utf-8").split("\t")

        if len(first_line) == 6: 
            # pyrange
            colnames = (["Chromosome", "Start", "End", "Strand", "percentMeth_Nanopore_5hmC", "percentMeth_Bisulphite_5hmC"])
        elif len(first_line) == 7:
            # long dataframe
            colnames.extend(["method", "strand", "readCount", "percentMeth"])
        elif len(first_line) == 8:
            if self.path.__contains__("two_mod"): 
                colnames.extend(["strand", "readCount_WGBS", "percentMeth_WGBS", "readCount_Nanopore", "percentMeth_Nanopore"])
            elif self.path.__contains__("mc") and not self.path.__contains__("hmc"):
                colnames.extend(["strand", "readCount_oxBS_5mC", "percentMeth_oxBS_5mC", "readCount_Nanopore_5mC", "percentMeth_Nanopore_5mC"])
            elif self.path.__contains__("hmc"):
                colnames.extend(["strand", "readCount_TAB_5hmC", "percentMeth_TAB_5hmC", "readCount_Nanopore_5hmC", "percentMeth_Nanopore_5hmC"])
            else:
                colnames.extend(["strand", "percentMeth_oxBS_5mC", "percentMeth_Nanopore_5mC", "percentMeth_TAB_5hmC", "percentMeth_Nanopore_5hmC"])
            
        else:
            raise ValueError("Check file headers.")
        return colnames
    
    def loadCheckpoint(self, colnames=None, path=None):
        """
        Loads the Checkpoint object TSV from either the attribute path or path provided.
        """
        if not path:
            if not self.path:
                raise ValueError("Missing required path string.")
            else: 
                path = self.path
        else:
            self.path = path

        if not colnames:
            colnames = self.getHeaders()

        dataframe = pd.read_csv(path, sep="\t", names=colnames)
        return dataframe