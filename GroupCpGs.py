import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pyranges as pr
from FeatureReferences import *
from CpGIntersects import CpGIntersects
from scipy import stats
   
class GroupedDF:
    """
    Dataframe-type objects where CpG positions are grouped. Child classes contain additional functionality. Contains interface for relevant Seaborn plotting functions.  
    """
    def __init__(self, df, cpg_threshold=None):
        self.df = df.loc[df.loc[:, "CpG_count"].ge(cpg_threshold)]
    
    def __ratioToMean(self, column: str):
        new_col = self.df[column].divide(self.df[column].mean())
        return new_col
    
    def __log2RatioWrapper(self, column: str, include_zeros=False):
        if include_zeros:
            add = 1
        else: 
            add = 0
        
        with np.errstate(divide="ignore"):
            log2_col = np.log2(
                np.add(self.__ratioToMean(column), add))
        return log2_col

    def enrichmentComparison(self, include_zeros=False):
        """
        Provides additional columns with log_2 scale scores showing enrichment relative to the mean for 5mC and 5hmC. 

        :param bool include_zeros: Whether groups with an average CpG modification of zero are kept. Adds 1 to the ratio calculation to avoid zero division. 
        """
        df = self.df.copy()
        df = df.assign(
            log2enrichment_5mC_Min=self.__log2RatioWrapper("percentMeth_5mC_Min", include_zeros), 
            log2enrichment_5mC_Bisulphite=self.__log2RatioWrapper("percentMeth_5mC_Bisulphite", include_zeros), 
            log2enrichment_5hmC_Min=self.__log2RatioWrapper("percentMeth_5hmC_Min", include_zeros),            
            log2enrichment_5hmC_Bisulphite=self.__log2RatioWrapper("percentMeth_5hmC_Bisulphite", include_zeros)
            )
        
        return df

    def methodComparison(self):
        """
        Adds "Average" and "Difference" to the dataframe, displaying the average level of enrichment and difference between method enrichment levels respectively.
        
        TO BE DEPRECATED/REPAIRED
        """
        df = self.enrichmentComparison()

        df["Average"] = df[["Log2FromMean_TAB", "Log2FromMean_Nanopore"]].mean(axis=1)
        df["Difference"] = np.subtract(df["Log2FromMean_TAB"], df["Log2FromMean_Nanopore"])
        
        return df
        
class FeatureAndGene(GroupedDF):
    """
    Dataframe-like objects where CpG sites are grouped by gene, genomic feature, or CpG island. 
    """
    def __init__(self, df, cpg_threshold=None):
        super().__init__(df)
        self.cpg_threshold = cpg_threshold

    def asLongDf(self):
        """
        Converts the DF from a wide-form to a longer form. 
        """
        cdf = self.df.copy()

        if "percentMeth_TAB_5hmC" in cdf.columns:
            cdf.rename(columns={"percentMeth_TAB_5hmC" : "percentMeth_TAB"}, inplace=True)

        stubs = ["percentMeth"]
        indices = ["Name", "group_start", "group_end"]
        
        return pd.wide_to_long(cdf, stubs, indices, "method", sep="_", suffix="\D+").reset_index()
    
class tiledGroup(GroupedDF):
    """
    Dataframe-like objects where CpG sites are grouped by genomic window or tile. 

    TO BE DEPRECATED/REPAIRED.
    """
    def __init__(self, df, cpg_threshold):
        super().__init__(df, cpg_threshold)

    def positiveControlGroupDF(self, number_target_tiles):
        """
        Returns the positive control group - entries enriched in both methods.
        """
        df = super().methodComparison()
        df = df.nlargest(50, "Average")
        df = df.nsmallest(number_target_tiles, "Difference")
        return tiledGroup(df, self.cpg_threshold)
    
    def NegativeControlGroupDF(self, number_target_tiles):
        """
        Returns the negative control group - entries not enriched in either method.
        """
        df = super().methodComparison()
        df = df.nsmallest(50, "Average")
        df = df.nsmallest(number_target_tiles, "Difference")
        return tiledGroup(df, self.cpg_threshold)

    def NanoporePositiveGroupDF(self, number_target_tiles):
        """
        Returns the Nanopore positive test group - entries enriched only in Nanopore.
        """
        df = super().methodComparison()
        df = df.loc[df["Log2FromMean_TAB"] <= 0]
        df = df.nlargest(number_target_tiles, "Log2FromMean_Nanopore")
        return tiledGroup(df, self.cpg_threshold)
    
    def TabPositiveGroupDF(self, number_target_tiles):
        """
        Returns the TAB positive test group - entries enriched only in TAB.
        """
        df = super().methodComparison()
        df = df = df.loc[df["Log2FromMean_Nanopore"] <= 0]
        df = df.nlargest(number_target_tiles, "Log2FromMean_TAB")
        return tiledGroup(df, self.cpg_threshold)
    
    def getSequences(self):
        pyrange = pr.PyRanges(self.df)
        ref_fa = './data/Reference_data/mm39.fa'

        sequences = pr.get_sequence(pyrange, ref_fa)
        sequences.name = "sequence"

        pyrange = pyrange.insert(sequences)       
        return tiledGroup(pyrange.as_df(), self.cpg_threshold)
    
    def getIGVcoords(self):
        df = self.df
        coordinates = []
        for line, value in df.iterrows():
            chrom = value[0]
            s = str(value[1])
            e = str(value[2])
            first_half = ":".join([chrom, s])
            coord = "-".join([first_half, e])
            coordinates.append(coord)
        df["coordinates"] = coordinates
        return tiledGroup(df, self.cpg_threshold)
    
    def __reorderDF(self):
        df = self.df
        return df[["coordinates", "CpG_count", "percentMeth_Nanopore", "percentMeth_TAB", "Log2FromMean_TAB", "Log2FromMean_Nanopore", "Average", "Difference",  "sequence"]]

    def exportTests(self, number_target_tiles):
        wr = pd.ExcelWriter('/u/n/doh28/Documents/Nanopore_HMC/primer_design_regions.xlsx')

        list_of_groups = [self.NanoporePositiveGroupDF(number_target_tiles), self.TabPositiveGroupDF(number_target_tiles), self.positiveControlGroupDF(number_target_tiles), self.NegativeControlGroupDF(number_target_tiles)]
        processed_groups = [group.getSequences().getIGVcoords().__reorderDF() for group in list_of_groups]
        list_of_names = ["Nanopore_positive", "TAB_positive", "Positive_ctrl", "Negative_ctrl"]

        for i in np.arange(0, 4):
            processed_groups[i].to_excel(wr, list_of_names[i], index=False)

        wr.close()
        return 
    
    def tileWithLogCols(self):
        df = super().enrichmentComparison()
        return tiledGroup(df, self.cpg_threshold)
    
    def asCpGIntersect(self):
        df = self.df
        return CpGIntersects(df)
