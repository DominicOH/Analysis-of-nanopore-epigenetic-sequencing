import pandas as pd
import pyranges as pr
from scipy.spatial import distance
from scipy.cluster import hierarchy
import seaborn as sns
import matplotlib.pyplot as plt

class ModkitExtract():
    """
    Class for working with the tsv output of Modkit Extract. 
    """
    def __init__(self, path: str):
        self.path = path
        self._cpg_table = None 
        self.read_table = self.__read_table()

    @property
    def cpg_table(self):
        if self._cpg_table == None:
           cpg_table = self.__cpg_table()
           self._cpg_table = cpg_table
           return cpg_table
        else: 
            return self._cpg_table
        
    def __read_table(self):
        read_table = (pd.read_table(self.path, sep="\t")
                      .pivot(index=["chrom", "read_id", "ref_position", "ref_strand", "read_length"],
                             columns="mod_code", values="mod_qual")
                      .reset_index()
                      .assign(c = lambda r: round(1 - (r["h"] + r["m"]), 6))  # classify modified bases - modcalls discretised as c, m, and h 
                      .assign(classification = lambda r: r.loc[:, ("h", "m", "c")].idxmax(axis=1))
                      .drop(columns=["c", "m", "h"]))  

        return read_table
        
    def __cpg_table(self):
        cpg_t = self.read_table.query("ref_position != -1") # remove unmapped bases

        # change to pyranges default formatting
        cpg_t["Start"] = cpg_t["ref_position"]
        cpg_t.eval("End = Start + 1", inplace=True)

        # remove unnecessary columns
        cpg_t = cpg_t.rename(columns={"chrom" : "Chromosome", "ref_strand" : "Strand"})
        cpg_t = cpg_t.loc[:, ("Chromosome", "Start", "End", "Strand", "read_id", "classification")]
        return CpGTable(cpg_t)
        
class CpGTable():
    def __init__(self, df, include_bed=None):
        self.df = df
        self._include_bed = include_bed

    def set_include_bed(self, include):
        if type(include) == str:
            self._include_bed = pr.read_bed(include)
        elif type(include) == pr.PyRanges:
            self._include_bed = include
        elif type(include) == pd.DataFrame:
            try: 
                include_pr = pr.PyRanges(include)
                self._include_bed = include_pr
            except:
                raise ValueError(f"Could not coerce {include} to PyRanges.")
    
    @property
    def include_bed(self):
        return self._include_bed.as_df()
        
    def select_gene(self, gene_name):
        assert self._include_bed is not None, "Need to attach a target region bedfile with '.set_include_bed()"

        df = self.df
        
        include_bed = self._include_bed
        read_pr = pr.PyRanges(df)
        annotated_df = read_pr.join(include_bed, False, suffix="_Gene").as_df()

        annotated_df = annotated_df.query(f"Name == '{gene_name}'")
        new_read_table = annotated_df.loc[:, ("Chromosome", "Start", "End", "Strand", "read_id", "classification")]
        return GeneCpGTable(new_read_table, gene_name)
    
class GeneCpGTable(CpGTable):
    def __init__(self, df, gene_name, include_bed=None):
        super().__init__(df, include_bed)
        self.gene_name = gene_name
        self._read_matrix = self.__generate_read_matrix()
        self._linkage_matrix = self.__generate_linkage_matrix()
        self._cluster_table = self.__generate_flat_clusters()
    
    def __generate_read_matrix(self, 
                               minimum_read_proportion: float = 0.1, 
                               min_cpg_proportion: float = 0.15):
        df = self.df 
        read_matrix = df.pivot(index="read_id", columns="Start", values="classification")
        read_matrix = read_matrix.replace(["c", "m", "h"], [0, 1, 2])

        
        # first: remove CpGs present in fewer than the minimum count of reads
        total_reads = len(read_matrix.index)
        read_matrix = read_matrix.dropna(thresh=minimum_read_proportion*total_reads, 
                      axis="columns") 
        
        # next: remove reads containing fewer than the minimum count of CpG sites
        total_sites = len(read_matrix.columns)
        read_matrix = read_matrix.dropna(thresh=min_cpg_proportion*total_sites, 
                      axis="index") 
        
        return read_matrix
    
    def __generate_linkage_matrix(self):
        read_matrix = self._read_matrix
        distance_matrix = distance.pdist(read_matrix, "hamming")
        linkage_matrix = hierarchy.linkage(distance_matrix, "average")
        return linkage_matrix
    
    def __generate_flat_clusters(self):
        # Note that this is optimised for DMRs with two detected haplotypes
        linkage_matrix = self._linkage_matrix
        cluster_table = hierarchy.fcluster(linkage_matrix, 2, "maxclust")
        return cluster_table
    
    def __update_matrices(self, minimum_read_proportion, min_cpg_proportion):
        self._read_matrix = self.__generate_read_matrix(minimum_read_proportion, min_cpg_proportion)
        self._linkage_matrix = self.__generate_linkage_matrix()
        self._cluster_table = self.__generate_flat_clusters()
        return 
    
    def clustermap(self, 
                   minimum_read_proportion: float = 0.1, 
                   min_cpg_proportion: float = 0.15,
                   quiet: bool = False):
        """
        Plot clusters of reads within the read table using a seaborn clustermap. Image is figure level. For axes level plotting use 'heatmap()'.

        minimum_read_proportion : float 
            CpG sites must be present in at least this proportion of all reads in the dataframe (default: 0.1)
        min_cpg_proportion : float 
            Reads must cover at least this proportion of all CpG sites in the dataframe (default: 0.15).
        """
        if minimum_read_proportion != 0.1 or min_cpg_proportion != 0.15:
            self.__update_matrices(minimum_read_proportion, min_cpg_proportion)

        cm = sns.clustermap(self._read_matrix.fillna(-1), 
                            mask=self._read_matrix.isnull(), 
                            row_linkage=self._linkage_matrix, 
                            col_cluster=False, 
                            method="average", metric="hamming",
                            xticklabels="auto", yticklabels=False,
                            cmap=sns.color_palette("Blues", 3),
                            cbar_kws={"location" : "top",
                                      "orientation" : "horizontal",
                                      "ticks" : [0, 1, 2]})
        
        cm.ax_cbar.set_position([0.40, 0.85, 0.3, 0.03])
        cm.ax_cbar.set_xticks([0.33, 1, 1.66])
        cm.ax_cbar.set_xticklabels(["C", "5mC", "5hmC"])
        cm.ax_cbar.set_title("Modification type", fontdict={"fontsize" : 10})
        cm.ax_heatmap.set_ylabel("Read")
        cm.ax_heatmap.set_xlabel("CpG position")

        if quiet:
            plt.close()
            return cm
        return cm
    
    def heatmap(self, 
                minimum_read_proportion: float = 0.1, 
                min_cpg_proportion: float = 0.15, 
                fontsize: int = 5,
                ax=None):
        
        """
        Plot clusters of reads within the read table using a seaborn heatmap. The axes-level equivalent of 'clustermap()'.

        minimum_read_proportion : float 
            CpG sites must be present in at least this proportion of all reads in the dataframe (default: 0.1)
        min_cpg_proportion : float 
            Reads must cover at least this proportion of all CpG sites in the dataframe (default: 0.15).
        """
        data2d = self.clustermap(minimum_read_proportion, 
                                 min_cpg_proportion,
                                 quiet=True).data2d
        mask = data2d == -1

        if not ax:
            fig, ax = plt.subplots()
        
        hm = sns.heatmap(data2d, mask=mask,
                xticklabels=False, yticklabels=False,
                cmap=sns.color_palette("YlGnBu"),
                cbar=False,
                ax=ax)
    
        right = ax.get_xlim()[1]*0.9
        left = ax.get_xlim()[0] + 0.14*ax.get_xlim()[1]

        cluster_props = self.cluster_extract().groupby(["cluster", "read_id"]).count().reset_index()["cluster"].value_counts()

        chromosome = self.df["Chromosome"].values[0]

        ax.set_xticks([left, right], 
            labels=[chromosome + ": " + str(data2d.columns[0]), str(data2d.columns[-1])], rotation="horizontal", fontsize=fontsize)
        ax.tick_params("x", bottom=False)
    
        ax.set_title(f"{self.gene_name}\nU = {cluster_props['Unmethylated']}; M = {cluster_props['Methylated']}", loc="center", fontsize=fontsize)
        ax.set_ylabel(None)
        ax.set_xlabel(None)

        return hm
    
    def cluster_extract(self, 
                        minimum_read_proportion: float = None,
                        min_cpg_proportion: float = None):
        """
        Clusters reads within the table based on modification state.
        """
        
        if minimum_read_proportion or min_cpg_proportion:
            self.__update_matrices(minimum_read_proportion, min_cpg_proportion)

        read_matched_with_cluster = pd.DataFrame(zip(list(self._read_matrix.T.columns), self._cluster_table),
            columns=["read_id", "cluster"]) 
        
        annotated_table = pd.merge(self.df, read_matched_with_cluster, "inner")

        c2_meth = annotated_table.groupby("cluster")["classification"].value_counts()[(2, "m")]
        c1_meth = annotated_table.groupby("cluster")["classification"].value_counts()[(1, "m")]

        if c1_meth > c2_meth:
            output_df = annotated_table.replace({"cluster" : {1 : "Methylated", 2 : "Unmethylated"}})
        else: 
            output_df = annotated_table.replace({"cluster" : {2 : "Methylated", 1 : "Unmethylated"}})

        return output_df

    def get_unmethylated_readIDs(self,
                                 minimum_read_proportion = None,
                                 min_cpg_proportion = None
                                 ):
        
        """
        Extracts read IDs from the unmethylated allele. 
        """
        gene_table = self.cluster_extract(minimum_read_proportion, min_cpg_proportion)
        cluster_groups = gene_table.groupby("cluster")
        read_array = cluster_groups.get_group("Unmethylated")["read_id"].unique()

        read_ids = pd.Series([id for id in read_array])

        return read_ids
    
    def get_methylated_readIDs(self,
                               minimum_read_proportion = None,
                               min_cpg_proportion = None
                               ):
        
        """
        Extracts read IDs from the methylated allele. 
        """
        gene_table = self.cluster_extract(minimum_read_proportion, min_cpg_proportion)
        cluster_groups = gene_table.groupby("cluster")
        read_array = cluster_groups.get_group("Methylated")["read_id"].unique()

        read_ids = pd.Series([id for id in read_array])

        return read_ids