import pandas as pd
from scipy.spatial import distance
from scipy.cluster import hierarchy
import seaborn as sns

class ModBaseAssemblage:
    """
    A pandas DataFrame object built to contain modified bases 
    """
    def __init__(self, df: pd.DataFrame):
        self.df = df
    
    def __pivot_mods(self):
        df = self.df
        index = ["Name", "Chromosome", "regionStart", "regionEnd", "readStart", "readEnd", "readID", "refPos", "qPos", "strand"]
        pivoted_df = df.pivot(index=index, columns="modBase", values="modScore").reset_index()

        return pivoted_df
    
    def aggregated_score_table(self): 
        pivoted_df = self.__pivot_mods()

        score_table = pivoted_df.assign(c = lambda r: 255 - r["m"] - r["h"])
        score_table = score_table.assign(classification = score_table.loc[:, ("h", "m", "c")].idxmax(axis=1))
        score_table["refPos"] = score_table.apply(lambda row: row["refPos"] if row["strand"] == "+" else row["refPos"] - 1, axis=1)

        return score_table

    def read_matrix(self, min_reads_per_cpg=None, min_cpgs_per_read=None, gene_name=None):
        df = self.aggregated_score_table().replace(["c", "m", "h"], [0, 1, 2])

        if gene_name: 
            matrix = df.query(f"Name == '{gene_name}' & refPos > 0")
        else: 
            matrix = df.query("refPos > 0")

        matrix = matrix.pivot(index="readID", columns=["refPos"], values="classification")
        return DMR_Matrix(matrix, min_reads_per_cpg, min_cpgs_per_read)
    
class DMR_Matrix:
    def __init__(self, df, min_reads_per_cpg=None, min_cpgs_per_read=None):
        self.df = self.dropReads(df, min_reads_per_cpg, min_cpgs_per_read)
        self._clusters = None
        self._linkage_matrix = None

    def dropReads(self, df, min_reads_per_cpg=None, min_cpgs_per_read=None):
        filtered_df = df.copy()

        # first: remove CpGs present in fewer than the minimum count of reads
        if not min_reads_per_cpg:
            filtered_df = filtered_df.dropna(thresh=0.5*len(filtered_df.index), axis="columns")  
        else: 
            filtered_df = filtered_df.dropna(thresh=min_reads_per_cpg, axis="columns") 
        
        # next: remove reads containing fewer than the minimum count of CpG sites.
        if not min_cpgs_per_read:
            filtered_df = filtered_df.dropna(thresh=0.5*len(filtered_df.columns), axis="index")
        else: 
            filtered_df = filtered_df.dropna(thresh=min_cpgs_per_read, axis="index") 

        return filtered_df
    
    @property
    def linkage_matrix(self):
        if not self._linkage_matrix:
            p_distances = distance.pdist(self.df, "hamming")
            linkage = hierarchy.linkage(p_distances, "average")
            self._linkage_matrix = linkage
        return self._linkage_matrix

    @property
    def clusters(self):
        if not self._clusters:
            linkage = self.linkage_matrix
            clusters = hierarchy.fcluster(linkage, 2, "maxclust")
            self._clusters = clusters
        return self._clusters

    def demo_clustermap(self):
        cm = sns.clustermap(self.df.fillna(-1), 
                            mask=self.df.isnull(), 
                            row_linkage=self.linkage_matrix,
                            col_cluster=False, cmap=sns.color_palette("Blues", 3), 
                            method="average", metric="hamming")
        return cm

    def demo_phased_reads_as_df(self):
        clusters = self.clusters
        read_ids = self.df.T.columns
        phased_reads_df = pd.DataFrame(zip(list(read_ids), clusters),
                                    columns=["readID", "Cluster"])
                                    
        return phased_reads_df

    def output_readIDs(self, outname):
        df = self.demo_phased_reads_as_df()
        
        for cluster in [1, 2]:
            cluster_df = df.query(f"Cluster == {cluster}")
            cluster_df.to_csv(f"./data_tables/DMR_analysis/{outname}_c{cluster}.txt", header=False, index=False, columns=["readID"])
        return 
    
    def as_modDF(self, 
                 dmr_df: ModBaseAssemblage,
                 merge_sites=False
                 ):
        pivoted_dmr_df = dmr_df.aggregated_score_table().query("refPos >= 1")
        clusters = self.demo_phased_reads_as_df()

        merged_with_clustered_reads = pivoted_dmr_df.merge(clusters, "inner", "readID").rename(columns={
            "refPos" : "Start"
        }).assign(End = lambda row: row["Start"] + 1)        
        
        if not merge_sites:
            return merged_with_clustered_reads.loc[:, ("Chromosome", "Start", "End", "readID", "classification", "Cluster")]
        
        count_of_modbases = merged_with_clustered_reads.groupby([
            "Start", "End", "Chromosome", "Cluster", "classification"
            ]).size().reset_index().rename(
                    columns={0 : "size"}).pivot(
                            index=["Chromosome", "Start", "End", "Cluster"], 
                            columns="classification", values="size").fillna(0).reset_index(col_level="Chromosome")

        count_of_modbases = count_of_modbases.assign(
            readCount = lambda row: row["c"] + row["h"] + row["m"])
        
        count_of_modbases = count_of_modbases.assign(
            percentMeth_C = lambda row: row["c"] / row["readCount"],
            percentMeth_5mC = lambda row: row["m"] / row["readCount"],
            percentMeth_5hmC = lambda row: row["h"] / row["readCount"]
        )

        count_of_modbases = count_of_modbases.loc[:, ("Chromosome", "Start", "End", "c", "m", "h", "percentMeth_C", "percentMeth_5mC", "percentMeth_5hmC", "readCount", "Cluster")]
        return count_of_modbases