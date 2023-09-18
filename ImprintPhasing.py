import pandas as pd
from scipy.spatial import distance
from scipy.cluster import hierarchy
import functools
import seaborn as sns

class ModBaseAssemblage:
    """
    A pandas DataFrame object built to contain modified bases 
    """
    def __init__(self, df: pd.DataFrame):
        self.df = df
        self._score_table = None
    
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
    
    @property
    def score_table(self):
        if self._score_table is None:
            self._score_table = self.aggregated_score_table()
        return self._score_table

    def read_matrix(self, gene_name=None, minimum_read_count_proportion=0.5, min_cpg_count_proportion=0.5):
        df = self.aggregated_score_table().replace(["c", "m", "h"], [0, 1, 2])

        if gene_name: 
            matrix = df.query(f"Name == '{gene_name}' & refPos > 0")
        else:  # currently doesn't work
            matrix = df.query("refPos > 0")

        matrix = matrix.pivot(index="readID", columns=["refPos"], values="classification")
        return DMR_Matrix(matrix, self.score_table, minimum_read_count_proportion, min_cpg_count_proportion)

class DMR_Matrix:
    def __init__(self, df, score_table=None, minimum_read_count_proportion=0.5, min_cpg_count_proportion=0.5):
        self.df = self.drop_reads(df, minimum_read_count_proportion, min_cpg_count_proportion)
        self.__original_score_table = score_table
        self._clusters = None
        self._linkage_matrix = None

    def drop_reads(self, df, minimum_read_count_proportion=0.5, min_cpg_count_proportion=0.5):
        filtered_df = df.copy()

        # first: remove CpGs present in fewer than the minimum count of reads
        filtered_df = filtered_df.dropna(thresh=minimum_read_count_proportion*len(filtered_df.index), axis="columns") 
    
        # next: remove reads containing fewer than the minimum count of CpG sites
        filtered_df = filtered_df.dropna(thresh=min_cpg_count_proportion*len(filtered_df.columns), axis="index") 

        return filtered_df
    
    @property
    def linkage_matrix(self):
        if self._linkage_matrix is None:
            p_distances = distance.pdist(self.df, "hamming")
            linkage = hierarchy.linkage(p_distances, "average")
            self._linkage_matrix = linkage
        return self._linkage_matrix

    @property
    def clusters(self):
        if self._clusters is None:
            linkage = self.linkage_matrix
            clusters = hierarchy.fcluster(linkage, 2, "maxclust")
            self._clusters = clusters
        return self._clusters

    def demo_clustermap(self):
        df = self.df
        cm = sns.clustermap(df.fillna(-1), 
                            mask=df.isnull(), 
                            row_linkage=self.linkage_matrix, col_cluster=False, 
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
    
    def as_modDF(self, merge_sites=False, keep_borders=False):
        score_table = self.__original_score_table.query("refPos >= 1")
        clusters = self.demo_phased_reads_as_df()

        unmerged_read_table = score_table.merge(clusters, "inner", "readID").rename(columns={
            "refPos" : "Start"
        }).assign(End = lambda row: row["Start"] + 1) 

        if not merge_sites:
            unmerged_read_output = unmerged_read_table.loc[:, ("Chromosome", "Start", "End", "regionStart", "regionEnd", "readID", "classification", "Cluster")]

            if not keep_borders:
                unmerged_read_output = unmerged_read_output.query("Start >= regionStart & End <= regionEnd")
            
            c1_meth = unmerged_read_output.groupby("Cluster")["classification"].value_counts()[(1, "m")]
            c2_meth = unmerged_read_output.groupby("Cluster")["classification"].value_counts()[(2, "m")]

            if c1_meth > c2_meth:
                output_df = unmerged_read_output.replace({"Cluster" : {1 : "Methylated", 2 : "Unmethylated"}})
            else: 
                output_df = unmerged_read_output.replace({"Cluster" : {2 : "Methylated", 1 : "Unmethylated"}})
       
            return output_df
        else: 
            grouped_read_table = unmerged_read_table.groupby([
                "Start", "End", "Chromosome", "Cluster", "classification"
                ]).size().reset_index().rename(
                        columns={0 : "size"}).pivot(
                                index=["Chromosome", "Start", "End", "Cluster"], 
                                columns="classification", values="size").fillna(0).reset_index(col_level="Chromosome")

            grouped_read_table = grouped_read_table.assign(
                readCount = lambda row: row["c"] + row["h"] + row["m"])
            
            grouped_read_table = grouped_read_table.assign(
                percentMeth_C = lambda row: row["c"] / row["readCount"],
                percentMeth_5mC = lambda row: row["m"] / row["readCount"],
                percentMeth_5hmC = lambda row: row["h"] / row["readCount"]
            )
            grouped_read_table = grouped_read_table.loc[:, ("Chromosome", "Start", "End", "c", "m", "h", "percentMeth_C", "percentMeth_5mC", "percentMeth_5hmC", "readCount", "Cluster")]
            return grouped_read_table