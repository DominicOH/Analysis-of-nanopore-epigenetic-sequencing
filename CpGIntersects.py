import subprocess
import pandas as pd
import pyranges as pr
from FeatureReferences import *

def geneRefPyRange():
    gene_ref_path = './feature_references/revised/GENCODE_Basic_mm39_Genes_merged.bed'
    df = pd.read_csv(gene_ref_path, sep="\t", names=["Chromosome", "Start", "End", "Name", "Strand"])
    return pr.PyRanges(df).unstrand()

def geneFeatureRefPyRange():
    gene_feature_list = subprocess.check_output(["ls", "./feature_references/revised/gene_features/name_adjusted/"]).decode("utf-8").split("\n") 
    gene_feature_list.pop(-1)

    df_list = []
    for file in gene_feature_list:
        path = "./feature_references/revised/gene_features/name_adjusted/" + file
        feature_tsv = Features(path)
        feature_df = feature_tsv.toDF()
        df_list.append(feature_df)

    feature_reference_df = pd.concat(df_list).drop(columns=["Score"])
    return pr.PyRanges(feature_reference_df)

def repeatTypeRefPyRange():
    repeat_ref_path = "./feature_references/revised/repeats/UCSC_rmsk_mm39_Repeat_IDd.bed"
    df = pd.read_csv(repeat_ref_path, sep="\t", names=["Chromosome", "Start", "End", "feature_type", "Name"])
    return pr.PyRanges(df).unstrand()

def CGIrefPyRange():
    cgi_feature_list = subprocess.check_output(["ls", "./feature_references/revised/cgi/named/"]).decode("utf-8").split("\n") 
    cgi_feature_list.pop(-1)

    cgi_df_list = []
    for file in cgi_feature_list:
        path = "./feature_references/revised/cgi/named/" + file
        cgi_tsv = CGIs(path)
        cgi_df = cgi_tsv.toDF()
        cgi_df_list.append(cgi_df)

    cgi_reference_df = pd.concat(cgi_df_list)
    return pr.PyRanges(cgi_reference_df)

class CpGIntersects(pr.PyRanges):
    """
    Main class for feature/gene level comparison. Inherits from PyRange. 
    """
    def __init__(self, df, target_mod="5hmC"):
        super().__init__(df, df)
        self.target_mod = target_mod

    def intersectGenes(self):
        """
        Intersects CpGs with genes. Based on gene start/end coordinates in the GENCODE Basic reference build. 
        """
        gene_ref_pr = geneRefPyRange()
        df_with_genes = self.join(gene_ref_pr, slack=0).as_df()
        df_with_genes["feature_type"] = "Gene"
        return  df_with_genes
    
    def intersectFeatures(self):
        """
        Intersects CpGs with genomic features. Output is a dataframe-type object. 
        """
        feature_ref = geneFeatureRefPyRange().unstrand()
        df_with_features = self.join(feature_ref, slack=0).as_df()
        categories = ["Intergenic", "Repeat", "Promoter", "5UTR", "TSS", "Intron", "Exon", "3UTR", "TTS"]
        df_with_features["feature_type"] = pd.Categorical(df_with_features["feature_type"], categories)

        return  df_with_features
    
    def intersectRepeatTypes(self):
        """
        Intersects CpGs with repetitive features in DNA. Output is a dataframe-type object. 
        """
        repeat_ref = repeatTypeRefPyRange().unstrand()
        df_labelled_repeats = self.join(repeat_ref, slack=0).as_df()
        categories = ["LINE", "SINE", "Simple_repeat", "LTR", "DNA", "Retroposon", "Low_complexity", "Satellite"]
        df_labelled_repeats["feature_type"] = pd.Categorical(df_labelled_repeats["feature_type"], categories)

        return  df_labelled_repeats
    
    def intersectCpGIslands(self):
        """
        Intersects CpGs with CpG islands. Islands are broken into island feature (i.e.: shelf, shore). 
        """
        cgi_ref = CGIrefPyRange().unstrand()
        df_with_cgis = self.join(cgi_ref, slack=0).as_df()
        categories = ["Open sea", "Upstream shelf", "Upstream shore", "CGI", "Downstream shore", "Downstream shelf"]
        df_with_cgis["feature_type"] = pd.Categorical(df_with_cgis["feature_type"], categories)

        return  df_with_cgis
    
    def groupByGenomicWindow(self, window_size):
        """
        Groups CpGs based according to Xkb windows ("tiles"), using the average (mean) hydroxymethlyation of CpGs within those windows. Output is distinct from the grouping function below as the chromosomal coordinates are actually what defines each cluster. 
        """
        tiled_pr = self.tile(window_size, strand=False).cluster(slack=-1, strand=False)
        cluster_pr = tiled_pr.merge(slack=-1, count=True, strand=False)

        # requires columns: ["percentMeth_5hmC_Min" & "percentMeth_5hmC_Bisulphite"] | ["percentMeth_5mC_Min" & "percentMeth_5mC_Bisulphite"]
        
        if self.target_mod == "5hmC":
            cluster_pr = cluster_pr.insert(
                tiled_pr.apply(f=lambda df: df.groupby("Cluster")[["percentMeth_5hmC_Bisulphite", "percentMeth_5hmC_Min"]].mean(), as_pyranges=False, strand=False))
        elif self.target_mod == "5mC":
            cluster_pr = cluster_pr.insert(
                tiled_pr.apply(f=lambda df: df.groupby("Cluster")[["percentMeth_5mC_Bisulphite", "percentMeth_5mC_Min"]].mean(), as_pyranges=False, strand=False))

        if "percentMeth_5hmC_Bisulphite" and "percentMeth_5hmC_Min" in tiled_pr.as_df().columns:
            cluster_pr = cluster_pr.insert(
                tiled_pr.apply(
                f=lambda df: df.groupby("Cluster")[["percentMeth_5hmC_Bisulphite", "percentMeth_5hmC_Min"]].mean(), 
                as_pyranges=False, strand=False)
                )
        elif "percentMeth_5hmC_Bisulphite" and "percentMeth_5hmC_Prom" in tiled_pr.as_df().columns:
            cluster_pr = cluster_pr.insert(
                tiled_pr.apply(
                f=lambda df: df.groupby("Cluster")[["percentMeth_5hmC_Bisulphite", "percentMeth_5hmC_Prom"]].mean(), 
                as_pyranges=False, strand=False)
                ) 
            # the following header names are out of date as of this current version
        elif "percentMeth_Bisulphite_5hmC" and "percentMeth_Nanopore_5hmC" in tiled_pr.as_df().columns:
            cluster_pr = cluster_pr.insert(
                tiled_pr.apply(
                f=lambda df: df.groupby("Cluster")[["percentMeth_Bisulphite_5hmC", "percentMeth_Nanopore_5hmC"]].mean(), 
                as_pyranges=False, strand=False)
                )  
        elif "percentMeth_TAB" and "percentMeth_Nanopore" in tiled_pr.as_df().columns:
            cluster_pr = cluster_pr.insert(
                tiled_pr.apply(
                f=lambda df: df.groupby("Cluster")[["percentMeth_TAB", "percentMeth_Nanopore"]].mean(), 
                as_pyranges=False, strand=False)
                )
        else: 
            raise ValueError(f"Check column name consistency: {tiled_pr.as_df().columns}")
        
        return cluster_pr.as_df()
    
    def group(self, intersect_with):
        """
        Groups CpGs based on intersects.
        """
        if intersect_with == "other":
            intersect_df = self.df
        elif intersect_with == "features" or intersect_with == "Features":
            intersect_df = self.intersectFeatures()
        elif intersect_with == "CGI" or intersect_with == "islands":
            intersect_df = self.intersectCpGIslands()
        elif intersect_with == "genes" or intersect_with == "Genes":
            intersect_df = self.intersectGenes()
        elif intersect_with == "repeats" or intersect_with == "Repeats":
            intersect_df = self.intersectRepeatTypes()
        else: 
            raise ValueError("Please input appropriate element to intersect with: ['genes', 'features', 'repeats', or 'CGI']")
        
        groupby_df = intersect_df.groupby(["Name", "feature_type", "Start_b", "End_b"], 
                                          observed=True).agg(
            {"percentMeth_Nanopore_5hmC" : "mean",
             "percentMeth_TAB_5hmC" : "mean",
             "Start" : "count"}
             ).reset_index()
        
        groupby_df.rename(columns={"Start" : "CpG_count",
                                   "Start_b" : "group_start",
                                   "End_b" : "group_end",
                                   "percentMeth_TAB_5hmC" : "percentMeth_TAB",
                                   "percentMeth_Nanopore_5hmC" : "percentMeth_Nanopore"}, 
                                   inplace=True)

        return  groupby_df