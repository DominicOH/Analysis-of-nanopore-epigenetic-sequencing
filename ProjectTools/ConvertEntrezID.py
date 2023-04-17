import pandas as pd
import mygene

def convertENSID():
    ens_df = pd.read_csv("./feature_references/genes/GENCODE_Basic_mm39_promotersSubtracted.bed", sep="\t",
                     names=["Chromosome", "Start", "End", "ENSID", "Type"])
    ens_df["ENSID"] = ens_df["ENSID"].str.split(".", expand=True)[0]

    mg = mygene.MyGeneInfo()
    query_df = pd.DataFrame(mg.querymany(ens_df["ENSID"], scopes="ensembl.transcript", fields="symbol"))
    output_df = query_df.merge(ens_df.rename(columns={"ENSID" : "query"}), "inner", "query")[["Chromosome", "Start", "End", "symbol", "query"]]

    return output_df.to_csv("./feature_references/genes/GENCODE_Basic_mm39_entrezConverted.bed", sep="\t", header=None, index=False)

if __name__ == "__main__":
    convertENSID()