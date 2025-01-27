import pandas as pd
import subprocess
import concurrent.futures

def read_table(path, usecols: list=None, test_run: bool=False):
    default_usecols = ["Chromosome", "Start", "End"]

    if test_run:
        nrows=1000000
    else: 
        nrows=None

    if usecols:
        if type(usecols) == list:
            default_usecols.extend(usecols)
        else: 
            default_usecols.append(usecols)

    df = pd.read_table(path, sep="\t", 
                       usecols=default_usecols,
                       nrows=nrows)
    return df

def fetch_modbeds(dirpath, usecols=None, test_run: bool = False):
    path_ls = subprocess.check_output(["ls", dirpath]).decode("utf-8").split("\n")
    path_ls.pop(-1)

    with concurrent.futures.ThreadPoolExecutor(len(path_ls)) as read_executor:
        tables = read_executor.map(lambda path: read_table(dirpath + path, usecols, test_run), path_ls)
        tables = [table.rename(columns={"N_mC" : "N_5mC", "N_hmC" : "N_5hmC"}) for table in tables]
    
    return tables

def calculate_percentages(df: pd.DataFrame, cols) -> pd.DataFrame:
    if type(cols) == list:
        for col in cols:
            df.eval(f"percentMeth_{col.split('_')[1]} = ({col}/readCount)*100", inplace=True)
    else:
        df.eval(f"percentMeth_{cols.split('_')[1]} = ({cols}/readCount)*100", inplace=True)

    return df

def merge_positions(dfs, calculate_percentage=False, cols=None, drop=None) -> pd.DataFrame:
    merged = (pd.concat(dfs)
              .groupby(["Chromosome", "Start", "End"], observed=False, as_index=False, sort=False)
              .sum(numeric_only=True))

    # Replace with two functions/methods that do this
    if calculate_percentage:
        merged = calculate_percentages(merged, cols)

    if drop:
        merged = merged.drop(columns=drop)

    return merged
