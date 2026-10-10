from pathlib import Path
from collections import Counter
import os
import time

import pandas as pd


GENES = {
    "BRCA1": "DNA repair (homologous recombination)",
    "BRCA2": "DNA repair (homologous recombination)",
    "PALB2": "BRCA2 partner, DNA repair",
    "ATM": "DNA damage response kinase",
    "CHEK2": "Cell-cycle checkpoint kinase",
    "TP53": "Tumour suppressor",
}

# ClinVar column name -> our column name
COLUMNS = {
    "#AlleleID": "allele_id",
    "VariationID": "variation_id",
    "GeneSymbol": "gene",
    "Type": "variant_type",
    "ClinicalSignificance": "significance_raw",
    "Assembly": "assembly",
    "Chromosome": "chromosome",
    "Start": "start_pos",
    "Stop": "stop_pos",
    "ReviewStatus": "review_status",
    "LastEvaluated": "last_evaluated",
    "PhenotypeList": "phenotypes",
}


def extract(path: Path, assembly="GRCh38", chunksize=200_000):
    """Reads in only the needed columns and rows in chunks"""
    header = pd.read_csv(path, sep="\t", nrows=0)
    missing = set(COLUMNS) - set(header.columns)
    if missing:
        raise ValueError(f"Input file is missing expected columns: {sorted(missing)}")
    
    parts = []
    for chunk in pd.read_csv(path, sep="\t", usecols=list(COLUMNS), chunksize=chunksize, low_memory=False):
        chunk = chunk[(chunk["Assembly"] == assembly) & (chunk["GeneSymbol"].isin(GENES))]
        chunk = chunk.rename(columns = COLUMNS)
        parts.append(chunk)
    return pd.concat(parts, ignore_index=True)


def classify_significance(raw):
    """streamlines the significance categories"""
    for index, row in raw.iterrows():
        if ("Conflicting" in row["significance_raw"]):
            raw.at[index, "significance_raw"] = "Conflicting"
        elif (row["significance_raw"].startswith("Pathogenic")):
            raw.at[index, "significance_raw"] = "Pathogenic"
        elif (row["significance_raw"].startswith("Likely pathogenic")):
            raw.at[index, "significance_raw"] = "Likely pathogenic"
        elif ("Uncertain" in row["significance_raw"]):
            raw.at[index, "significance_raw"] = "VUS"
        elif (row["significance_raw"].startswith("Benign")):
            raw.at[index, "significance_raw"] = "Benign"
        elif (row["significance_raw"].startswith("Likely benign")):
            raw.at[index, "significance_raw"] = "Likely benign"
        else:
            raw.at[index, "significance_raw"] = "Other"
    raw.rename(columns={"significance_raw": "significance"}, inplace=True)
    return raw


def review_stars(raw):
    """turns the value of review status into the star number on the webpage"""
    for index, row in raw.iterrows():
        if ("practice guideline" in row["review_status"]):
            raw.at[index, "review_status"] = "4"
        elif ("reviewed by expert panel" in row["review_status"]):
            raw.at[index, "review_status"] = "3"
        elif ("multiple submitters, no conflicts" in row["review_status"]):
            raw.at[index, "review_status"] = "2"
        elif ("criteria provided" in row["review_status"] and "no assertion" not in row["review_status"]):
            raw.at[index, "review_status"] = "1"
        else:
            raw.at[index, "review_status"] = "0"
    return raw


def transform(raw):
    """transforms the dataframe to filter for false entries, duplicates
    and streamlines review status and significance classes"""
    raw = classify_significance(raw)
    raw = review_stars(raw)
    raw["last_evaluated"] = pd.to_datetime(raw["last_evaluated"], format="%b %d, %Y", errors="coerce")
    raw["start_pos"] = pd.to_numeric(raw["start_pos"], errors="coerce")
    raw["stop_pos"] = pd.to_numeric(raw["stop_pos"], errors="coerce")
    #raw = raw[raw["stop_pos"] < raw["start_pos"]]
    raw.dropna(axis= 0, subset= ["stop_pos", "start_pos"], inplace= True)
    raw.drop_duplicates(subset= ["variation_id", "assembly"], inplace= True)
    return raw


def validate(data):
    """Checks the dataframe for significant errors (wrong genes, wrong review status)"""
    if (data.empty):
        raise ValueError("Dataframe is empty")
    
    header = data.columns.values
    header = Counter(header)
    if (header["variation_id"] == 2):
        raise ValueError("variation_id exists twice")
    if (header["assembly"] == 2):
        raise ValueError("assembly exists twice")
    
    review_value = set(data["review_status"])
    review_value.discard("0")
    review_value.discard("1")
    review_value.discard("2")
    review_value.discard("3")
    review_value.discard("4")
    if (bool(review_value)):
        raise ValueError("review_status is out of range")
    
    genes_df = set(data["gene"])
    genes_list = list(GENES.keys())
    for value in genes_df:
        if (value not in genes_list):
            raise ValueError(f"wrong gene in gene list: {value}")
    for value in genes_list:
        if (value not in genes_df):
            raise ValueError(f"gene is missing: {value}")
    
    

cwd = os.getcwd()
cwd = cwd + "\\data\\variant_summary.txt.gz"
start = time.time()
db = extract(cwd)
db = transform(db)
end = time.time()
total_time = (end - start)/60
print(total_time)
validate(db)
#print(db["review_status"])
#print(db["significance"].value_counts())
#print(db["start_pos"].dtype)