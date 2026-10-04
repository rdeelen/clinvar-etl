from pathlib import Path
import os

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

cwd = os.getcwd()
cwd = cwd + "\\data\\variant_summary.txt.gz"
db = extract(cwd)
print(db)
print(len(db))