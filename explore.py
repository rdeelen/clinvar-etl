import pandas as pd




df = pd.read_csv("data/variant_summary.txt.gz", sep="\t", nrows=1000)
print(df.columns.tolist())
print(df.head())
print(df["ClinicalSignificance"].value_counts().head(15))