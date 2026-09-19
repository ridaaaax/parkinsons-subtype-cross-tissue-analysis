import pandas as pd
from pathlib import Path

# Location of the supplementary Table S3 file
data_file = Path("data/brain_signature/Table S3.xlsx")

print("Loading Table S3...")

df = pd.read_excel(data_file)

print("Table loaded successfully.")

print("\n--- TABLE DIMENSIONS ---")
print("Rows:", df.shape[0])
print("Columns:", df.shape[1])

print("\n--- COLUMN NAMES ---")
print(df.columns.tolist())

print("\n--- FIRST 5 ROWS ---")
print(df.head())

# Columns we need for our brain subtype signature
signature = df[
    [
        "Gene",
        "logFC",
        "P.Value",
        "adj.P.Val",
        "HGNC",
        "Entrez",
        "description"
    ]
].copy()

print("\n--- BRAIN SIGNATURE ---")
print("Number of genes:", len(signature))

print(signature.to_string(index=False))

# Create results folder if it does not exist
Path("results").mkdir(exist_ok=True)

# Save clean signature
output_file = "results/brain_signature_42_genes.csv"

signature.to_csv(
    output_file,
    index=False
)

print("\n--- COMPLETE ---")
print("Brain signature saved to:")
print(output_file)