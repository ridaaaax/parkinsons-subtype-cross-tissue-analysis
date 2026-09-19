import pandas as pd

signature = pd.read_csv(
    "results/brain_signature_42_genes.csv"
)

muscle = pd.read_csv(
    "data/GSE128177_featureCounts.txt.gz",
    sep="\t",
    compression="gzip"
)

print("\n--- BRAIN SIGNATURE IDs ---")
print(signature["Gene"].head(10).to_string(index=False))

print("\n--- MUSCLE GENE IDs ---")
print(muscle["Geneid"].head(10).to_string(index=False))

print("\n--- DATA TYPES ---")
print("Brain Gene:", signature["Gene"].dtype)
print("Muscle Geneid:", muscle["Geneid"].dtype)

print("\n--- CLEANED IDs ---")

brain_ids = (
    signature["Gene"]
    .astype(str)
    .str.strip()
    .str.split(".")
    .str[0]
)

muscle_ids = (
    muscle["Geneid"]
    .astype(str)
    .str.strip()
    .str.split(".")
    .str[0]
)

print("Brain:")
print(brain_ids.head(10).to_string(index=False))

print("\nMuscle:")
print(muscle_ids.head(10).to_string(index=False))

print("\n--- OVERLAP TEST ---")

overlap = set(brain_ids) & set(muscle_ids)

print("Brain genes:", len(set(brain_ids)))
print("Muscle genes:", len(set(muscle_ids)))
print("Matching genes:", len(overlap))

print("\nMatching IDs:")
print(list(overlap))