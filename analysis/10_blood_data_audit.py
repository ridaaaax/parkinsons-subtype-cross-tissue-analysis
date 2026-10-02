import pandas as pd
from pathlib import Path

# -----------------------------------
# 1. Load blood expression matrix
# -----------------------------------

data_file = Path("data/GSE165082_PD-CC.counts.txt.gz")

print("Loading GSE165082 blood expression matrix...")

df = pd.read_csv(
    data_file,
    sep="\t",
    compression="gzip"
)

print("Data loaded successfully.")

# -----------------------------------
# 2. Dataset dimensions
# -----------------------------------

print("\n--- DATASET DIMENSIONS ---")
print("Rows:", df.shape[0])
print("Columns:", df.shape[1])
print("Number of genes:", df.shape[0])
print("Number of samples:", df.shape[1] - 1)

# -----------------------------------
# 3. Sample names
# -----------------------------------

print("\n--- SAMPLE NAMES ---")

sample_names = list(df.columns[1:])

for sample in sample_names:
    print(sample)

# -----------------------------------
# 4. Gene ID check
# -----------------------------------

print("\n--- GENE ID CHECK ---")

print("First 5 gene IDs:")
print(df.iloc[:, 0].head().to_string(index=False))

print(
    "Duplicated gene IDs:",
    df.iloc[:, 0].duplicated().sum()
)

# -----------------------------------
# 5. Missing values
# -----------------------------------

print("\n--- MISSING VALUES ---")

expression_data = df.iloc[:, 1:]

missing_values = expression_data.isna().sum().sum()

print("Total missing values:", missing_values)

# -----------------------------------
# 6. Zero values
# -----------------------------------

print("\n--- ZERO VALUES ---")

zero_count = (expression_data == 0).sum().sum()

print("Total zero entries:", zero_count)

# -----------------------------------
# 7. Integer count check
# -----------------------------------

print("\n--- VALUE TYPE CHECK ---")

fractional_count = (
    (expression_data % 1) != 0
).sum().sum()

print(
    "Fractional/non-integer entries:",
    fractional_count
)

# -----------------------------------
# 8. Negative values
# -----------------------------------

negative_count = (
    (expression_data < 0)
).sum().sum()

print("Negative entries:", negative_count)

# -----------------------------------
# 9. Sample totals
# -----------------------------------

print("\n--- SAMPLE TOTALS ---")

sample_totals = expression_data.sum()

for sample, total in sample_totals.items():
    print(f"{sample}: {total:.2f}")

# -----------------------------------
# 10. Save audit table
# -----------------------------------

audit = pd.DataFrame({
    "sample": sample_names,
    "library_total": [
        sample_totals[sample]
        for sample in sample_names
    ]
})

output_dir = Path("results")
output_dir.mkdir(exist_ok=True)

output_file = (
    output_dir /
    "GSE165082_data_audit.csv"
)

audit.to_csv(
    output_file,
    index=False
)

print("\n--- AUDIT COMPLETE ---")
print("Audit table saved to:")
print(output_file)