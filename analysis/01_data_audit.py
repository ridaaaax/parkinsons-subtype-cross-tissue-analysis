import pandas as pd
from pathlib import Path


# -----------------------------
# 1. Locate the data file
# -----------------------------

data_file = Path("data/GSE128177_featureCounts.txt.gz")


# -----------------------------
# 2. Load the expression matrix
# -----------------------------

print("Loading GSE128177 expression matrix...")

df = pd.read_csv(
    data_file,
    sep="\t",
    compression="gzip"
)

print("Data loaded successfully.")


# -----------------------------
# 3. Basic dimensions
# -----------------------------

print("\n--- DATASET DIMENSIONS ---")

print("Rows:", df.shape[0])
print("Columns:", df.shape[1])

print("Number of genes:", df.shape[0])
print("Number of samples:", df.shape[1] - 1)


# -----------------------------
# 4. Check column names
# -----------------------------

print("\n--- SAMPLE NAMES ---")

sample_names = list(df.columns[1:])

for sample in sample_names:
    print(sample)


# -----------------------------
# 5. Check gene IDs
# -----------------------------

print("\n--- GENE ID CHECK ---")

print("First 5 gene IDs:")
print(df["Geneid"].head().to_string(index=False))

print("Duplicated gene IDs:", df["Geneid"].duplicated().sum())


# -----------------------------
# 6. Check missing values
# -----------------------------

print("\n--- MISSING VALUES ---")

missing_values = df.isna().sum().sum()

print("Total missing values:", missing_values)


# -----------------------------
# 7. Check zero values
# -----------------------------

print("\n--- ZERO VALUES ---")

expression_data = df.iloc[:, 1:]

zero_count = (expression_data == 0).sum().sum()

print("Total zero entries:", zero_count)


# -----------------------------
# 8. Check fractional values
# -----------------------------

print("\n--- VALUE TYPE CHECK ---")

fractional_count = (
    (expression_data % 1) != 0
).sum().sum()

print("Fractional/non-integer entries:", fractional_count)


# -----------------------------
# 9. Library-size totals
# -----------------------------

print("\n--- SAMPLE TOTALS ---")

sample_totals = expression_data.sum()

for sample, total in sample_totals.items():
    print(f"{sample}: {total:.2f}")


# -----------------------------
# 10. Define sample groups
# -----------------------------

groups = {
    "M-M-01": "Parkinson_disease",
    "M-M-04": "Parkinson_disease",
    "M-M-08": "Parkinson_disease",
    "M-M-12": "Parkinson_disease",
    "M-M-16": "Parkinson_disease",
    "M-M-20": "Parkinson_disease",
    "M-M-24": "Parkinson_disease",
    "M-M-27": "Parkinson_disease",
    "M-M-30": "Parkinson_disease",
    "M-M-33": "Parkinson_disease",
    "M-M-36": "Parkinson_disease",
    "M-M-39": "Parkinson_disease",

    "M-M-02": "Old_adult",
    "M-M-06": "Old_adult",
    "M-M-10": "Old_adult",
    "M-M-14": "Old_adult",
    "M-M-18": "Old_adult",
    "M-M-22": "Old_adult",
    "M-M-25": "Old_adult",
    "M-M-28": "Old_adult",
    "M-M-31": "Old_adult",
    "M-M-34": "Old_adult",
    "M-M-37": "Old_adult",
    "M-M-40": "Old_adult",

    "M-M-03": "Young_adult",
    "M-M-07": "Young_adult",
    "M-M-11": "Young_adult",
    "M-M-15": "Young_adult",
    "M-M-19": "Young_adult",
    "M-M-23": "Young_adult",
    "M-M-26": "Young_adult",
    "M-M-29": "Young_adult",
    "M-M-32": "Young_adult",
    "M-M-35": "Young_adult",
    "M-M-38": "Young_adult",
    "M-M-41": "Young_adult"
}


# -----------------------------
# 11. Check group assignment
# -----------------------------

print("\n--- SAMPLE GROUPS ---")

for sample in sample_names:

    if sample in groups:
        print(f"{sample}: {groups[sample]}")
    else:
        print(f"{sample}: GROUP NOT FOUND")


# -----------------------------
# 12. Create an audit table
# -----------------------------

audit = pd.DataFrame({
    "sample": sample_names,
    "group": [groups.get(sample, "UNKNOWN") for sample in sample_names],
    "library_total": [sample_totals[sample] for sample in sample_names]
})


# -----------------------------
# 13. Save the audit table
# -----------------------------

output_file = Path("results")

output_file.mkdir(exist_ok=True)

audit.to_csv(
    output_file / "GSE128177_data_audit.csv",
    index=False
)

print("\n--- AUDIT COMPLETE ---")

print("Audit table saved to:")
print("results/GSE128177_data_audit.csv")