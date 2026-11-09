import os
import gzip
import pandas as pd

INPUT_FILE = "data/GSE165082_family.soft.gz"
OUTPUT_FILE = "results/GSE165082_blood_metadata_audit.csv"

os.makedirs("results", exist_ok=True)

sample_data = {}
current_sample = None

with gzip.open(INPUT_FILE, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        line = line.strip()

        if line.startswith("^SAMPLE = "):
            current_sample = line.split("=", 1)[1].strip()
            sample_data[current_sample] = {}

        elif current_sample and line.startswith("!Sample_"):
            if " = " not in line:
                continue

            field, value = line.split(" = ", 1)
            field = field.replace("!Sample_", "")

            sample_data[current_sample].setdefault(field, [])
            sample_data[current_sample][field].append(value.strip())

metadata = pd.DataFrame.from_dict(sample_data, orient="index")
metadata.index.name = "GEO_accession"

for column in metadata.columns:
    metadata[column] = metadata[column].apply(
        lambda x: " | ".join(x) if isinstance(x, list) else x
    )

metadata.to_csv(OUTPUT_FILE)

print("Samples extracted:", len(metadata))
print("\nAvailable metadata fields:")
print(metadata.columns.tolist())

print("\nAge / sex / RIN-related fields:")

found = False

for column in metadata.columns:
    values = metadata[column].fillna("").astype(str)

    if any(term in column.lower() for term in ["age", "sex", "rin"]):
        found = True
        print("\nFIELD:", column)
        print(values.to_string())

    elif values.str.contains(
        r"\bage\b|\bsex\b|\bRIN\b",
        case=False,
        regex=True
    ).any():
        found = True
        print("\nFIELD:", column)
        print(values.to_string())

if not found:
    print("No age, sex or RIN fields detected.")

print("\nSaved:", OUTPUT_FILE)