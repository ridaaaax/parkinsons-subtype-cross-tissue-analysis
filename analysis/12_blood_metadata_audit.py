import gzip
import os
import pandas as pd


INPUT_FILE = "data/GSE165082_family.soft.gz"
OUTPUT_FILE = "results/GSE165082_blood_metadata_audit.csv"

os.makedirs("results", exist_ok=True)


print("Reading GEO family SOFT file...")


# --------------------------------------------------
# Read sample-level metadata from SOFT
# --------------------------------------------------

sample_data = {}
current_sample = None

with gzip.open(INPUT_FILE, "rt", encoding="utf-8", errors="replace") as f:

    for line in f:

        line = line.rstrip("\n")

        # Start of a new sample block
        if line.startswith("!Sample_geo_accession"):
            parts = line.split("\t")

            if len(parts) >= 2:
                current_sample = parts[1].strip()
                sample_data[current_sample] = {}

        elif current_sample is not None and line.startswith("!Sample_"):

            parts = line.split("\t")

            field = parts[0].replace("!Sample_", "")

            values = [x.strip() for x in parts[1:]]

            sample_data[current_sample][field] = " | ".join(values)


# --------------------------------------------------
# Convert to dataframe
# --------------------------------------------------

metadata = pd.DataFrame.from_dict(
    sample_data,
    orient="index"
)

metadata.index.name = "GEO_accession"


print("\nNumber of GEO samples found:", len(metadata))


# --------------------------------------------------
# Print useful fields
# --------------------------------------------------

useful_fields = [
    "geo_accession",
    "title",
    "characteristics_ch1",
    "source_name_ch1",
    "organism_ch1",
    "library_source",
    "library_strategy",
    "instrument_model",
    "description"
]

print("\nAvailable metadata fields:")

for column in metadata.columns:
    print("-", column)


# --------------------------------------------------
# Search metadata for age / sex / RIN
# --------------------------------------------------

search_terms = ["age", "sex", "rin"]

print("\nSearching metadata for age / sex / RIN...\n")

matching_columns = []

for column in metadata.columns:

    column_text = metadata[column].fillna("").astype(str)

    combined_text = " ".join(column_text.tolist()).lower()

    matched = [
        term for term in search_terms
        if term in column.lower() or term in combined_text
    ]

    if matched:
        matching_columns.append(column)
        print(f"Possible match: {column} -> {matched}")


# --------------------------------------------------
# Save complete metadata
# --------------------------------------------------

metadata.to_csv(OUTPUT_FILE)

print("\nMetadata audit saved to:")
print(OUTPUT_FILE)


# --------------------------------------------------
# Print relevant metadata values
# --------------------------------------------------

if matching_columns:

    print("\nRelevant metadata:\n")

    for column in matching_columns:

        print("=" * 70)
        print("FIELD:", column)
        print("=" * 70)

        for accession, value in metadata[column].items():
            print(accession, ":", value)

else:

    print("\nNo age/sex/RIN information was detected in the SOFT metadata.")