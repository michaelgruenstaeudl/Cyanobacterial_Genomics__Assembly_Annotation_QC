import pandas as pd
import matplotlib.pyplot as plt
import argparse
import sys
import os

# ---------------------------
# Parse command-line arguments
# ---------------------------
parser = argparse.ArgumentParser(
    description="Generate a donut chart from a QIIME2 taxonomy table with fallback to lowest assigned rank."
)

parser.add_argument(
    "--manifest",
    required=True,
    help="Path to manifest.tsv (must contain 'sample-id' column)"
)

parser.add_argument(
    "--table",
    required=True,
    help="Path to taxonomy table (e.g., genus-table.tsv)"
)

parser.add_argument(
    "--threshold",
    required=True,
    type=float,
    help="Relative abundance threshold as a decimal fraction (e.g., 0.03 for 3%)"
)

args = parser.parse_args()

manifest_path = args.manifest
table_path = args.table
threshold = args.threshold

# ---------------------------
# Validate input files
# ---------------------------
if not os.path.exists(manifest_path):
    sys.exit(f"Error: manifest file not found: {manifest_path}")

if not os.path.exists(table_path):
    sys.exit(f"Error: taxonomy table not found: {table_path}")

if not 0 <= threshold <= 1:
    sys.exit("Error: --threshold must be between 0 and 1 (e.g., 0.03 for 3%).")

# ---------------------------
# Read sample name from manifest.tsv
# ---------------------------
manifest = pd.read_csv(manifest_path, sep="\t")
manifest.columns = manifest.columns.str.strip()

if "sample-id" not in manifest.columns:
    sys.exit(
        f"Error: 'sample-id' column not found in {manifest_path}. "
        f"Available columns: {list(manifest.columns)}"
    )

sample_name = str(manifest.loc[0, "sample-id"]).strip()

# ---------------------------
# Read taxonomy table
# ---------------------------
df = pd.read_csv(table_path, sep="\t", skiprows=1)
df.columns = df.columns.str.strip()

if "#OTU ID" not in df.columns:
    sys.exit(
        f"Error: '#OTU ID' column not found in {table_path}. "
        f"Available columns: {list(df.columns)}"
    )

df = df.rename(columns={"#OTU ID": "Taxon"})

if sample_name not in df.columns:
    sys.exit(
        f"Error: Sample '{sample_name}' not found in {table_path}. "
        f"Available columns: {list(df.columns)}"
    )

df[sample_name] = pd.to_numeric(df[sample_name], errors="coerce").fillna(0)

# ---------------------------
# Taxonomy parsing logic
# ---------------------------
RANK_NAMES = {
    "k": "Kingdom",
    "p": "Phylum",
    "c": "Class",
    "o": "Order",
    "f": "Family",
    "g": "Genus",
    "s": "Species",
}

def extract_lowest_assigned_taxon(taxon):
    taxon = str(taxon).strip()

    if taxon.lower().startswith("unassigned"):
        return "Unassigned"

    parts = [p.strip() for p in taxon.split(";")]

    assigned = []
    genus_value = None

    for part in parts:
        if "__" not in part:
            continue

        prefix, value = part.split("__", 1)
        prefix = prefix.strip()
        value = value.strip()

        if not value or value == "_":
            continue

        assigned.append((prefix, value))

        if prefix == "g":
            genus_value = value

    if genus_value:
        return genus_value

    if assigned:
        lowest_prefix, lowest_value = assigned[-1]
        rank_name = RANK_NAMES.get(lowest_prefix, lowest_prefix)
        return f"{lowest_value} ({rank_name})"

    return "Unassigned"

df["AssignedTaxon"] = df["Taxon"].apply(extract_lowest_assigned_taxon)

# ---------------------------
# Aggregate counts
# ---------------------------
taxon_counts = df.groupby("AssignedTaxon")[sample_name].sum()
taxon_counts = taxon_counts[taxon_counts > 0]

# Remove unassigned taxa
taxon_counts = taxon_counts[taxon_counts.index != "Unassigned"]

if taxon_counts.empty:
    sys.exit("Error: No positive counts found after removing unassigned taxa.")

rel = taxon_counts / taxon_counts.sum()

# ---------------------------
# Group minor taxa
# ---------------------------
major = rel[rel >= threshold].copy()
minor = rel[rel < threshold].sum()

if minor > 0:
    major.loc["Other"] = minor

major = major.sort_values(ascending=False)

# ---------------------------
# Plot donut chart
# ---------------------------
fig, ax = plt.subplots(figsize=(8, 8))

wedges, texts, autotexts = ax.pie(
    major,
    labels=None,
    autopct=lambda p: f"{p:.1f}%" if p >= threshold * 100 else "",
    startangle=90,
    wedgeprops={"width": 0.4, "edgecolor": "white", "linewidth": 1}
)

centre_circle = plt.Circle((0, 0), 0.60, fc="white")
ax.add_artist(centre_circle)

ax.text(0, 0, sample_name, ha="center", va="center", fontsize=12)

ax.legend(
    wedges,
    major.index,
    title="Assigned taxon",
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    frameon=False
)

ax.set_title(f"Taxonomic composition of {sample_name}")
plt.tight_layout()

# Save figures
plt.savefig(f"{sample_name}_donut.png", dpi=300, bbox_inches="tight")
plt.savefig(f"{sample_name}_donut.svg", bbox_inches="tight")

# ---------------------------
# Print summary table
# ---------------------------
summary = pd.DataFrame({
    "Count": taxon_counts,
    "RelativeAbundance": rel
}).sort_values("RelativeAbundance", ascending=False)

print(summary)
