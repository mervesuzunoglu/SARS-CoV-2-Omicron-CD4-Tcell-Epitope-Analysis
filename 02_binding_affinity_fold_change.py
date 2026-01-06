# Fold change analysis, Binding Strength Change Classification, Preserved vs. Unique Peptide–Allele Pairs

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# --- SETTINGS ---
PERCENTILE_COLUMN = 'netmhciipan_el percentile'
PEPTIDE_COL = 'peptide'
ALLELE_COL = 'allele'
STRONG_BINDER_THRESHOLD = 10

# --- File names ---
variant_files = {
    "BA.2.86": "ba286.tsv",
    "JN.1": "jn1.tsv",
    "LP.8.1": "lp81.tsv",
    "NB.1.8.1": "nb181.tsv",
    "XEC": "xec.tsv"
}
ancestral_file = "ancestral.tsv"

# --- Custom colors ---
custom_palette = {
    "Increased": "firebrick",   # red
    "Decreased": "royalblue",
    "Unchanged": "dimgray"      # gray
}
preserved_palette = ["dimgray", "firebrick"]

# --- Load ancestral data ---
ancestral_df = pd.read_csv(ancestral_file, sep="\t")
ancestral_df = ancestral_df[ancestral_df[PERCENTILE_COLUMN] <= STRONG_BINDER_THRESHOLD]
ancestral_df["pair"] = ancestral_df[PEPTIDE_COL] + "_" + ancestral_df[ALLELE_COL]
ancestral_pairs = ancestral_df.set_index("pair")[PERCENTILE_COLUMN]

# --- Initialize summary list ---
summary_data = []
variant_pair_sets = []

sns.set(style="whitegrid")

# --- Fold change analysis ---
for variant, file in variant_files.items():
    variant_df = pd.read_csv(file, sep="\t")
    variant_df = variant_df[variant_df[PERCENTILE_COLUMN] <= STRONG_BINDER_THRESHOLD]
    variant_df["pair"] = variant_df[PEPTIDE_COL] + "_" + variant_df[ALLELE_COL]
    variant_pairs = variant_df.set_index("pair")[PERCENTILE_COLUMN]
    variant_pair_sets.append(set(variant_pairs.index))

    # Match with ancestral
    shared = ancestral_pairs.index.intersection(variant_pairs.index)
    fold_changes = variant_pairs[shared] / ancestral_pairs[shared]

    # Classify fold change
    def classify(fc):
        if np.isclose(fc, 1.0, atol=0.1):
            return "Unchanged"
        elif fc < 1:
            return "Decreased"
        else:
            return "Increased"

    categories = fold_changes.apply(classify)

    # Count categories
    counts = categories.value_counts().to_dict()
    summary_data.append({
        "Variant": variant,
        "Shared_Pairs": len(shared),
        "Increased": counts.get("Increased", 0),
        "Decreased": counts.get("Decreased", 0),
        "Unchanged": counts.get("Unchanged", 0)
    })

    # Plot histogram of fold changes
    plt.figure(figsize=(6, 4))
    sns.histplot(np.log2(fold_changes), bins=30, kde=True, color="darkred")
    plt.axvline(0, color="black", linestyle="--")
    plt.title(f"Log2 Fold Change: {variant} vs. Ancestral")
    plt.xlabel("Log2(Fold Change of Binding Percentile)")
    plt.ylabel("Peptide–Allele Pair Count")
    plt.tight_layout()
    plt.show()

# --- Summary Table ---
summary_df = pd.DataFrame(summary_data)

# --- Plot summary bar chart ---
variant_order = ["BA.2.86", "JN.1", "LP.8.1", "NB.1.8.1", "XEC"]

summary_df_melted = summary_df.melt(
    id_vars=["Variant"],
    value_vars=["Increased", "Decreased", "Unchanged"],
    var_name="Change", value_name="Count"
)
plt.figure(figsize=(8, 5))
sns.barplot(
    data=summary_df_melted,
    x="Variant", y="Count", hue="Change",
    order=variant_order,
    palette=custom_palette
)
plt.title("Binding Strength Change Classification")
plt.ylabel("Peptide–Allele Pair Count")
plt.tight_layout()
plt.show()

# --- Preserved epitope–allele pairs across all variants + ancestral ---
common_preserved_pairs = set(ancestral_pairs.index).intersection(*variant_pair_sets)
print(f"\nNumber of preserved peptide–allele pairs across all: {len(common_preserved_pairs)}")

# --- Plot preserved vs. unique counts ---
preserved_counts = {
    "Preserved in All": len(common_preserved_pairs),
    "Unique": len(set.union(*variant_pair_sets, set(ancestral_pairs.index))) - len(common_preserved_pairs)
}
plt.figure(figsize=(5, 4))
sns.barplot(x=list(preserved_counts.keys()), y=list(preserved_counts.values()), palette=preserved_palette)
plt.title("Preserved vs. Unique Peptide–Allele Pairs")
plt.ylabel("Count")
plt.xticks(rotation=0)
plt.tight_layout()
plt.show()

# --- Save summary table ---
summary_df.to_csv("summary_table.tsv", sep="\t", index=False)
