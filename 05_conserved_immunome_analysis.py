# ============================================================
# CD4 T CELL CONSERVED EPITOPE–HLA PAIR ANALYSIS PIPELINE
# ============================================================
# ANALYSES
# 1. Strict conserved pairs (without threshold filtering)
#    - Identify conserved epitope–HLA pairs across all variants
#    - Plot per-allele and per-locus conserved counts
# 2. Threshold-based conserved analysis
#    - Filter by NetMHCIIpan percentile thresholds
#    - Identify conserved epitope–HLA pairs across all variants
#    - Heatmap visualization
#    - Plot per-allele and per-locus conserved counts
# ============================================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math

variant_files = {
    "Ancestral": "ancestral.tsv",
    "BA.2.86": "ba286.tsv",
    "JN.1": "jn1.tsv",
    "LP.8.1": "lp81.tsv",
    "NB.1.8.1": "nb181.tsv",
    "XEC": "xec.tsv"
}
# ============================================================
# 1. STRICT CONSERVED PAIRS (WITHOUT THRESHOLD FILTERING)
# ============================================================

variant_pairs = {}
for variant, filepath in variant_files.items():
    df = pd.read_csv(filepath, sep="\t")
    df["peptide"] = df["peptide"].str.strip()
    df["allele"] = df["allele"].str.strip()
    variant_pairs[variant] = set(zip(df["peptide"], df["allele"]))

# Find conserved pairs across all variants
conserved_pairs = set.intersection(*variant_pairs.values())
print(f"Number of conserved epitope–HLA pairs: {len(conserved_pairs)}")

# Save conserved pairs (TSV + Excel)
conserved_df = pd.DataFrame(list(conserved_pairs), columns=["peptide", "allele"])
conserved_df.to_csv("conserved_epitope_HLA_pairs.tsv", sep="\t", index=False)
conserved_df.to_excel("conserved_epitope_HLA_pairs.xlsx", index=False)

# Per-allele counts
allele_counts = {}
for pep, allele in conserved_pairs:
    allele_counts[allele] = allele_counts.get(allele, 0) + 1

allele_df = pd.DataFrame(list(allele_counts.items()), columns=["allele", "count"])
allele_df = allele_df.sort_values(by="count", ascending=False)

# Save table
allele_df.to_csv("conserved_epitope_counts_per_allele.tsv", sep="\t", index=False)

plt.figure(figsize=(10, 6))
bars = plt.barh(
    allele_df["allele"],
    allele_df["count"],
    color="firebrick",
    edgecolor="dimgray"
)
plt.gca().invert_yaxis()
plt.xlabel("Number of Conserved Epitopes", fontsize=12)
plt.ylabel("HLA Allele", fontsize=12)
plt.title("Conserved Epitope–HLA Pairs per Allele (≤10 percentile)", fontsize=12, color="black")
plt.yticks(fontsize=10)  # allele list font size
plt.grid(axis="x", linestyle="--", alpha=0.7, color="dimgray", zorder=0)
for bar in bars:
    bar.set_zorder(3)
plt.tight_layout()
plt.show()

# Per-locus counts
locus_counts = {}
for pep, allele in conserved_pairs:
    locus = allele.split("-")[1][:2]  # DP, DQ, DR
    locus_counts[locus] = locus_counts.get(locus, 0) + 1

locus_df = pd.DataFrame(list(locus_counts.items()), columns=["locus", "count"])
locus_df = locus_df.sort_values(by="count", ascending=False)

# Save table
locus_df.to_csv("conserved_epitope_counts_per_locus.tsv", sep="\t", index=False)

plt.figure(figsize=(6, 5))
bars = plt.bar(
    locus_df["locus"],
    locus_df["count"],
    color="firebrick",
    edgecolor="dimgray"
)

# Add numbers inside bars
for bar in bars:
    height = bar.get_height()
    plt.text(
        bar.get_x() + bar.get_width()/2., height - (0.25 * height),
        f"{int(height)}", ha="center", va="bottom",
        color="white", fontsize=11
    )

plt.ylabel("Number of Conserved Epitopes", fontsize=12)
plt.xlabel("HLA Locus", fontsize=12)
plt.title("Per-locus Conserved Epitope Counts (≤10 percentile)", fontsize=12, color="black")
plt.grid(axis="y", linestyle="--", alpha=0.7, color="dimgray", zorder=0)
for bar in bars:
    bar.set_zorder(3)
plt.tight_layout()
plt.show()

# ============================================================
# 2. THRESHOLD-BASED CONSERVED ANALYSIS
# ============================================================

def analyze_conserved(threshold, label, chunk_size=50, show_heatmap=False):
    print(f"\n=== Analysis for threshold: ≤{threshold} ===")

    # Load filtered epitope–HLA sets
    variant_pairs = {}
    for variant, filepath in variant_files.items():
        df = pd.read_csv(filepath, sep="\t")
        df["peptide"] = df["peptide"].str.strip()
        df["allele"] = df["allele"].str.strip()
        df = df[df["netmhciipan_el percentile"] <= threshold]
        variant_pairs[variant] = set(zip(df["peptide"], df["allele"]))

    # Conserved across all variants
    conserved_pairs = set.intersection(*variant_pairs.values())
    print(f"Number of conserved epitope–HLA pairs (≤{threshold}): {len(conserved_pairs)}")

    # 1. Heatmap of conserved pairs
    if show_heatmap:
        heatmap_df = pd.DataFrame(index=[f"{pep}_{hla}" for pep, hla in conserved_pairs])
        for variant, pairs in variant_pairs.items():
            heatmap_df[variant] = heatmap_df.index.map(
                lambda x: tuple(x.split("_")) in pairs
            )
        heatmap_data = heatmap_df.astype(int)

        n_chunks = math.ceil(len(heatmap_data) / chunk_size)
        for i in range(n_chunks):
            chunk = heatmap_data.iloc[i*chunk_size:(i+1)*chunk_size]
            if chunk.empty:
                continue
            plt.figure(figsize=(10, 0.25*len(chunk)))
            sns.heatmap(
                chunk,
                cmap=["firebrick", "lightgray"],
                cbar=False,
                linewidths=0.3,
                linecolor="white"
            )
            plt.title(f"Conserved Epitope–HLA Pairs (≤{threshold} percentile) [Part {i+1}/{n_chunks}]")
            plt.xlabel("Variants")
            plt.ylabel("Peptide–HLA Pairs")
            plt.yticks(rotation=0, fontsize=6)
            plt.tight_layout()
            plt.show()

    # 2. Per-allele counts 
    allele_counts = {}
    for pep, allele in conserved_pairs:
        allele_counts[allele] = allele_counts.get(allele, 0) + 1

    allele_df = pd.DataFrame(list(allele_counts.items()), columns=["allele", "count"])
    allele_df = allele_df.sort_values(by="count", ascending=False)

    plt.figure(figsize=(10, 6))
    plt.barh(allele_df["allele"], allele_df["count"], color="firebrick", edgecolor="dimgray", zorder=2)
    plt.gca().invert_yaxis()
    plt.xlabel("Number of Conserved Epitopes", fontsize=12)
    plt.ylabel("HLA Allele", fontsize=12)
    plt.title(f"Conserved Epitope–HLA Pairs per Allele (≤{threshold} percentile)", fontsize=12)
    plt.yticks(fontsize=10)
    plt.grid(axis="x", linestyle="--", alpha=0.6, color="dimgray", zorder=0)
    plt.tight_layout()
    plt.show()

    # 3. Per-locus counts
    locus_counts = {}
    for pep, allele in conserved_pairs:
        locus = allele.split("-")[1][:2]  # DP, DQ, DR
        locus_counts[locus] = locus_counts.get(locus, 0) + 1

    locus_df = pd.DataFrame(list(locus_counts.items()), columns=["locus", "count"])
    locus_df = locus_df.sort_values(by="count", ascending=False)

    plt.figure(figsize=(6, 5))
    bars = plt.bar(locus_df["locus"], locus_df["count"], color="firebrick", edgecolor="dimgray", zorder=2)

    # Add numbers inside bars
    for bar in bars:
        height = bar.get_height()
        plt.text(
            bar.get_x() + bar.get_width()/2., height - (0.3 * height),
            f"{int(height)}", ha="center", va="bottom",
            color="white", fontsize=11
        )

    plt.ylabel("Number of Conserved Epitopes", fontsize=12)
    plt.xlabel("HLA Locus", fontsize=12)
    plt.title(f"Per-locus Conserved Epitope Counts (≤{threshold} percentile)", fontsize=12)
    plt.grid(axis="y", linestyle="--", alpha=0.6, color="dimgray", zorder=0)
    plt.tight_layout()
    plt.show()

# Run analyses
analyze_conserved(5, "≤5", show_heatmap=False)
analyze_conserved(2, "≤2", show_heatmap=False)
analyze_conserved(1, "≤1", show_heatmap=True)  
