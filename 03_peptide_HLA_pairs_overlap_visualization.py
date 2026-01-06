# ============================================================
# SARS-CoV-2 Epitope Analysis + Visualization Pipeline
# ============================================================
# Combined analyses:
# 1. Shared peptide–HLA pairs across variants (UpSet)
# 2. Per-variant histograms (epitope distribution per allele)
# 3. Grouped comparison across variants
# 4. Preserved peptide–HLA pairs across all variants
# 5. Mutation overlap of epitope–HLA pairs per variant
# 6. Affected epitope–HLA pair lists
# 7. Pie charts of affected vs unaffected epitopes
# 8. Heatmap of mutation overlap across variants
# 9. UpSet plot of conserved / unique / lost / shared sets
# ============================================================

%pip install upsetplot

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from upsetplot import UpSet, from_memberships

# ============================================================
# COMMON SETTINGS
# ============================================================
files = {
    "Ancestral": "ancestral.tsv",
    "BA.2.86": "ba286.tsv",
    "JN.1": "jn1.tsv",
    "NB.1.8.1": "nb181.tsv",
    "LP.8.1": "lp81.tsv",
    "XEC": "xec.tsv",
}

colors = {
    "Ancestral": "#552D50",
    "BA.2.86": "#CB555D",
    "JN.1": "#F3AC67",
    "NB.1.8.1": "#99BDB7",
    "LP.8.1": "#1C646B",
    "XEC": "#A53852",
}

STRONG_BINDER_THRESHOLD = 10
PERCENTILE_COL = "netmhciipan_el percentile"
PEPTIDE_COL = "peptide"
ALLELE_COL = "allele"
START_COL = "start"
END_COL = "end"

mutations_df = pd.read_excel("variant_mutations.xlsx")

# ============================================================
# 1. LOAD DATA FOR VARIANT COMPARISON + UPSET ANALYSIS
# ============================================================
variant_sets = {}   # For UpSet plots
allele_dfs = []     # For per-variant and grouped comparisons

for variant, filepath in files.items():
    if not os.path.exists(filepath):
        print(f"File {filepath} not found, skipping.")
        continue

    df = pd.read_csv(filepath, sep="\t")

    # --- For UpSet plot ---
    variant_sets[variant] = set(df[PEPTIDE_COL] + "_" + df[ALLELE_COL])

    # --- For histograms and comparison ---
    allele_counts = df.groupby(ALLELE_COL)[PEPTIDE_COL].nunique().reset_index()
    allele_counts["Variant"] = variant
    allele_dfs.append(allele_counts)

# ============================================================
# 2. UPSET PLOT – Shared peptide–allele pairs across variants
# ============================================================
all_pairs = set().union(*variant_sets.values())
memberships = [[v for v, s in variant_sets.items() if pair in s] for pair in all_pairs]
data = from_memberships(memberships)

plt.figure(figsize=(18, 6))
UpSet(data, subset_size="count", show_counts=True, sort_categories_by="-input").plot()
plt.suptitle("Shared Peptide–HLA Allele Pairs Across Variants", fontsize=12)
plt.tight_layout()
plt.show()

# ============================================================
# 3. PER-VARIANT HISTOGRAMS (Epitope counts per allele)
# ============================================================
for allele_counts in allele_dfs:
    variant = allele_counts["Variant"].iloc[0]
    counts = allele_counts.groupby(ALLELE_COL)[PEPTIDE_COL].sum().sort_values(ascending=False)

    plt.figure(figsize=(12, 6))
    ax = counts.plot(kind="bar", color="firebrick")
    plt.title(f"Distribution of Strong-Binding Epitopes per HLA Class II Allele\n{variant}", fontsize=14)
    plt.xlabel("HLA Class II Allele", fontsize=12)
    plt.ylabel("Number of Unique Epitopes", fontsize=12)
    plt.xticks(rotation=90)
    ax.set_axisbelow(True)
    plt.grid(True, color="lightgray", linestyle="-", linewidth=0.7)
    plt.tight_layout()
    plt.show()

# ============================================================
# 4. GROUPED BAR CHART (Comparison across variants)
# ============================================================
all_data = pd.concat(allele_dfs, ignore_index=True)
pivot_data = all_data.pivot(index=ALLELE_COL, columns="Variant", values=PEPTIDE_COL).fillna(0)

# Sort alleles by total epitope count
pivot_data["Total"] = pivot_data.sum(axis=1)
pivot_data = pivot_data.sort_values("Total", ascending=False).drop(columns="Total")

ax = pivot_data.plot(
    kind="bar",
    figsize=(14, 7),
    width=0.8,
    color=[colors[col] for col in pivot_data.columns]
)
plt.title("Comparison of Strong-Binding Epitope Counts Across Variants per HLA-II Allele", fontsize=14)
plt.xlabel("HLA Class II Allele", fontsize=12)
plt.ylabel("Number of Unique Epitopes", fontsize=12)
plt.xticks(rotation=90)
plt.legend(title="Variant", bbox_to_anchor=(1.05, 1), loc="upper left")
ax.set_axisbelow(True)
ax.grid(True, color="lightgray", linestyle="-", linewidth=0.7)
plt.tight_layout()
plt.savefig("epitope_variant_comparison.png", dpi=300, bbox_inches="tight")
plt.show()

# ============================================================
# 5. PRESERVED PEPTIDE–HLA PAIRS ACROSS ALL VARIANTS
# ============================================================
pair_sets = {}
for name, file in files.items():
    df = pd.read_csv(file, sep="\t")
    df = df[df[PERCENTILE_COL] <= STRONG_BINDER_THRESHOLD]
    df["pair"] = df[PEPTIDE_COL] + "_" + df[ALLELE_COL]
    pair_sets[name] = set(df["pair"])

preserved_all = set.intersection(*pair_sets.values())
print(f"Number of peptide–HLA pairs preserved in ancestral + all variants: {len(preserved_all)}")

preserved_df = pd.DataFrame(list(preserved_all), columns=["Peptide_HLA_Pair"])
preserved_df[['Peptide', 'HLA']] = preserved_df['Peptide_HLA_Pair'].str.split('_', n=1, expand=True)
preserved_path = "preserved_pairs_all_variants.xlsx"
preserved_df.to_excel(preserved_path, index=False)
print(f"Saved preserved pairs as {preserved_path}")

# ============================================================
# 6. MUTATION OVERLAP + AFFECTED PAIRS PER VARIANT
# ============================================================
def check_overlap(epitope_start, epitope_end, mutation_list):
    for mut in mutation_list:
        match = re.search(r"(\d+)", str(mut))
        if match:
            pos = int(match.group(1))
            if epitope_start <= pos <= epitope_end:
                return True
    return False

variant_output_files = {}

for variant, file in files.items():
    df = pd.read_csv(file, sep="\t")
    df = df[df[PERCENTILE_COL] <= STRONG_BINDER_THRESHOLD].copy()

    if variant != "Ancestral":
        mutation_list = mutations_df[variant].dropna().tolist()
        df["MutationOverlap"] = df.apply(
            lambda row: check_overlap(row[START_COL], row[END_COL], mutation_list), axis=1
        )
    else:
        df["MutationOverlap"] = False

    variant_full_path = f"{variant}_mutation_overlap.xlsx"
    df.to_excel(variant_full_path, index=False)
    variant_output_files[variant] = variant_full_path
    print(f"Saved {variant_full_path} with mutation overlaps")

    affected_df = df[df["MutationOverlap"] == True].copy()
    if not affected_df.empty:
        affected_path = f"{variant}_affected_epitope_HLA_pairs.xlsx"
        affected_df.to_excel(affected_path, index=False)
        print(f"Saved {affected_path} ({len(affected_df)} affected pairs)")
    else:
        print(f"No affected epitopes found for {variant}")

# ============================================================
# 7. PIE CHARTS OF AFFECTED VS UNAFFECTED EPITOPES
# ============================================================
def make_autopct(values):
    def my_autopct(pct):
        return '{:.1f}%'.format(pct)
    return my_autopct

for variant, excel_file in variant_output_files.items():
    df = pd.read_excel(excel_file)
    counts = df["MutationOverlap"].value_counts()
    affected = counts.get(True, 0)
    unaffected = counts.get(False, 0)

    plt.figure(figsize=(5, 5))
    wedges, texts, autotexts = plt.pie(
        [affected, unaffected],
        labels=[f"Affected ({affected})", f"Unaffected ({unaffected})"],
        colors=["firebrick", "dimgray"],
        autopct=make_autopct([affected, unaffected]),
        startangle=90
    )
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('regular')

    plt.title(f"{variant}: Strong-Binding Epitopes Affected by Mutations")
    plt.tight_layout()
    plt.show()

# ============================================================
# 8. HEATMAP OF MUTATION OVERLAP ACROSS VARIANTS
# ============================================================
heatmap_df = pd.DataFrame()

for variant, file in variant_output_files.items():
    if variant == "Ancestral":
        continue
    df = pd.read_excel(file)
    df["Peptide_HLA"] = df[PEPTIDE_COL] + "_" + df[ALLELE_COL]
    df["MutationOverlap"] = df["MutationOverlap"].fillna(False)
    agg = df.groupby("Peptide_HLA")["MutationOverlap"].any()
    heatmap_df[variant] = agg

heatmap_df = heatmap_df.fillna(False)
heatmap_data = heatmap_df.astype(int)
heatmap_data = heatmap_data[heatmap_data.sum(axis=1) > 0]

plt.figure(figsize=(12, max(6, 0.2*len(heatmap_data))))
sns.heatmap(
    heatmap_data,
    cmap=["lightgray", "firebrick"],
    cbar=False,
    linewidths=0.5,
    linecolor='white'
)
plt.title("Variant-Defining Mutation Overlap with Strong-Binding Epitope–HLA Pairs")
plt.xlabel("Variants")
plt.ylabel("Peptide–HLA Pairs")
plt.yticks(rotation=0)
plt.tight_layout()
plt.show()

# ============================================================
# 9. UPSET PLOT – Conserved / Unique / Lost / Shared sets
# ============================================================
variant_sets = {}
for name, filepath in files.items():
    df = pd.read_csv(filepath, sep="\t")
    variant_sets[name] = set(df[PEPTIDE_COL] + "_" + df[ALLELE_COL])

all_pairs = set().union(*variant_sets.values())

memberships = []
for pair in all_pairs:
    present_in = [variant for variant, s in variant_sets.items() if pair in s]
    memberships.append(present_in)

data = from_memberships(memberships)

categories = []
for present_in in memberships:
    if set(files.keys()) == set(present_in):
        categories.append("Conserved")
    elif len(present_in) == 1 and present_in[0] != "Ancestral":
        categories.append("Unique")
    elif "Ancestral" in present_in and len(present_in) < len(files):
        categories.append("Lost")
    else:
        categories.append("Shared")

data = pd.DataFrame(data)
data["Category"] = categories

plt.figure(figsize=(14, 7))
up = UpSet(
    data,
    subset_size="count",
    show_counts=True,
    sort_categories_by="-input",
)
up.plot()
plt.suptitle("UpSet Plot: Unique vs Conserved vs Lost Epitope Sets", fontsize=14)
plt.tight_layout()
plt.show()
