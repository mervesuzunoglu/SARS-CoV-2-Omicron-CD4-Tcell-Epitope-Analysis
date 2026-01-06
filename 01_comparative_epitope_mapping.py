# ============================================
# SARS-CoV-2 Epitope Analysis Pipeline
# Analyses included:
# 1. Epitope overlap (variants only / with ancestral)
# 2. Epitope preservation vs ancestral
# 3. Shared mutations across variants
# 4. Impact of variant-defining mutations on ancestral epitopes 
# 5. Variant-specific epitopes
# 6. Positional mapping of mutated/conserved epitopes
# ============================================

!pip install upsetplot

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from upsetplot import from_contents, UpSet
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import re

# ============================================
# HELPER FUNCTIONS
# ============================================

def plot_heatmap(matrix_dict, title, cmap, cbar_label="Shared Epitopes"):
    variants = list(matrix_dict.keys())
    n = len(variants)
    shared_matrix = pd.DataFrame(0, index=variants, columns=variants)
    for i in range(n):
        for j in range(n):
            shared_matrix.iloc[i, j] = len(matrix_dict[variants[i]].intersection(matrix_dict[variants[j]]))
    plt.figure(figsize=(8, 6))
    sns.heatmap(shared_matrix, annot=True, fmt="d", cmap=cmap, cbar_kws={'label': cbar_label})
    plt.yticks(rotation=0, va="center")
    plt.xticks(rotation=0, ha="right")
    plt.title(title, fontsize=12)
    plt.tight_layout()
    plt.show()

def parse_mutation(mutation):
    mutation = str(mutation).strip()
    positions = []
    if re.match(r"[A-Z]\d+[A-Z]", mutation):
        pos = int(re.findall(r"\d+", mutation)[0])
        positions.append(pos)
    elif "del" in mutation or "Δ" in mutation:
        nums = re.findall(r"\d+", mutation)
        if len(nums) == 1:
            positions.append(int(nums[0]))
        elif len(nums) == 2:
            start, end = map(int, nums)
            positions.extend(range(start, end + 1))
    elif "ins" in mutation:
        nums = re.findall(r"\d+", mutation)
        if nums:
            positions.append(int(nums[0]))
    return positions

def mutation_hits_epitope(mutation, epitopes_df):
    match = re.search(r"(\d+)", str(mutation))
    if not match:
        return False
    pos = int(match.group(1))
    for _, row in epitopes_df.iterrows():
        if row["Mapped Start Position"] <= pos <= row["Mapped End Position"]:
            return True
    return False

def get_mutated_epitopes(epitopes, mutations):
    mutated, unmutated = [], []
    for _, row in epitopes.iterrows():
        start, end, seq = row["Mapped Start Position"], row["Mapped End Position"], row["Epitope Sequence"]
        is_mutated = any(start <= int(re.search(r"(\d+)", mut).group(1)) <= end for mut in mutations if re.search(r"(\d+)", mut))
        if is_mutated:
            mutated.append(seq)
        else:
            unmutated.append(seq)
    return set(mutated), set(unmutated)

# ============================================
# COLORMAP
# ============================================
heatmap_cmap = mcolors.LinearSegmentedColormap.from_list("firebrick_gray", ["#FFF3F2", "#9C0000"])

# ============================================
# 1. EPITOPE OVERLAP (UPSET + HEATMAP)
# ============================================
def analyze_epitope_overlap(epitope_file, include_ancestral=False):
    df = pd.read_excel(epitope_file)
    if include_ancestral:
        print("\n=== Overlap including Ancestral ===")
    else:
        print("\n=== Overlap across Variants Only ===")

    variants = df.columns.tolist()
    epitope_dict = {v: set(df[v].dropna().astype(str).tolist()) for v in variants}

    # UpSet plot
    contents = from_contents(epitope_dict)
    plt.figure(figsize=(12, 6))
    UpSet(contents, subset_size='count', show_counts=True).plot()
    title = "Overlap of Experimental Epitopes"
    if include_ancestral:
        title += " (Including Ancestral)"
    else:
        title += " (Variants Only)"
    plt.suptitle(title, fontsize=12)
    plt.show()

    # Heatmap
    plot_heatmap(epitope_dict, f"Shared Epitopes Across Variants{' + Ancestral' if include_ancestral else ''}", heatmap_cmap)

# ============================================
# 2. EPITOPE PRESERVATION VS ANCESTRAL
# ============================================
def analyze_epitope_preservation(epitope_file, show_unmatched=False):
    df = pd.read_excel(epitope_file)
    df.columns = ["Ancestral", "BA.2.86", "JN.1", "LP.8.1", "NB.1.8.1", "XEC"]

    epitope_sets = {col: set(df[col].dropna().astype(str).str.strip()) for col in df.columns}
    ancestral = epitope_sets["Ancestral"]

    print("\n=== Epitope Preservation vs Ancestral ===")
    for variant in ["BA.2.86", "JN.1", "LP.8.1", "NB.1.8.1", "XEC"]:
        variant_set = epitope_sets[variant]
        overlap = ancestral.intersection(variant_set)
        unmatched = variant_set - ancestral

        overlap_count = len(overlap)
        total_variant = len(variant_set)
        percentage = (overlap_count / total_variant * 100) if total_variant > 0 else 0
        print(f"{overlap_count} of {total_variant} {variant} epitopes found in ancestral ({percentage:.2f}%)")

        if show_unmatched:
            if unmatched:
                print(f"Unmatched epitopes in {variant} ({len(unmatched)}):")
                for ep in unmatched:
                    print(f"  - {ep}")
            else:
                print(f"All {variant} epitopes are conserved in ancestral.")
        print("-" * 60)

# ============================================
# 3. SHARED MUTATIONS ACROSS VARIANTS
# ============================================
def analyze_shared_mutations(mutation_file):
    df = pd.read_excel(mutation_file)
    mutation_dict = {col: set(df[col].dropna().astype(str).str.strip().tolist()) for col in df.columns}
    plot_heatmap(mutation_dict, "Pairwise Shared Spike Mutations Across Variants", heatmap_cmap, cbar_label="Shared Mutations")

# ============================================
# 4. IMPACT OF VARIANT-DEFINING MUTATIONS ON ANCESTRAL EPITOPES
# ============================================
def analyze_variant_defining_mutations(epitope_file, mutation_file, output_excel="epitope_mutation_status.xlsx"):
    # Load data
    epitopes_df = pd.read_excel(epitope_file)
    mutations_df = pd.read_excel(mutation_file)
    epitopes_df.columns = ["Epitope Sequence", "Mapped Start Position", "Mapped End Position"]
    variants = list(mutations_df.columns)

    results = []
    mutated_epitopes_by_variant = {}
    unmutated_epitopes_by_variant = {}

    # Collect mutated/unmutated epitopes
    for variant in variants:
        variant_mutations = mutations_df[variant].dropna().tolist()
        mut_positions = set()
        for mut in variant_mutations:
            mut_positions.update(parse_mutation(mut))

        mutated, unmutated = [], []
        for _, row in epitopes_df.iterrows():
            start, end = int(row["Mapped Start Position"]), int(row["Mapped End Position"])
            seq = row["Epitope Sequence"]
            ep_range = set(range(start, end + 1))
            if ep_range & mut_positions:
                mutated.append(seq)
            else:
                unmutated.append(seq)

        total = len(mutated) + len(unmutated)
        results.append({
            "Variant": variant,
            "Mutated %": len(mutated) / total * 100 if total else 0,
            "Unmutated %": len(unmutated) / total * 100 if total else 0,
            "Total Epitopes": total,
            "Mutated Count": len(mutated),
            "Unmutated Count": len(unmutated)
        })

        mutated_epitopes_by_variant[variant] = mutated
        unmutated_epitopes_by_variant[variant] = unmutated

    results_df = pd.DataFrame(results)
    print("\n=== Variant-Defining Mutation Impact ===")
    print(results_df)

    # Export mutated/unmutated epitope lists
    with pd.ExcelWriter(output_excel) as writer:
        for variant in variants:
            df_out = pd.DataFrame({
                "Mutated Epitopes": pd.Series(mutated_epitopes_by_variant[variant]),
                "Unmutated Epitopes": pd.Series(unmutated_epitopes_by_variant[variant])
            })
            df_out.to_excel(writer, sheet_name=variant, index=False)
    print(f"\nMutated/unmutated epitope lists exported to '{output_excel}'")

    # ============================================
    # 4a. STACKED BAR PLOT OF MUTATED/UNMUTATED PERCENTAGES
    # ============================================

    fig, ax = plt.subplots(figsize=(9,6))
    plot_data = results_df.set_index("Variant")[["Unmutated %", "Mutated %"]]
    plot_data.plot(kind="bar", stacked=True, ax=ax, color=["#949494", "#B22222"], edgecolor="white")

    for i, (unmut, mut) in enumerate(zip(results_df["Unmutated Count"], results_df["Mutated Count"])):
        ax.text(i, results_df.loc[i, "Unmutated %"]/2, str(unmut), ha="center", va="center", color="white", fontsize=10)
        ax.text(i, results_df.loc[i, "Unmutated %"] + results_df.loc[i, "Mutated %"]/2, str(mut), ha="center", va="center", color="white", fontsize=10)

    ax.set_ylabel("Percentage of Epitopes (%)")
    ax.set_title("Impact of Variant-Defining Mutations on Ancestral Spike CD4+ T Cell Epitopes")
    plt.xticks(rotation=45)
    plt.legend(title="Epitope Status")
    plt.tight_layout()
    plt.show()

# ============================================
# 5. HEATMAP: EPITOPE-RELEVANT MUTATIONS
# ============================================

    all_mutations = []
    for v in variants:
        all_mutations.extend(mutations_df[v].dropna().astype(str).tolist())
    all_mutations = sorted(set(all_mutations))

    presence_matrix = pd.DataFrame(0, index=all_mutations, columns=variants)
    for v in variants:
        muts = mutations_df[v].dropna().astype(str).tolist()
        for m in muts:
            if mutation_hits_epitope(m, epitopes_df):
                presence_matrix.loc[m, v] = 1

    cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", ["whitesmoke", "firebrick"])
    plt.figure(figsize=(10, len(all_mutations) * 0.3))
    sns.heatmap(
        presence_matrix,
        cmap=cmap,
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "Epitope-Relevant Mutation"},
        vmin=0, vmax=1
    )
    plt.title("Epitope-Relevant Mutations Across Variants")
    plt.xlabel("Variants")
    plt.ylabel("Mutations")
    plt.tight_layout()
    plt.savefig("epitope_relevant_mutations_heatmap.png", dpi=300, bbox_inches="tight")
    plt.show()

# ============================================
# 6. POSITIONAL MAPPING OF EPITOPES
# ============================================
def analyze_positional_mapping(epitope_file, variant_epitope_file, mutation_file, output_image="positional_mapping.png"):
    epitopes_df = pd.read_excel(epitope_file)
    epitopes_df.columns = epitopes_df.columns.str.strip()
    variant_epitopes_df = pd.read_excel(variant_epitope_file)
    variant_epitopes_df.columns = variant_epitopes_df.columns.str.strip()
    mutations_df = pd.read_excel(mutation_file)
    mutations_df.columns = mutations_df.columns.str.strip()
    variants = list(mutations_df.columns)

    mutated_present_dict = {}
    for variant in variants:
        mut_list = mutations_df[variant].dropna().astype(str).tolist()
        mutated, unmutated = get_mutated_epitopes(epitopes_df, mut_list)
        experimental_set = set(variant_epitopes_df[variant].dropna().astype(str).tolist())
        mutated_present_dict[variant] = mutated & experimental_set

    # Export
    output_df = pd.DataFrame(dict([(k, pd.Series(list(v))) for k, v in mutated_present_dict.items()]))
    output_df.to_excel("mutated_but_present_epitopes.xlsx", index=False)
    print("Exported mutated-but-present epitopes to 'mutated_but_present_epitopes.xlsx'")

    # Positional mapping visualization (with RBD highlight)

    plt.figure(figsize=(14, 7))

    # Define colors
    colors = {
        "conserved": "mediumseagreen",     # conserved epitopes
        "present": "gold",                 # mutated but still present
        "absent": "firebrick"              # mutated and absent
    }

    # Legend elements
    legend_elements = [
        Line2D([0], [0], color=colors["conserved"], lw=6, label="Conserved"),
        Line2D([0], [0], color=colors["present"], lw=6, label="Mutated but Present"),
        Line2D([0], [0], color=colors["absent"], lw=6, label="Mutated and Absent"),
    ]

    # --- RBD region highlight (aa 319–541) ---
    plt.axvspan(319, 541, color="lightcoral", alpha=0.2, label="RBD region")

    for i, variant in enumerate(variants):
        exp_epitopes = set(variant_epitopes_df[variant].dropna().astype(str).tolist())
        mutated_set, _ = get_mutated_epitopes(epitopes_df, mutations_df[variant].dropna().astype(str).tolist())

        for _, row in epitopes_df.iterrows():
            start, end, seq = row["Mapped Start Position"], row["Mapped End Position"], row["Epitope Sequence"]

            # Classification rules
            if seq not in mutated_set and seq in exp_epitopes:
                color = colors["conserved"]
            elif seq in mutated_set and seq in exp_epitopes:
                color = colors["present"]
            elif seq in mutated_set and seq not in exp_epitopes:
                color = colors["absent"]
            else:
                continue  # skip epitopes absent in both sets

            # Draw horizontal bars for each epitope
            plt.plot([start, end], [i, i], color=color, linewidth=20)

    # Y-axis labels = variants
    plt.yticks(range(len(variants)), variants)

    # Axis labels and title
    plt.xlabel("Ancestral Spike Protein Position (aa)")
    plt.ylabel("Variants")
    plt.title("Positional Mapping of Epitopes Across Variants", fontsize=14, weight="regular")

    # Legend placement (merged RBD region info automatically if same label)
    plt.legend(
        handles=legend_elements + [Line2D([0], [0], color="lightcoral", lw=6, alpha=0.5, label="RBD region")],
        loc='upper center',
        bbox_to_anchor=(0.5, 1.15),
        ncol=4,
        frameon=False
    )

    # Save + show
    plt.savefig("positional_mapping_with_RBD.png", dpi=300, bbox_inches="tight")
    plt.show()

# ============================================
# MAIN EXECUTION
# ============================================
if __name__ == "__main__":
    # Epitope overlaps
    analyze_epitope_overlap("variant_epitopes_only.xlsx", include_ancestral=False)
    analyze_epitope_overlap("variant_epitopes.xlsx", include_ancestral=True)

    # Epitope preservation and variant-specific
    analyze_epitope_preservation("variant_epitopes.xlsx", show_unmatched=True)

    # Shared mutations across variants
    analyze_shared_mutations("variant_mutations.xlsx")

    # Impact of variant-defining mutations
    analyze_variant_defining_mutations("epitopes.xlsx", "variant_mutations.xlsx")

    # Positional mapping
    analyze_positional_mapping("epitopes.xlsx", "variant_epitopes.xlsx", "variant_mutations.xlsx")
