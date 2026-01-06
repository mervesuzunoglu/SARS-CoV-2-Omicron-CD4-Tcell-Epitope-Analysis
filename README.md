# SARS-CoV-2-Omicron-CD4-Tcell-Epitope-Analysis
This repository contains a Python-based bioinformatics pipeline developed to analyze the impact of SARS-CoV-2 mutations on CD4+ T cell epitopes, predict HLA Class II binding affinities using NetMHCIIpan 4.1, and visualize epitope conservation across variants (BA.2.86, JN.1, LP.8.1, NB.1.8.1, XEC)

**🔬 Research Workflow Summary: From Mutation to Conservation**

This repository provides a multi-stage bioinformatics framework developed to evaluate the extent to which emergent Omicron sub-variants escape or preserve CD4+ T cell immunity relative to the ancestral strain.

The analytical process is structured in a cascading logic across five distinct pipelines:

Mutational Mapping (Pipeline 01): The raw mutational landscape is characterized, identifying specific epitopes altered by amino acid substitutions in lineages such as XEC and JN.1.

Affinity Quantification (Pipeline 02): The biological impact of identified mutations is measured by calculating the fold-change in binding affinity (percentile ranks), determining if viral evolution results in the attenuation of HLA-II binding.

Cross-Variant Integration (Pipeline 03): The shared immunological footprint across the Omicron phylogeny is visualized through UpSet plots and multi-variant intersections, allowing for the identification of overlapping epitope sets.

Structural Localization (Pipeline 04): Epitope loss is projected onto specific Spike protein structural domains. This analysis determines if immune evasion is concentrated within the Receptor Binding Domain (RBD) or if the S2 subunit maintains stability in T cell recognition.

Core Immunome Synthesis (Pipeline 05): The Universal Conserved Core is isolated, identifying high-confidence peptide–HLA pairs that remain functionally intact across all analyzed variants. These stable targets are highlighted for their potential relevance in pan-sarbecovirus vaccine design.
