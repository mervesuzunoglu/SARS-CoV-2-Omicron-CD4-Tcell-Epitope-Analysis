# SARS-CoV-2-Omicron-CD4-Tcell-Epitope-Analysis
This repository contains a Python-based bioinformatics pipeline developed to analyze the impact of SARS-CoV-2 mutations on CD4+ T cell epitopes, predict HLA Class II binding affinities using NetMHCIIpan 4.1, and visualize epitope conservation across variants (BA.2.86, JN.1, LP.8.1, NB.1.8.1, XEC)

**🔬 Research Workflow Summary: From Mutation to Conservation**
This repository provides a multi-stage bioinformatics framework designed to answer a central thesis question: To what extent do emergent Omicron sub-variants escape or preserve CD4+ T cell immunity compared to the ancestral strain?

The five pipelines work in a cascading logic:

Identification (Pipeline 01): We first map the raw mutational landscape, identifying exactly which epitopes are altered by amino acid substitutions in variants like XEC and JN.1.

Quantification (Pipeline 02): We then measure the "biological cost" of those mutations, calculating the fold-change in binding affinity (IC 
50
​	
  or percentile ranks) to see if the virus is actively weakening HLA-II binding.

Integration (Pipeline 03): Using UpSet plots and multi-variant intersections, we visualize the shared "immunological footprint" across the entire Omicron family, moving beyond simple one-to-one comparisons.

Localization (Pipeline 04): We specifically project these changes onto the Spike Protein structural domains. This reveals whether immune pressure is concentrated on critical areas like the RBD or if the S2 subunit remains a "safe haven" for T cell recognition.

Synthesis (Pipeline 05): Finally, we isolate the Conserved Core Immunome. This identifies the high-confidence peptide–HLA pairs that remain stable across all evolutionary shifts, serving as potential targets for future pan-sarbecovirus vaccines.
