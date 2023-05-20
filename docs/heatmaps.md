---
layout: heatmaps
permalink: /RBD-heatmaps/
---

---

*NOTE data are preliminary*

### Overview

You can use this tool to explore the experimentally determined impacts of amino acid mutations on animal ACE2 binding in SARS-CoV-2 variant receptor-binding domains (RBD). 

#### Instruction

To use this tool, select the metric that you wish to display in the heatmap (change in affinty (delta-log10 $$K_A,app$$) for binding to ACE2 from white-tailed deer, Syrian hamster, little brown bat, cat, or human) by selecting that metric in the corresponding drop down menu. Hover over individual mutations to see exact numerical details.

#### Technical Details

The impact on ACE2 receptor-binding of every single amino-acid mutation in SARS-CoV-2 RBDs, as determined by high-throughput FACS-seq assays. Wildtype amino acids are indicated by an 'x', and gray squares indicate missing mutations from each library. The number of internally replicated barcodes with which a mutation was measured is visible as `Barcode Count` in the tooltips, where higher numbers indicate higher-confidence measurements.


### Data

Raw data  can be found [here](https://github.com/tstarrlab/SARS-CoV-2-RBD_DMS_animal-ACE2/blob/main/results/final_variant_scores/final_variant_scores.csv). The code used to make these plots can be found [here](https://github.com/tstarrlab/SARS-CoV-2-RBD_DMS_animal-ACE2/blob/main/RBD-Heatmaps-Interactive-Visualization.ipynb).
