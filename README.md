# TMRvsBioLogic
This repository is associated with the manuscript: "Simple threshold-based Boolean rules fall short in capturing biological regulatory network dynamics".
## Contributor
- Priyotosh Sil


This repository contains 6 folders which are described below.

### 1. computational 
This folder contains the IMRs and BMRs corresponding to each sign combination. It also includes a list of one representative from each equivalence class of NCFs.

### 2. TMR_to_BF
This folder contains code for mapping TMRs (IMRs and BMRs) to the corresponding EFs (TMR_to_BF.ipynb).

### 3. Enrichment_analysis
This folder is associated with the enrichment analysis of TMRs in the three empirical datasets (BBM, MCBF and Harris)
- input: This folder contains the counts and theoretical fractions of different subtypes of TMRs within EThF (Expected_fraction_data_EF.tsv). The Empirical_dataset folder includes all three reference datasets used in our analysis. These files provide information about which classes each BF falls into.
- src: This folder contains code for performing enrichment analysis and generating plots.
- output: This folder contains the empirical fractions of various subtypes in the datasets.

<img src="schematic_fig_1_main.png">
