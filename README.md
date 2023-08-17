## Immuno-oncology_nanostring
Example scripts for NanoString immuno-oncology assay data processing and Tumor Inflammation Signature (TIS) calculation. 

This includes a self-contained workflow with example test data and random TIS weights (see Notes below, actual weights must be obtained from patent filing). These scripts can then be further used develop a workflow for your own specific data, or applied to the supplementary data included with our paper.

For more information see the [preprint](https://www.researchsquare.com/article/rs-2655807/v1) for our paper _"Immuno-oncologic profiling of pediatric CNS tumors reveals major clinical significance of the tumor immune microenvironment"_ by Levine, Nobre, et al.

#### System requirements
- Tested on macOS Big Sur 11.7.7 (x86_64-apple-darwin17.0 (64-bit))
- R version 4.1.0 (2021-05-18)
- Packages: Biobase_2.54.0, BiocGenerics_0.40.0, NanoStringQCPro_1.24.0, here_1.0.1, janitor_2.1.0,
ggpubr_0.4.0, tidyverse_1.3.1 (includes forcats_0.5.1, stringr_1.4.0., dplyr_1.0.10, purrr_0.3.4, readr_2.1.1, tidyr_1.1.4, tibble_3.1.6, ggplot2_3.3.5)

###### Installation guide: 
- Download GitHub package using git clone https://github.com/adrianblevine/Immuno-oncology_nanostring
- Run from main folder
- Install time: 1 minute


#### Workflow:
1. load\_RCC\_QC.R - Load raw RCC files and QC by geometric mean of housekeeping genes.  
2. normalize.R - Functions for normalization by housekeeping genes. 
3. calculate_TIS.R - Tumor Inflammation Signature calculation from normalized gene counts (such as the supplementary data provided with our paper for pediatric brain tumors). This could also be applied to any matrix of gene expression counts from other modalities, such as RNA sequencing.

###### Expected run time: 1-2 minutes on normal desktop computer

###### Expected output: 
- results/eset_raw.rds : expression set of normalized counts
- results/eset_norm.rds : expression set of raw NanoString counts
- results/normalized_nanostring_counts.tsv : TSV file of normalized counts
- results/TIS_scores.tsv : TIS scores; note that these DO NOT use the actual weights, which must be obtained from the patent filing from the original authors (see Notes below)
- results/TIS_boxplot.pdf : example box plot of TIS scores by tumor type


#### Notes:
- example_tis_weights_file.csv - an example file is provided for how to format TIS weights, however the actual gene specific weights must be obtained from the patent filing WO2016094377 in claim 21c  (available at https://patents.google.com/patent/WO2016094377A1).   
