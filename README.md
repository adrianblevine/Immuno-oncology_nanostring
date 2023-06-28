## Immuno-oncology_nanostring
Example scripts for NanoString immuno-oncology assay dataprocessing and Tumor Inflammation Signature (TIS) calculation. Note that this is not a self-contained workflow but rather representative scripts that could be used to develop a workflow for your own specific data, or applied to the supplementary data included with our paper.

For more information see the [preprint](https://www.researchsquare.com/article/rs-2655807/v1) for our paper _"Immuno-oncologic profiling of pediatric CNS tumors reveals major clinical significance of the tumor immune microenvironment"_ by Levine, Nobre, et al.

##### Workflow:
1. load\_RCC\_QC.R - Load raw RCC files and QC by geometric mean of housekeeping genes.  
2. normalize.R - Functions for normalization by housekeeping genes. 
3. calculate_TIS.R - Tumor Inflammation Signature calculation from normalized gene counts (such as the supplementary data provided with our paper for pediatric brain tumors). This could also be applied to any matrix of gene expression counts from other modalities, such as RNA sequencing.

##### Notes:
- example_tis_weights_file.csv - an example file is provided for how to format TIS weights, however the actual gene specific weights must be obtained from the patent filing WO2016094377 in claim 21c  (available at https://patents.google.com/patent/WO2016094377A1).   
