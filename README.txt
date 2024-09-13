This repository contains code for N-NMF, structured in a way to reproduce all results in our article "Normal tissue transcriptional signatures for tumor-type-agnostic phenotype prediction" - C. Weistuch et al.


All necessary inputs (data), relevant outputs, and functions are included in this repository except for the following exceptions which were too large for github:

1. CCLE_expression.csv
2. CCLE_mutations.csv
3. CCLE_gene_cn.csv
4. TCGA_BRCA_TPM.csv (as well as the corresponding files for PAAD and COAD).

Users can contact the lead author at weistucc@mskcc.org to access these files if needed.

Users should ensure that their data is TPM-normalized and not otherwise standardized.  Additional formatting and setup is detailed in Archetype_Test_Script.m


Files:

1. Archetype_Analysis_Script.m - reproduces all results from the manuscript.  This should only take a few minutes or so.
2. genesets.v7.5.1.txt - contains the pathway information and gene list queried from MSigDB for our analysis.
3. genesets.v2023.2.Hs.txt- contains the pathways information and gene list for all MSigDB Hallmark pathways.
4. GTEx_by_tissue.txt - contains the median TPM-normalized gene expression of each tissue in GTEx.
5. CCLE_expression.csv - contains TPM-normalized gene expression for the CCLE samples
6. CCLE_Drugs.csv - contains the drug sensitivities obtained from CCLE.
7. CCLE_mutations.csv - contains the mutation calls available with CCLE.
8. CCLE_gene_cn.csv - contains the gene-level CNA data available with CCLE.
9. sample_info_CCLE.csv - contains categorical information about each CCLE cell line.
10. TCGA_BRCA_TPM.csv - contains TPM-normalized gene expression for the breast cancer samples from TCGA.  (Analogous files for colorectal (COAD) and pancreatic (PAAD) data).
11. TCGA_BRCA_clin.csv - contains the clinical information (survival) for the breast cancer samples from TCGA.  (Analogous files for colorectal (COAD) and pancreatic (PAAD) data).
12. GSE155800_GIST_gene_TPM_matrix.txt - contains the TPM-normalized data for the GIST analysis.
13. GSE186901_palbo_geo_tpm_mat.txt - contains the TPM-normalized data for the palbociclib/MBC analysis.
14.QN_TPM_bc_nonlog_median.csv - contains the TPM-normalized data for the site-specific metastatic breast cancer analysis.
15. Impact.csv - contains the list of MSK-IMPACT genes
16. CNA_Heatmap.R - plots the copy number heatmap (Fig 5C)
17. bins.txt - contains the chromosome of each CNA bin (used in CNA_Heatmap.R)
18. gene_data.csv - contains the list of genes in each CNA bin
19. CNA_Table.Bin.txt - contains the table of bin-level Spearman correlations (used in CNA_Heatmap.R)
20. Plot_CNA.pdf - contains the sample output of CNA_Heatmap.R




21. NMF.m - computes the archetypes scores [H] and weights [W] from input data [V] and a specified number of archetypes [k].
22. NMF_New_Weights.m - computes new archetype scores [H] from input weights [W], data [V], and the appropriate number of archetypes [k].
23. Archetype_Computer.m - computes the archetype scores [H] associated with our trained archetypes from new data [Data].
24. Archetype_Test_Script.m - runs a test example and details how users should format their data.



Supplemental functions:
25. Clusterfunc.m - optimally permutes matrix entries to improve visualization.
26. MatSurv.m - computes and visualizes KM survival curves.

Creed et al., (2020). MatSurv: Survival analysis and visualization in MATLAB. Journal of Open Source Software, 5(46), 1830, https://doi.org/10.21105/joss.01830


Additional Outputs (in /Outputs):

1. Archetypes.csv - contains the trained archetype weights used in our analysis (780 genes x 6 archetypes).
2. S.csv - contains gene standard deviations used for data normalization normalization.
3. Entrez.csv - contains the Entrez Gene IDs for each of the 780 genes used.
4. Genes.csv - contains the gene names associated with 4.
5. Test_Archetype_Scores.csv - contains the computed archetype scores associated with the example provided in Archetype_Test_Script.m


Instructions:

1. To run the scripts, set the repository as the current MATLAB directory.  Then run Archetype_Analysis_Script.m to reproduce all of the analysis from the manuscript.
2. To run the example and learn how to format your data, run Archetype_Test_Script.m
3. To run the trained archetype model on new data, used Archetype_Computer.m
