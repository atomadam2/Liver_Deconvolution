# Liver_Deconvolution
These are the scripts used for my paper on liver deconvolution using single-cell RNA-seq data

This manuscripts consists of three parts:
1. Combining scRNA-seq data from healthy liver with scRNA-seq data from PBMCs.
2. Deconvolution of bulk RNA-seq data from chronic liver disease patients
3. WGCNA of bulk-RNA-seq data with clinical data and deconvolution results

Part 1: Combining scRNA-seq data from healthy liver with scRNA-seq data from PBMCs.

Data was obtained from:

scRNA-seq of five Healthy Livers - MacParland et al 2018 - GSE115469
These data were obtained from GSE as a raw count matrix .csv

scRNA-seq of 4 Healhy and 4 AH PBMCs - Kim et al 2020 - PRJNA596980
These data were generated by our lab, and analyzed as previously described
See git repository - atomadam2/PBMC_AH_LPS_scRNA-seq for details and code for alignment using CellRanger
After alignment, all data was uploaded into R using Seurat

LiverDeconv_scCombine.R - Script with code for pulling in all scRNA-seq data and combining
This script will create an R object containing all data. This RObj is needed for all other scripts.
Note: This will contain all PBMC data we generated from these samples, including cell treted with LPS. Subsequent analyses will remove this data for simplicity.

LiverDeconv_scFigures.R - Script with code for making all figures describing scRNA-seq data 
All parts of Figure 1 from paper are generated in this script

Part 2: Deconvolution of bulk RNA-seq data from chronic liver disease patients 

Reduce scRNA-seq data for deconvolution

LiverDeconv_Reduction.R - Script with code reducing scRNA-seq data
Data are reduced by:
Removing LPS treated cells
Creating single separate clusters for CD4 T-cells, CD8 T-cells, NK-cells
Removing clusters with very small numbers of cells (as described in paper)
Subsampling all clusters to just 200 cells/cluster.

