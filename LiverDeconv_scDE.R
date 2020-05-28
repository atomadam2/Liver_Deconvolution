
########################################
##                                    ##
##   1. Install programs needed       ##
##                                    ##
########################################

install.packages('Seurat')

# There version of Seurat used here is 3.1.4

library(Seurat)
# Needs python
library(reticulate)
py_config()
# create a new environment 
conda_create("r-reticulate")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

### You will need this to update R ###
#install.packages("installr")
#install.packages("stringr")
#library(installr)
#install.packages("dplyr")
#install.packages("cowplot")
#updateR()

########################################
##                                    ##
##         2. Load packages           ##
##                                    ##
########################################

library(Seurat)
library(cowplot)
library(dplyr)
library(monocle3)
library(sctransform)
library(Matrix)
library(ggplot2)
library(ggdendro)
library(plyr)
library(reshape2)
library(scales)
library(viridis)
library(viridisLite)
library(heatmaply)

# Set a working directory to have all files

setwd("~/XXXXXX/")
getwd()

########################################
##                                    ##
##       3. Clustering, 2 ways        ##
##                                    ##
########################################

# Get RObj from the other script: LiverDeconv_scCombine.R

load("XXXXX/PBMC_Liver.integrated_20200310.Robj")

DefaultAssay(Immune.integrated) <- "integrated"

Immune.integrated <- RunPCA(Immune.integrated, verbose = FALSE)
Immune.integrated <- RunUMAP(Immune.integrated, dims = 1:30)
Immune.integrated <- FindNeighbors(Immune.integrated)
Immune.integrated <- FindClusters(Immune.integrated,  resolution = 2.0)

# These figures were not published but are helpfu to look at
plots <- DimPlot(Immune.integrated, combine = FALSE)
plots <- lapply(X = plots, 
                FUN = function(x) x + theme(legend.position = "top") + 
                  guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3))))
CombinePlots(plots)

DimPlot(Immune.integrated, reduction = "umap", label = TRUE, label.size = 3)
DimPlot(Immune.integrated, reduction = "umap", label = TRUE, label.size = 3,split.by = "type")

# Label all clusters by predicted identity
# This is the reduced names, with T-cells combined
Immune.combined <- RenameIdents(Immune.integrated, 
                                `0`="NK-cell",`1`="NK-cell", 
                                `2`="CD4-T-cell",`3`="CD4-T-cell",`4`="Mp1", 
                                `5`="CD4-T-cell",`6`="CD8-T-cell",`7`="Hep1",
                                `8`="CD4-T-cell",`9`="CD8-T-cell",
                                `10`="CD4-T-cell",`11`="Mp2",`12`="CD8-T-cell",
                                `13`="Hep2",`14`="Hep3",
                                `15`="CD4-T-cell",`16`="CD4-T-cell",
                                `17`="Mp3",`18`="NonInf-Mp",
                                `19`="Hep4",`20`="CD4-T-cell",
                                `21`="Plasma",`22`="B-cell",
                                `23`="LSEC1",`24`="Mp4",
                                `25`="Hep5",`26`="Mp5",
                                `27`="LSEC2",`28`="Chol",
                                `29`="CD4-T-cell",`30`="Eryth",`31`="Mp6"
)
DimPlot(Immune.combined, label = TRUE, label.size = 4)

# Organize labels the way I want it
Idents(Immune.combined) <- factor(Idents(Immune.combined), levels = c("NK-cell","CD4-T-cell","CD8-T-cell",
                                                                      "Mp1", "Mp2","Mp3","Mp4","Mp5","Mp6","NonInf-Mp",
                                                                      "Hep1","Hep2","Hep3","Hep4","Hep5",
                                                                      "LSEC1","LSEC2",
                                                                      "Plasma","B-cell",
                                                                      "Chol","Eryth"))

# Create new metadata label
Immune.combined$type_clus <- paste(Immune.combined$type, 
                                   Idents(Immune.combined), sep = "-")
Idents(Immune.combined) <- Immune.combined$type_clus

# This for loop will create separate files for DE genes for all comparisons for all clusters:
# Healthy Liver vs Healthy PBMC
# Healthy PBMC vs AH PBMC
# AH PBMC vs Healthy Liver
# Does not include liver specific cells (Hep, Chol, LSEC, Plasma cells)

DefaultAssay(Immune.combined) <- "SCT"
for (spec_cell in c("B-cell","NK-cell","CD4-T-cell","CD8-T-cell",
                    "NonInf-Mp","Mp1", "Mp2","Mp3","Mp4","Mp5"
)){
  
  Healthy <- (paste( "Healthy-", spec_cell, sep = ""))
  Alcohol <- (paste("Alcohol-", spec_cell, sep = ""))
  Resident <- (paste("NA-", spec_cell, sep = ""))
  
  HvA <- FindMarkers(Immune.combined, ident.1 = Healthy, ident.2 = Alcohol, verbose = FALSE, test.use = "LR")
  HvA_Comp_File <- (paste(spec_cell, "_Healthy_Alcohol_ResvPerip_scRNA.txt", sep = ""))
  write.table(HvA, file=HvA_Comp_File)
  
  HvL <- FindMarkers(Immune.combined, ident.1 = Healthy, ident.2 = Resident, verbose = FALSE, test.use = "LR")
  HvL_Comp_File <- (paste(spec_cell, "_Healthy_Resident_ResvPerip_scRNA.txt", sep = ""))
  write.table(HvL, file=HvL_Comp_File)  
  
  AvL <- FindMarkers(Immune.combined, ident.1 = Alcohol, ident.2 = Resident, verbose = FALSE, test.use = "LR")
  AvL_Comp_File <- (paste(spec_cell, "_Alcohol_Resident_ResvPerip_scRNA.txt", sep = ""))
  write.table(AvL, file=AvL_Comp_File)
}

# If there are issues, feel free to contact Adam - @atomadam2