
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

# These figures were not used in the paper but are helpful to look at
plots <- DimPlot(Immune.integrated, combine = FALSE)
plots <- lapply(X = plots, 
                FUN = function(x) x + theme(legend.position = "top") + 
                  guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3))))
CombinePlots(plots)

# This analysis was not used but helpful
Immune.integrated <- BuildClusterTree(object = Immune.integrated, slot = "scale.data") 
Tool(object = Immune.integrated, slot = 'BuildClusterTree')
PlotClusterTree(Immune.integrated)

DimPlot(Immune.integrated, reduction = "umap", label = TRUE, label.size = 3)
DimPlot(Immune.integrated, reduction = "umap", label = TRUE, label.size = 3,split.by = "type")

markers.to.plot <- c("ERGIC1","ARAP2","PSMA6",
                     "AGMAT","DUSP6","TPT1",
                     "IL7R","HIST2H2AC"
)
DotPlot(Immune.integrated, features = rev(markers.to.plot), cols = c("blue","green","red"), dot.scale = 6, 
        split.by = "type", assay = 'SCT' ) + RotatedAxis()
DotPlot(Immune.integrated, features = rev(markers.to.plot),  dot.scale = 6, 
        assay = 'SCT' ) + RotatedAxis()

# Cluster identification

DefaultAssay(Immune.integrated) <- "RNA"

FeaturePlot(Immune.integrated, c("COL1A1", "ACTA2", "RAMP3", "s100A8", "GNLY", "ALB"))
VlnPlot(Immune.integrated, features = c("CYP3A7","SCD","SEC16B","BCHE","RPP25L","HPR","MGP","CCL14"), 
        pt.size = 0.2, ncol = 4)
VlnPlot(Immune.integrated, features = c("RAMP3","KRT7","ACTA2","S100A8","CD5L","CD2","GNLY","STMN1"), 
        pt.size = 0.2, ncol = 4)
VlnPlot(Immune.integrated, features = c("MS4A1","IGLC2","CD7","HBB"), 
        pt.size = 0.2, ncol = 4)
VlnPlot(Immune.integrated, features = c("CYP2A7","HMGCS1","SLBP","G6PC","HSD11B1","GSTA2","SPARCL1","CLEC1B"), 
        pt.size = 0.2, ncol = 4)
VlnPlot(Immune.integrated, features = c("INMT","KRT19","COL1A1","LYZ","MARCO","CD3D","PTGDS","HMGB2"), 
        pt.size = 0.2, ncol = 4)
VlnPlot(Immune.integrated, features = c("LTB","IGHG1","KLRB1","CA1"), 
        pt.size = 0.2, ncol = 4)
VlnPlot(Immune.integrated, features = c("GLUL","CYP2E1","ASS1","ASL","ALB","CYP2F2"), 
        pt.size = 0.2, ncol = 4)


# Label each Cluster

# Label is with reduced names
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

# Organize clusters the way I want
Idents(Immune.combined) <- factor(Idents(Immune.combined), levels = c("NK-cell","CD4-T-cell","CD8-T-cell",
                                                                      "Mp1", "Mp2","Mp3","Mp4","Mp5","Mp6","NonInf-Mp",
                                                                      "Hep1","Hep2","Hep3","Hep4","Hep5",
                                                                      "LSEC1","LSEC2",
                                                                      "Plasma","B-cell",
                                                                      "Chol","Eryth"))

# SAVE LABELS
Immune.combined$preserved <- Idents(Immune.combined)

# Make a new name with dashes (not underscore) of samples origins and cell type
Immune.combined$celltype_data <- paste(Immune.combined$type, 
                                       Immune.combined$stim, 
                                       Idents(Immune.combined), sep = "-")
Immune.combined$type_clus <- paste(Immune.combined$type, 
                                   Idents(Immune.combined), sep = "-")

# FIND ALL MARKERS - VERY SLOW
# These analyses were not used in the manuscript
pbmc.markers <- FindAllMarkers(Immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

write.table(pbmc.markers, file="ClusterMarkers_PBMCLiverscRNA_small.txt")

markers.to.plot <- c("CYP3A7","CYP2A7", "SCD","HMGCS1","SEC16B","SLBP","BCHE","G6PC",
                     "RPP25L","HSD11B1","HPR","GSTA2","MGP","SPARCL1","CCL14","CLEC1B",
                     "RAMP3","INMT","KRT7","KRT19","ACTA2","COL1A1","S100A8","LYZ",
                     "CD5L","MARCO","CD2","CD3D","GNLY","PTGDS","STMN1","HMGB2",
                     "MS4A1","LTB","IGLC2","IGHG1","CD7","KLRB1","HBB","CA1"
)

pdf("DotPlot_allcells.pdf", height = 10, width = 15, useDingbats=FALSE)
DotPlot(Immune.combined, features = rev(markers.to.plot), 
        cols = c("blue","red"), 
        dot.scale = 8, 
        assay = 'SCT'
) + RotatedAxis()
dev.off()

# Clustering of all cells with reduced cell labels
pdf("ClusteringFig2A_reduced.pdf", height = 5, width = 9, useDingbats=FALSE)
p <- DimPlot(Immune.combined, label = TRUE, label.size = 4)
p
dev.off()


########################################
##                                    ##
##             6. Subset              ##
##                                    ##
########################################

DefaultAssay(Immune_subset_v2) <- "SCT"

# Rename cells by mixing cell cluster with name

Immune.combined$cellstim_data <- paste(Immune.combined$stim, Idents(Immune.combined), sep = "-")

head(x = colnames(x = Immune.combined))
Immune.combined <- RenameCells(Immune.combined, add.cell.id = Immune.combined$cellstim_data)
head(x = colnames(x = Immune.combined))

# Subset to remove LPS

Immune.combined$metacomb <- paste(Immune.combined$type, 
                                  Immune.combined$stim, sep = "-")
Idents(Immune.combined) <-Immune.combined$metacomb
table(Idents(Immune.combined))
Immune_subset <- subset(Immune.combined, idents = c("Healthy-Basal","Alcohol-Basal","NA-NA"))

Idents(Immune_subset) <-Immune_subset$metacomb
table(Idents(Immune_subset))

Idents(Immune_subset) <-Immune_subset$preserved
table(Idents(Immune_subset))

# This is deprecated but the new option doesnt do it right
Immune_subset <- SubsetData(Immune_subset, 
                            ident.use=c("NK-cell","CD8-T-cell","CD4-T-cell",
                                        "Mp1", "Mp2","Mp3","Mp4","Mp5","Mp6","NonInf-Mp",
                                        "Hep1","Hep2","Hep3","Hep4","Hep5",
                                        "LSEC1","LSEC2","Plasma","B-cell","Chol"))


Idents(Immune_subset) <-Immune_subset$cellstim_data
table(Idents(Immune_subset))

Immune_subset_v2 <- SubsetData(Immune_subset, max.cells.per.ident = 200,
                               ident.remove=c("Basal-Chol",
                                              "Basal-LSEC2",
                                              "Basal-Hep5","Basal-Hep3","Basal-Hep1",
                                              "Basal-LSEC1","Basal-LSEC1",
                                              "Basal-Hep4","Basal-Hep2",
                                              "Basal-Plasma",
                                              "Basal-Mp6"))

table(Idents(Immune_subset_v2))
head(x = colnames(x = Immune_subset_v2))

expr <- GetAssayData(object = Immune_subset_v2, assay = "RNA", slot = "data")
expr <- as(Class = 'matrix', object = expr)
write.csv(x = expr, file = "PBMC_Liver_200_v2.csv", quote = FALSE)


# If there are issues, feel free to contact Adam - @atomadam2
