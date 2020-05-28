
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

# These figures were not published but are helpful to look at
plots <- DimPlot(Immune.integrated, combine = FALSE)
plots <- lapply(X = plots, 
                FUN = function(x) x + theme(legend.position = "top") + 
                  guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3))))
CombinePlots(plots)

# This is bulding a cluster tree. This analysis was also not used but helpful to look at
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

# These markers are from a few papers. This is not an exhaustive list of cell type markers, but a helpful place to start
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

# Each Cluster named separately
Immune.combined <- RenameIdents(Immune.integrated, 
                                `0`="NK-cell1",`1`="NK-cell2", 
                                `2`="CD4-T1",`3`="CD4-T2",`4`="Mp1", 
                                `5`="CD4-T3",`6`="CD8-T2",`7`="Hep1",
                                `8`="CD4-T4",`9`="CD8-T3",
                                `10`="CD4-T5",`11`="Mp2",`12`="CD8-T1",
                                `13`="Hep2",`14`="Hep3",
                                `15`="CD4-T6",`16`="CD4-T7",
                                `17`="Mp3",`18`="NonInf-Mp",
                                `19`="Hep4",`20`="CD4-T8",
                                `21`="Plasma",`22`="B-cell",
                                `23`="LSEC1",`24`="Mp4",
                                `25`="Hep5",`26`="Mp5",
                                `27`="LSEC2",`28`="Chol",
                                `29`="CD4-T9",`30`="Eryth",`31`="Mp6"
)
DimPlot(Immune.combined, label = TRUE, label.size = 4)

# Organize clusters the way I want
Idents(Immune.combined) <- factor(Idents(Immune.combined), levels = c("NK-cell1","NK-cell2",
                                                                      "CD8-T1","CD8-T2","CD8-T3",
                                                                      "CD4-T1","CD4-T2","CD4-T3","CD4-T4",
                                                                      "CD4-T5","CD4-T6","CD4-T7","CD4-T8","CD4-T9",
                                                                      "Mp1","Mp2","Mp3","Mp4","Mp5","Mp6","NonInf-Mp",
                                                                      "Hep1","Hep2","Hep3","Hep4","Hep5",
                                                                      "LSEC1","LSEC2",
                                                                      "Chol","B-cell","Plasma",
                                                                      "Eryth"))
table (Idents(Immune.combined))
table(Immune.combined$orig.ident)
table(Immune.combined$type)

VlnPlot(Immune.combined, features = c("HMGCS1", "SCD",  "G6PC", "BCHE"), 
        pt.size = 0.2, ncol =2)
FeaturePlot(Immune.combined,  c("MFSD2A","SOX9","KRT7","KRT19"))
FeaturePlot(Immune.combined,  c("GLUL", "ASS1", "CYP2E1", "CYP2F2"))
FeaturePlot(Immune.combined,  c("SCD", "HMGCS1", "G6PC", "BCHE"))


# Make Figures in paper including Supplement and First clustering analysis
DefaultAssay(Immune.combined) <- "SCT"

pdf("ZonationMarkersSuppFig1.pdf", height = 10, width = 10, useDingbats=FALSE)
FeaturePlot(Immune.combined,  c("HMGCS1", "SCD",  "G6PC", "BCHE"))
dev.off()

pdf("ZonationMarkersSuppFig1_LSEC.pdf", height = 10, width = 10, useDingbats=FALSE)
FeaturePlot(Immune.combined,  c("SPARCL1", "CLEC1B","SCD", "BCHE"))
dev.off()

pdf("ZonationMarkersSuppFig1_Hep.pdf", height = 10, width = 10, useDingbats=FALSE)
FeaturePlot(Immune.combined,  c("GHR", "HAMP","CYP2A7", "CYP2E1"))
dev.off()

pdf("Clustering_allFig1A.pdf", height = 5, width = 9, useDingbats=FALSE)
p <- DimPlot(Immune.combined, label = TRUE, label.size = 4)
p
dev.off()

pdf("ClusteringFig1B_split.pdf", height = 5, width = 15, useDingbats=FALSE)
p <- DimPlot(Immune.combined, reduction = "umap", label = FALSE, split.by = "type")
p
dev.off()

#Create table with all cell numbers for all samples
Immune.combined$celltype_origin <- paste(Immune.combined$type, 
                                         Immune.combined$stim, 
                                         Idents(Immune.combined), sep = "_")
cell.num <- as.data.frame(table(Immune.combined$celltype_origin))
write.table(cell.num, file="Cells_PBMC_Liver.txt")


# Do this with less cell labels (T-cells combined)
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

Idents(Immune.combined) <- factor(Idents(Immune.combined), levels = c("NK-cell","CD4-T-cell","CD8-T-cell",
                                                                      "Mp1", "Mp2","Mp3","Mp4","Mp5","Mp6","NonInf-Mp",
                                                                      "Hep1","Hep2","Hep3","Hep4","Hep5",
                                                                      "LSEC1","LSEC2",
                                                                      "Plasma","B-cell",
                                                                      "Chol","Eryth"))

# SAVE LABELS AND ESTABLISH NEW LABELS
Immune.combined$preserved <- Idents(Immune.combined)

# Make a new name with dashes (not underscore) of samples origins and cell type
Immune.combined$celltype_data <- paste(Immune.combined$type, 
                                       Immune.combined$stim, 
                                       Idents(Immune.combined), sep = "-")

Immune.combined$type_clus <- paste(Immune.combined$type, 
                                   Idents(Immune.combined), sep = "-")

# FIND ALL MARKERS - VERY SLOW
pbmc.markers <- FindAllMarkers(Immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

write.table(pbmc.markers, file="ClusterMarkers_PBMCLiverscRNA_small.txt")

# This figure was also not used in the final manuscript
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

# This is the second clustering analysis figure
pdf("ClusteringFig2A_reduced.pdf", height = 5, width = 9, useDingbats=FALSE)
p <- DimPlot(Immune.combined, label = TRUE, label.size = 4)
p
dev.off()

########################################
##                                    ##
##       4. Just Myeloid cells        ##
##                                    ##
########################################

# More specific analyses looking at the myeloid cells but subseting them
All_Mono <- subset(Immune.combined, idents = c("Mp1","Mp2","Mp3","Mp4","Mp5","NonInf-Mp"))
Idents(All_Mono) <- All_Mono$celltype_data

# ORganize clusters the way I want
All_Mono <- subset(All_Mono, idents = c('Healthy-Basal-Mp1','Healthy-Basal-Mp2','Healthy-Basal-Mp3','Healthy-Basal-Mp4',
                                        'Healthy-Basal-Mp5','Healthy-Basal-NonInf-Mp',
                                        'Alcohol-Basal-Mp1','Alcohol-Basal-Mp2','Alcohol-Basal-Mp3','Alcohol-Basal-Mp4',
                                        'Alcohol-Basal-Mp5','Alcohol-Basal-NonInf-Mp'))
Idents(All_Mono) <- factor(Idents(All_Mono), levels = c('Healthy-Basal-Mp1','Healthy-Basal-Mp2','Healthy-Basal-Mp3','Healthy-Basal-Mp4',
                                                        'Healthy-Basal-Mp5','Healthy-Basal-NonInf-Mp',
                                                        'Alcohol-Basal-Mp1','Alcohol-Basal-Mp2','Alcohol-Basal-Mp3','Alcohol-Basal-Mp4',
                                                        'Alcohol-Basal-Mp5','Alcohol-Basal-NonInf-Mp'))

# Change data used
DefaultAssay(All_Mono) <- "SCT"
Idents(All_Mono) <- All_Mono$preserved

# Find markers for every cluster compared to all remaining cells, report only the positive ones
Mono.markers <- FindAllMarkers(All_Mono, only.pos = TRUE, min.pct = 0.25, min.diff.pct = 0.2, logfc.threshold = 0.25)
allBest <- Mono.markers %>% group_by(cluster) %>% top_n(n = 8, wt = avg_logFC)

write.table(Mono.markers, file="ClusterMarkers_PBMCLiverscRNA_Mp.txt")

# Create the dotplot figure of myeloid specific genes
markers.to.plot <- c(allBest$gene)
markers.to.plot <- markers.to.plot[!duplicated(markers.to.plot)]

pdf("DotPlotFig1C.pdf", height = 10, width = 15, useDingbats=FALSE)
DotPlot(All_Mono, features = rev(markers.to.plot), 
        cols = c("blue","red"), 
        dot.scale = 8, 
        assay = 'SCT'
) + RotatedAxis()+ scale_colour_gradient2(low = "blue", mid = "white", high = "red")
dev.off()


########################################
##                                    ##
##    4. Just Hepatocytes cells       ##
##                                    ##
########################################

# These Hepatocyte analyses were not used in the manuscript

All_Hep <- subset(Immune.combined, idents = c("Hep1","Hep2","Hep3","Hep4","Hep5"))

Idents(All_Hep) <- All_Hep$celltype_data

All_Hep <- subset(All_Hep, idents = c('NA-NA-Hep1','NA-NA-Hep2','NA-NA-Hep3','NA-NA-Hep4',
                                      'NA-NA-Hep5'))


Idents(All_Hep) <- factor(Idents(All_Hep), levels = c('NA-NA-Hep1','NA-NA-Hep2','NA-NA-Hep3','NA-NA-Hep4',
                                                      'NA-NA-Hep5'))

DefaultAssay(All_Hep) <- "SCT"

Idents(All_Hep) <- All_Hep$preserved

# find markers for every cluster compared to all remaining cells, report only the positive ones
Hep.markers <- FindAllMarkers(All_Hep, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
allBest <- Hep.markers %>% group_by(cluster) %>% top_n(n = 8, wt = avg_logFC)

markers.to.plot <- c(allBest$gene)
markers.to.plot <- markers.to.plot[!duplicated(markers.to.plot)]

pdf("DotPlot_Hep.pdf", height = 10, width = 15, useDingbats=FALSE)
DotPlot(All_Hep, features = rev(markers.to.plot), 
        cols = c("blue","red"), 
        dot.scale = 8, 
        assay = 'SCT'
) + RotatedAxis()+ scale_colour_gradient2(low = "blue", mid = "white", high = "red")
dev.off()

# find markers for every cluster compared to all remaining cells, report only the positive ones
Hep.markers <- FindAllMarkers(All_Hep, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Hep.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

write.table(Hep.markers, file="ClusterMarkers_PBMCLiverscRNA_Hep.txt")

VlnPlot(All_Hep, features = c("LEF1","CCR7","GIMAP"), 
        pt.size = 0.2, ncol = 4)

markers.to.plot <- c(
  "TPT1","CCR7","GIMAP7","LEF1"
)

pdf("DotPlot_HepMarkers.pdf", height = 10, width = 15, useDingbats=FALSE)
DotPlot(All_Hep, features = rev(markers.to.plot), 
        cols = c("blue","red"), 
        dot.scale = 8, 
        assay = 'SCT'
) + RotatedAxis()
dev.off()


########################################
##                                    ##
##    5. Violin Plots of Expression   ##
##                                    ##
########################################

# Store the labels
Immune.combined$preserved <- Idents(Immune.combined)

# Make a new name with dashes (not underscore) of samples origins and cell type
Immune.combined$celltype_data <- paste(Immune.combined$type, 
                                       Immune.combined$stim, 
                                       Idents(Immune.combined), sep = "-")
Immune.combined$type_clus <- paste(Immune.combined$type, 
                                   Idents(Immune.combined), sep = "-")
Idents(Immune.combined) <- Immune.combined$type_clus

# ORganize the clusters the way I want
Idents(Immune.combined) <- factor(Idents(Immune.combined), levels = c('Healthy-Mp1','Healthy-Mp2','Healthy-Mp3','Healthy-Mp4',
                                                                      'Healthy-Mp5','Healthy-Mp6','Healthy-NonInf-Mp',
                                                                      'Healthy-Hep1','Healthy-Hep2','Healthy-Hep3','Healthy-Hep4','Healthy-Hep5',
                                                                      'Alcohol-Mp1','Alcohol-Mp2','Alcohol-Mp3','Alcohol-Mp4',
                                                                      'Alcohol-Mp5','Alcohol-Mp6','Alcohol-NonInf-Mp',
                                                                      'Alcohol-Hep1','Alcohol-Hep2','Alcohol-Hep3','Alcohol-Hep4','Alcohol-Hep5',
                                                                      'NA-Mp1','NA-Mp2','NA-Mp3','NA-Mp4','NA-Mp5','NA-Mp6','NA-NonInf-Mp',
                                                                      'NA-Hep1','NA-Hep2','NA-Hep3','NA-Hep4','NA-Hep5'))

# Make all violin plots used in paper, this is long and exhaustive. I should have made a function
DefaultAssay(Immune.combined) <- "SCT"

pdf("ViolinFig1D_CD1416.pdf", height = 7, width = 7, useDingbats=FALSE)
plots <- VlnPlot(Immune.combined, features = c("CD14", "FCGR3A"), 
                 pt.size = 0, combine = FALSE, idents = c('Healthy-Mp1','Healthy-Mp2',
                                                          'Healthy-Mp3','Healthy-Mp4',
                                                          'Healthy-Mp5',
                                                          'Healthy-NonInf-Mp',
                                                          'Alcohol-Mp1','Alcohol-Mp2',
                                                          'Alcohol-Mp3','Alcohol-Mp4',
                                                          'Alcohol-Mp5',
                                                          'Alcohol-NonInf-Mp',
                                                          'NA-Mp1','NA-Mp2','NA-Mp3',
                                                          'NA-Mp4','NA-Mp5',
                                                          'NA-NonInf-Mp'
                 ))
plots <- lapply(
  X = plots,
  FUN = function(p) p 
  + theme(legend.position = 'none') 
  + ggplot2::scale_fill_manual(values = c('black','black','black','black','black','black',
                                          'white','white','white','white','white','white',
                                          'blue','blue','blue','blue','blue','blue'))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()

pdf("ViolinFig1D_CD11BCD71.pdf", height = 7, width = 7, useDingbats=FALSE)
plots <- VlnPlot(Immune.combined, features = c("ITGAM", "TFRC"), 
                 pt.size = 0, combine = FALSE, idents = c('Healthy-Mp1','Healthy-Mp2',
                                                          'Healthy-Mp3','Healthy-Mp4',
                                                          'Healthy-Mp5',
                                                          'Healthy-NonInf-Mp',
                                                          'Alcohol-Mp1','Alcohol-Mp2',
                                                          'Alcohol-Mp3','Alcohol-Mp4',
                                                          'Alcohol-Mp5',
                                                          'Alcohol-NonInf-Mp',
                                                          'NA-Mp1','NA-Mp2','NA-Mp3',
                                                          'NA-Mp4','NA-Mp5',
                                                          'NA-NonInf-Mp'
                 ))
plots <- lapply(
  X = plots,
  FUN = function(p) p 
  + theme(legend.position = 'none') 
  + ggplot2::scale_fill_manual(values = c('black','black','black','black','black','black',
                                          'white','white','white','white','white','white',
                                          'blue','blue','blue','blue','blue','blue'))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()

pdf("ViolinFig1D_LYZMARCO.pdf", height = 7, width = 7, useDingbats=FALSE)
plots <- VlnPlot(Immune.combined, features = c("LYZ", "MARCO"), 
                 pt.size = 0, combine = FALSE, idents = c('Healthy-Mp1','Healthy-Mp2',
                                                          'Healthy-Mp3','Healthy-Mp4',
                                                          'Healthy-Mp5',
                                                          'Healthy-NonInf-Mp',
                                                          'Alcohol-Mp1','Alcohol-Mp2',
                                                          'Alcohol-Mp3','Alcohol-Mp4',
                                                          'Alcohol-Mp5',
                                                          'Alcohol-NonInf-Mp',
                                                          'NA-Mp1','NA-Mp2','NA-Mp3',
                                                          'NA-Mp4','NA-Mp5',
                                                          'NA-NonInf-Mp'
                 ))
plots <- lapply(
  X = plots,
  FUN = function(p) p 
  + ggplot2::scale_y_continuous(breaks = c(0,1,2,3))
  + theme(legend.position = 'none') 
  + ggplot2::scale_fill_manual(values = c('black','black','black','black','black','black',
                                          'white','white','white','white','white','white',
                                          'blue','blue','blue','blue','blue','blue'))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()

########################################
##                                    ##
##             6. Subset              ##
##                                    ##
########################################

# Subset the data to look at resident vs peripheral cell gene expression and compare

# Differential expression measurements were made in a different script: LiverDeconv_scDE.R

# Rename cells by mixing cell cluster with name
head(x = colnames(x = Immune.combined))
Immune.combined <- RenameCells(Immune.combined, add.cell.id = Immune.combined$celltype_data)
head(x = colnames(x = Immune.combined))

# Subset to remove LPS
Immune.combined$metacomb <- paste(Immune.combined$type, 
                                  Immune.combined$stim, sep = "-")
Idents(Immune.combined) <-Immune.combined$metacomb
table(Idents(Immune.combined))
Immune_subset <- subset(Immune.combined, idents = c("Healthy-Basal","Alcohol-Basal","NA-NA"))

# Get saved labels
Idents(Immune_subset) <-Immune_subset$preserved
table(Idents(Immune_subset))

# This is deprecated but the new option doesnt do it right
Immune_subset <- SubsetData(Immune_subset, 
                            ident.use=c("NK-cell","CD8-T-cell","CD4-T-cell",
                                        "Mp1", "Mp2","Mp3","Mp4","Mp5","Mp6","NonInf-Mp",
                                        "Hep1","Hep2","Hep3","Hep4","Hep5",
                                        "LSEC1","LSEC2","Plasma","B-cell","Chol"))

table(Idents(Immune_subset))
head(x = colnames(x = Immune_subset))

Idents(Immune_subset) <-Immune_subset$celltype_data
table(Idents(Immune_subset))

# Remove clusters with very few cells
Immune_subset_v1 <- SubsetData(Immune_subset, 
                               ident.remove=c("Alcohol-Basal-Chol","Healthy-Basal-Chol",
                                              "Healthy-Basal-LSEC2","Alcohol-Basal-LSEC2",
                                              "Alcohol-Basal-Hep5","Healthy-Basal-Hep3","Healthy-Basal-Hep1",
                                              "Alcohol-Basal-Hep3","Alcohol-Basal-Hep1","Alcohol-Basal-Hep2",
                                              "Healthy-Basal-LSEC1","Alcohol-Basal-LSEC1",
                                              "Healthy-Basal-Hep4","Alcohol-Basal-Hep4",
                                              "Healthy-Basal-Plasma","Alcohol-Basal-Plasma",
                                              "Healthy-Basal-Mp6","Alcohol-Basal-Mp6"))

table(Idents(Immune_subset_v1))

# ORganize clusters the way I want
Idents(Immune_subset_v1) <- factor(Idents(Immune_subset_v1), levels = c("Healthy-B-cell","Alcohol-B-cell","NA-B-cell",
                                                                        "Healthy-NK-cell","Alcohol-NK-cell","NA-NK-cell",
                                                                        "Healthy-CD4-T-cell","Alcohol-CD4-T-cell","NA-CD4-T-cell",
                                                                        "Healthy-CD8-T-cell","Alcohol-CD8-T-cell","NA-CD8-T-cell",
                                                                        "Healthy-Mp1","Alcohol-Mp1","NA-Mp1",
                                                                        "Healthy-Mp2","Alcohol-Mp2","NA-Mp2",
                                                                        "Healthy-Mp3","Alcohol-Mp3","NA-Mp3",
                                                                        "Healthy-Mp4","Alcohol-Mp4","NA-Mp4",
                                                                        "Healthy-Mp5","Alcohol-Mp5","NA-Mp5", 
                                                                        "Healthy-NonInf-Mp","Alcohol-NonInf-Mp","NA-NonInf-Mp",
                                                                        "NA-Hep1","NA-Hep2","NA-Hep3","NA-Hep4","NA-Hep5",
                                                                        "NA-LSEC1","NA-LSEC2","NA-Chol","NA-Plasma"))

DefaultAssay(Immune_subset_v1) <- "SCT"

# look at a few notable gennes (Hepatocyte specific genes)
Idents(Immune_subset_v1) <-Immune_subset_v1$type_clus
plots <- VlnPlot(Immune_subset_v1, features = c("SERPINA1","HP","ALB"),  
                 pt.size = 0, combine = FALSE )
CombinePlots(plots = plots, ncol = 1)

# Create all figures showing hepatocyte gene expression in liver cells

pdf("ViolinFSuppFig_ALB_HP.pdf", height = 7, width = 7, useDingbats=FALSE)
plots <- VlnPlot(Immune_subset_v1, features = c("ALB","HP"), 
                 pt.size = 0, combine = FALSE, idents = c("Healthy-B-cell","Alcohol-B-cell","NA-B-cell",
                                                          "Healthy-NK-cell","Alcohol-NK-cell","NA-NK-cell",
                                                          "Healthy-CD4-T-cell","Alcohol-CD4-T-cell","NA-CD4-T-cell",
                                                          "Healthy-CD8-T-cell","Alcohol-CD8-T-cell","NA-CD8-T-cell",
                                                          "Healthy-Mp1","Alcohol-Mp1","NA-Mp1",
                                                          "Healthy-Mp2","Alcohol-Mp2",
                                                          "Healthy-Mp3","Alcohol-Mp3","NA-Mp3",
                                                          "Healthy-Mp4","Alcohol-Mp4",
                                                          "Healthy-Mp5","Alcohol-Mp5","NA-Mp5", 
                                                          "Healthy-NonInf-Mp","Alcohol-NonInf-Mp","NA-NonInf-Mp",
                                                          "NA-Hep1","NA-Hep2","NA-Hep3","NA-Hep4","NA-Hep5",
                                                          "NA-LSEC1","NA-LSEC2","NA-Chol","NA-Plasma"
                 ))
plots <- lapply(
  X = plots,
  FUN = function(p) p 
  + ggplot2::scale_y_continuous(breaks = c(0,1,2,3))
  + theme(legend.position = 'none') 
  + ggplot2::scale_fill_manual(values = c('black','white','blue',
                                          'black','white','blue',
                                          'black','white','blue',
                                          'black','white','blue',
                                          'black','white','blue',
                                          'black','white',
                                          'black','white','blue',
                                          'black','white',
                                          'black','white','blue',
                                          'black','white','blue',
                                          'blue','blue','blue','blue','blue',
                                          'blue','blue','blue','blue'))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()

pdf("ViolinFSuppFig_TF_APOA2.pdf", height = 7, width = 7, useDingbats=FALSE)
plots <- VlnPlot(Immune_subset_v1, features = c("TF","APOA2"), 
                 pt.size = 0, combine = FALSE, idents = c("Healthy-B-cell","Alcohol-B-cell","NA-B-cell",
                                                          "Healthy-NK-cell","Alcohol-NK-cell","NA-NK-cell",
                                                          "Healthy-CD4-T-cell","Alcohol-CD4-T-cell","NA-CD4-T-cell",
                                                          "Healthy-CD8-T-cell","Alcohol-CD8-T-cell","NA-CD8-T-cell",
                                                          "Healthy-Mp1","Alcohol-Mp1","NA-Mp1",
                                                          "Healthy-Mp2","Alcohol-Mp2",
                                                          "Healthy-Mp3","Alcohol-Mp3","NA-Mp3",
                                                          "Healthy-Mp4","Alcohol-Mp4",
                                                          "Healthy-Mp5","Alcohol-Mp5","NA-Mp5", 
                                                          "Healthy-NonInf-Mp","Alcohol-NonInf-Mp","NA-NonInf-Mp",
                                                          "NA-Hep1","NA-Hep2","NA-Hep3","NA-Hep4","NA-Hep5",
                                                          "NA-LSEC1","NA-LSEC2","NA-Chol","NA-Plasma"
                 ))
plots <- lapply(
  X = plots,
  FUN = function(p) p 
  + ggplot2::scale_y_continuous(breaks = c(0,1,2,3))
  + theme(legend.position = 'none') 
  + ggplot2::scale_fill_manual(values = c('black','white','blue',
                                          'black','white','blue',
                                          'black','white','blue',
                                          'black','white','blue',
                                          'black','white','blue',
                                          'black','white',
                                          'black','white','blue',
                                          'black','white',
                                          'black','white','blue',
                                          'black','white','blue',
                                          'blue','blue','blue','blue','blue',
                                          'blue','blue','blue','blue'))
)
CombinePlots(plots = plots, ncol = 1)
dev.off()

plots <- VlnPlot(Immune_subset_v1, features = c("FGA","FGG"), 
                 pt.size = 0, combine = FALSE, idents = c("Healthy-B-cell","Alcohol-B-cell","NA-B-cell",
                                                          "Healthy-NK-cell","Alcohol-NK-cell","NA-NK-cell",
                                                          "Healthy-CD4-T-cell","Alcohol-CD4-T-cell","NA-CD4-T-cell",
                                                          "Healthy-CD8-T-cell","Alcohol-CD8-T-cell","NA-CD8-T-cell",
                                                          "Healthy-Mp1","Alcohol-Mp1","NA-Mp1",
                                                          "Healthy-Mp2","Alcohol-Mp2",
                                                          "Healthy-Mp3","Alcohol-Mp3","NA-Mp3",
                                                          "Healthy-Mp4","Alcohol-Mp4",
                                                          "Healthy-Mp5","Alcohol-Mp5","NA-Mp5", 
                                                          "Healthy-NonInf-Mp","Alcohol-NonInf-Mp","NA-NonInf-Mp",
                                                          "NA-Hep1","NA-Hep2","NA-Hep3","NA-Hep4","NA-Hep5",
                                                          "NA-LSEC1","NA-LSEC2","NA-Chol","NA-Plasma"
                 ))
plots <- lapply(
  X = plots,
  FUN = function(p) p 
  + ggplot2::scale_y_continuous(breaks = c(0,1,2,3))
  + theme(legend.position = 'none') 
  + ggplot2::scale_fill_manual(values = c('black','white','blue',
                                          'black','white','blue',
                                          'black','white','blue',
                                          'black','white','blue',
                                          'black','white','blue',
                                          'black','white',
                                          'black','white','blue',
                                          'black','white',
                                          'black','white','blue',
                                          'black','white','blue',
                                          'blue','blue','blue','blue','blue',
                                          'blue','blue','blue','blue'))
)
CombinePlots(plots = plots, ncol = 1)

# If there are issues, feel free to contact Adam - @atomadam2
