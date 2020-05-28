
########################################
##                                    ##
##   1. Install programs needed       ##
##                                    ##
########################################

install.packages('Seurat')

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
##        3. Load scRNA DATA         ##
##                                    ##
########################################

# Import Human PBMC Data

# This is a function to create the names of allfile paths for all PBMC data from CellRanger
# I am keeping the location and file names as I had it, but they will need to be modified depnding on how you ran everything

SeuratAll <- function(x){
  
  print(x)
  spec_sample <- x
  filename <- (paste("I:/Adam/Unix/scRNA-seq/PBMC_LPS_scRNA/ind_analysis",
                     spec_sample,
                     "outs/filtered_feature_bc_matrix/", sep = "/"))
  print(filename)
  
  testvec <- unlist(strsplit(spec_sample, "_"))
  patientval <- testvec[3]
  stimval <- testvec[4]
  if (patientval=="HC"){
    finalpat <- "Healthy"
  }
  else {
    finalpat <- "Alcohol"
  }
  if (stimval=="B"){
    finalstim <- "Basal"
  }
  else {
    finalstim <- "LPS"
  }
  
  SeuratSample <- Read10X(data.dir = filename)
  
  # Set up control
  Seurat_temp <- CreateSeuratObject(counts = SeuratSample , project = spec_sample)
  
  # store mitochondrial percentage in object meta data
  Seurat_temp <- PercentageFeatureSet(Seurat_temp, pattern = "^MT-", col.name = "percent.mt")
  Seurat_temp$type <- finalpat
  Seurat_temp$stim <- finalstim
  Seurat_temp$CellType <- "none"
  Seurat_temp$Method <- "none"
  Seurat_temp$Experiment <- "none"
  
  # look at numbers - seems similar +/- aggr
  VlnPlot(Seurat_temp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  plot1 <- FeatureScatter(Seurat_temp, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(Seurat_temp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  CombinePlots(plots = list(plot1, plot2))
  
  # run sctransform
  Seurat_temp <- SCTransform(Seurat_temp, vars.to.regress = "percent.mt", verbose = TRUE)
  # FILTER  
  Seurat_temp <- subset(Seurat_temp, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
  
  return(Seurat_temp)
}

# Go through all files and label it specifically
HC_B_01_Seurat <- SeuratAll("pbmc_01_HC_B")
HC_L_02_Seurat <- SeuratAll("pbmc_02_HC_L")
AH_B_03_Seurat <- SeuratAll("pbmc_03_AH_B")
AH_L_04_Seurat <- SeuratAll("pbmc_04_AH_L")

HC_B_05_Seurat <- SeuratAll("pbmc_05_HC_B")
HC_L_06_Seurat <- SeuratAll("pbmc_06_HC_L")
AH_B_07_Seurat <- SeuratAll("pbmc_07_AH_B")
AH_L_08_Seurat <- SeuratAll("pbmc_08_AH_L")

HC_B_09_Seurat <- SeuratAll("pbmc_09_HC_B")
HC_L_10_Seurat <- SeuratAll("pbmc_10_HC_L")
AH_B_11_Seurat <- SeuratAll("pbmc_11_AH_B")
AH_L_12_Seurat <- SeuratAll("pbmc_12_AH_L")

HC_B_13_Seurat <- SeuratAll("pbmc_13_HC_B")
HC_L_14_Seurat <- SeuratAll("pbmc_14_HC_L")
AH_B_15_Seurat <- SeuratAll("pbmc_15_AH_B")
AH_L_16_Seurat <- SeuratAll("pbmc_16_AH_L")

# Import Human Liver Data from GSE115469
# This is pre-analyzed data from GSE. The original fastq don't seem to be available in a usable format

# Get path to the GSE data for healthy liver as downloaded
raw_counts<-read.table(file=paste0("XXXXXX/GSE115469_Data.csv"), 
                       sep=",", header = TRUE, row.names=1, stringsAsFactors=FALSE)

head(raw_counts)

# Create a Seurat Object from .csv
Liverdata <- CreateSeuratObject(counts = raw_counts, min.cells = 3, min.features = 200, project = "Hs_Liver_scRNA")
Liverdata
Liverdata <- subset(Liverdata, subset = nFeature_RNA > 200) # This >200 is a standard cutoff, and important for liver

Liverdata[["percent.mt"]] <- PercentageFeatureSet(Liverdata, pattern = "^MT-") # Label mito content

head(Liverdata@meta.data, 5)

# All original samples names are stored in orig.ident (P1TLH, P2TLH, etc) *There were 5 livers sequenced*

# A few quality control plots - These are nit published but helpful to look at
VlnPlot(Liverdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(Liverdata, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Liverdata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# This is the normalization with SCTransform
Liverdata <- SCTransform(Liverdata, vars.to.regress = "percent.mt", verbose = TRUE)

# FILTER: based on RNA reads and mitochondrial content (This is a fairly strict mito cutoff)
Liverdata <- subset(Liverdata, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)

# Create one big Seurat object of everything
sample_all <- c(HC_B_01_Seurat,HC_L_02_Seurat,AH_B_03_Seurat,AH_L_04_Seurat,
                HC_B_05_Seurat,HC_L_06_Seurat,AH_B_07_Seurat,AH_L_08_Seurat,
                HC_B_09_Seurat,HC_L_10_Seurat,AH_B_11_Seurat,AH_L_12_Seurat,
                HC_B_13_Seurat,HC_L_14_Seurat,AH_B_15_Seurat,AH_L_16_Seurat, Liverdata)

########################################
##                                    ##
##         4. Integrate Data          ##
##                                    ##
########################################

# Use the  reduction = "rpca" for all cells everywhere
Immune.features <- SelectIntegrationFeatures(object.list = sample_all, nfeatures = 3000)
options(future.globals.maxSize = 4800 * 1024^2)

Immune.list <- PrepSCTIntegration(object.list = sample_all, 
                                  anchor.features = Immune.features, 
                                  verbose = TRUE)

Immune.list <- lapply(X = Immune.list, FUN = RunPCA, verbose = TRUE, features = Immune.features)


Immune.anchors <- FindIntegrationAnchors(object.list = Immune.list, 
                                         normalization.method = "SCT",
                                         reduction = "rpca", 
                                         anchor.features = Immune.features, 
                                         verbose = TRUE)
Immune.integrated <- IntegrateData(anchorset = Immune.anchors, normalization.method = "SCT",
                                   verbose = TRUE)

save(Immune.integrated,file="PBMC_Liver.integrated_20200310.Robj")

# If there are issues, feel free to contact Adam - @atomadam2
