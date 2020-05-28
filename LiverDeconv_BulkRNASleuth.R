
########################################
##                                    ##
##         1. Install Packages        ##
##                                    ##
########################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
BiocManager::install("rhdf5")
install.packages("rlang")
install.packages("Rcpp")
install.packages("digest")

install.packages("devtools")
devtools::install_github("pachterlab/sleuth")


########################################
##                                    ##
##         2. Load packages           ##
##                                    ##
########################################

library("sleuth")
library(ggplot2)
library(psych)
library(statmod)
library(dplyr)
library(readr)
library(plyr)
library(stringr)
library(reshape2)
library(Matrix)

setwd("~/XXXXXX/")
getwd()

########################################
##                                    ##
##     3. Create Sleuth Object        ##
##                                    ##
########################################

# HC - Healthy Control
# EAH - Early AH
# AHL - AH with liver failure
# ExAH - Explant tissue AH
# NAFLD - Non-alcoholic Fatty Liver Disease
# HCV - Hepatitis C 
# AC - Hepatitis C with Cirrhosis

# Get path to all alignment directories
sample_id <- dir(file.path("I:/Adam/Unix/Pitt_AH_Liver/align_b/"))
sample_id

# Get all paths to final results
# This needs to go to a directory that contains all output from kallsito
kal_dirs <- file.path("XXXXX/align_b", sample_id, "abundance.h5")
kal_dirs
kal_dirs <- kal_dirs[-c(22, 43, 54)] # Get rid of 538, 518, and 550, all outliers


# Get metadata, may be able to incorporate more
# The metadata file for all patients

s2c <- read.table(file.path("XXXX/Compiled_Names_v3.txt"), 
                  header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = SRA, AFFECTION_STATUS, Sex_2, Disease_2, Disease_3, Disease_4, Disease_5)
s2c 
#s2c <- s2c[c(1:10),] # For subsets

# Add paths-to-files to metadata
s2c <- dplyr::mutate(s2c, path = kal_dirs)

# Get Gene names from biomart
# This can be fussy
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
# This will generate a matrix of all Enemble IDs (genes and transcripts) with true gene names

# Step 1 - Initialize the sleuth object
so <- sleuth_prep(s2c, 
                  target_mapping = t2g, 
                  aggregation_column = 'ext_gene',        # Required for aggregation of transcripts
                  gene_mode = TRUE,                       # Tells the function that you want genes and not transcripts
                  extra_bootstrap_summary = TRUE,         # Required for DE analysis
                  read_bootstrap_tpm = TRUE)

# Step 2 - PCA
plot_pca(so,   color_by = 'Disease_4')

# What are the drivers of the PCA
plot_loadings(so, pc_input = 1) 

# ALB, FAM83A, HP, FGA, and APOB are driving PC1
plot_bootstrap(so,  target_id = "ALB", units = "scaled_reads_per_base", color_by = "Disease_4")

# What are the drivers of the PCA
plot_loadings(so, pc_input = 2) 

# FAM83A, ALB, MT-CO1, SERPINA1, FGA are driving PC2
plot_bootstrap(so,  target_id = "MT-CO1", units = "scaled_reads_per_base", color_by = "Disease_4")

# Other variables
plot_pca(so, color_by = 'Sex_2') 

# Step 3 - Save Sleuth object
sleuth_save(so, "Pitt_RNA_Sleuth.Robj")

# Step 4 - Save all TPM values
tpm_matrix <- sleuth_to_matrix(so, 'obs_norm', 'tpm')
write.table(tpm_matrix, file='Pitt_AH_TPM.txt')


########################################
##                                    ##
##               Heatmap              ##
##                                    ##
########################################

# This will create the Heatmap
# Requires the entire Sleuth Object

# Load my Sleuth Object
so <- sleuth_load('Pitt_RNA_Sleuth.Robj')
plot_bootstrap(so,  target_id = "ALB", units = "scaled_reads_per_base", color_by = "Disease_4")

### I edited this function ###
# Run the following
# trace("plot_transcript_heatmap",edit=TRUE)
# Make the following change in the script
# cluster_cols=FALSE
##############################
 
# Heatmap
WNT<-c("WNT1","WNT2","WNT2B","WNT3","WNT3A",
       "WNT4","WNT5A","WNT5B","WNT6",
       "WNT7A","WNT7B","WNT8A","WNT8B",
       "WNT9A","WNT9B","WNT10A","WNT10B",
       "WNT11","WNT16"
       )
pdf("HeatmapWNT_Fig3E.pdf", height = 4, width = 12, useDingbats=FALSE)
plot_transcript_heatmap(so, transcripts=WNT, units = "scaled_reads_per_base", trans = "log",
                        scale ="row", cluster_rows = TRUE, #cluster_cols = FALSE,
                        cluster_transcripts = FALSE, annotation_cols = c("Disease_4"),
                        offset = 1, color_high = "blue", color_mid = "lightblue",
                        color_low = "white", x_axis_angle = 50)
dev.off()

# Heatmap
FZD<-c("FZD1","FZD2","FZD3","FZD4","FZD5","FZD6",
       "FZD7","FZD8","FZD9","FZD10"
)
pdf("HeatmapFZD_Fig3E.pdf", height = 2.5, width = 12, useDingbats=FALSE)
plot_transcript_heatmap(so, transcripts=FZD, units = "tpm", trans = "log", 
                        scale ="row", cluster_rows = TRUE, #cluster_cols = FALSE,
                        cluster_transcripts = FALSE, annotation_cols = c("Disease_4"),
                        offset = 1, color_high = "blue", color_mid = "lightblue",
                        color_low = "white", x_axis_angle = 50)
dev.off()


########################################
##                                    ##
##               Boxplots             ##
##                                    ##
########################################

# This can be run without running anything above
# Load TPM File
Pitt_AH_TPM <- read.csv("XXXX/Pitt_AH_TPM.txt", row.names=1, sep="")

# Load metatdata
# Get metadata, may be able to incorporate more
# I:/Adam/Unix/Pitt_AH_Liver/other_files/Compiled_Names_v3.txt

s2c <- read.table(file.path("XXXX/Compiled_Names_v3.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = SRA, Disease_4)
s2c 
row.names(s2c) <- s2c$sample


# SOX9
keep <- as.character(c('SOX9'))
matrix_subbed <- Pitt_AH_TPM[keep,]
matrix_subbed <- na.omit(matrix_subbed)
matrix_subbed2 <- as.data.frame(t(matrix_subbed))

# Join metadata with gene expression
finalData <- cbind(s2c, matrix_subbed2)

# Make figure
# Melt Data
Gene_all <- melt(finalData)
# Change order
levels(Gene_all$Disease_4)
# Specify levels in the order you want
Gene_all$Disease_4 <- factor(Gene_all$Disease_4, levels = c("HC", "EAH", "AHL", "ExAH", "NAFLD", "HCV", "AC"))
p <- ggplot(Gene_all, aes(x=Disease_4, y=value)) + geom_boxplot() +
  theme_bw() + 
  theme(                                 # eliminates background, gridlines, and chart border
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  theme(axis.line = element_line(colour = "black", 
                                 size = 1, linetype = "solid"))+
  ylab("TPM")+
  theme(axis.text.x = element_text(face="bold", color="#000000", size=14),
        axis.text.y = element_text(face="bold", color="#000000", size=14),
        axis.title.x = element_blank(),
        #        axis.title.y = element_blank(),
        legend.text = element_text( size = 14)) +
  theme(axis.title=element_text(size=20), plot.title = element_text(size=24), 
        legend.title=element_text(size=16)) + 
  labs(list(title = "Cell Types", y="TPM")) +
  scale_x_discrete(labels=c("AC" = "HCV_Cirr"))+
  #scale_fill_manual(values = c('grey45','red','grey90','pink')) + 
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.75), dotsize = 0.4)
p

pdf("Boxplot_SOX9_Fig.pdf", height = 4, width = 6, useDingbats=FALSE)
p
dev.off()


# KRT7
keep <- as.character(c('KRT7'))
matrix_subbed <- Pitt_AH_TPM[keep,]
matrix_subbed <- na.omit(matrix_subbed)
matrix_subbed2 <- as.data.frame(t(matrix_subbed))

# Join metadata with gene expression
finalData <- cbind(s2c, matrix_subbed2)

# Make figure
# Melt Data
Gene_all <- melt(finalData)
# Change order
levels(Gene_all$Disease_4)
# Specify levels in the order you want
Gene_all$Disease_4 <- factor(Gene_all$Disease_4, levels = c("HC", "EAH", "AHL", "ExAH", "NAFLD", "HCV", "AC"))
p <- ggplot(Gene_all, aes(x=Disease_4, y=value)) + geom_boxplot() +
  theme_bw() + 
  theme(                                 # eliminates background, gridlines, and chart border
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  theme(axis.line = element_line(colour = "black", 
                                 size = 1, linetype = "solid"))+
  ylab("TPM")+
  theme(axis.text.x = element_text(face="bold", color="#000000", size=14),
        axis.text.y = element_text(face="bold", color="#000000", size=14),
        axis.title.x = element_blank(),
        #        axis.title.y = element_blank(),
        legend.text = element_text( size = 14)) +
  theme(axis.title=element_text(size=20), plot.title = element_text(size=24), 
        legend.title=element_text(size=16)) + 
  labs(list(title = "Cell Types", y="TPM")) +
  scale_x_discrete(labels=c("AC" = "HCV_Cirr"))+
  #scale_fill_manual(values = c('grey45','red','grey90','pink')) + 
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.75), dotsize = 0.4)
p

pdf("Boxplot_KRT7Fig.pdf", height = 4, width = 6, useDingbats=FALSE)
p
dev.off()


# KRT19
keep <- as.character(c('KRT19'))
matrix_subbed <- Pitt_AH_TPM[keep,]
matrix_subbed <- na.omit(matrix_subbed)
matrix_subbed2 <- as.data.frame(t(matrix_subbed))

# Join metadata with gene expression
finalData <- cbind(s2c, matrix_subbed2)

# Make figure
# Melt Data
Gene_all <- melt(finalData)
# Change order
levels(Gene_all$Disease_4)
# Specify levels in the order you want
Gene_all$Disease_4 <- factor(Gene_all$Disease_4, levels = c("HC", "EAH", "AHL", "ExAH", "NAFLD", "HCV", "AC"))
p <- ggplot(Gene_all, aes(x=Disease_4, y=value)) + geom_boxplot() +
  theme_bw() + 
  theme(                                 # eliminates background, gridlines, and chart border
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  theme(axis.line = element_line(colour = "black", 
                                 size = 1, linetype = "solid"))+
  ylab("TPM")+
  theme(axis.text.x = element_text(face="bold", color="#000000", size=14),
        axis.text.y = element_text(face="bold", color="#000000", size=14),
        axis.title.x = element_blank(),
        #        axis.title.y = element_blank(),
        legend.text = element_text( size = 14)) +
  theme(axis.title=element_text(size=20), plot.title = element_text(size=24), 
        legend.title=element_text(size=16)) + 
  labs(list(title = "Cell Types", y="TPM")) +
  scale_x_discrete(labels=c("AC" = "HCV_Cirr"))+
  #scale_fill_manual(values = c('grey45','red','grey90','pink')) + 
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.75), dotsize = 0.4)
p

pdf("Boxplot_KRT19Fig5B.pdf", height = 4, width = 6, useDingbats=FALSE)
p
dev.off()


# WNT2
keep <- as.character(c('WNT2'))
matrix_subbed <- Pitt_AH_TPM[keep,]
matrix_subbed <- na.omit(matrix_subbed)
matrix_subbed2 <- as.data.frame(t(matrix_subbed))

# Join metadata with gene expression
finalData <- cbind(s2c, matrix_subbed2)

# Make figure
# Melt Data
Gene_all <- melt(finalData)
# Change order
levels(Gene_all$Disease_4)
# Specify levels in the order you want
Gene_all$Disease_4 <- factor(Gene_all$Disease_4, levels = c("HC", "EAH", "AHL", "ExAH", "NAFLD", "HCV", "AC"))
p <- ggplot(Gene_all, aes(x=Disease_4, y=value)) + geom_boxplot() +
  theme_bw() + 
  theme(                                 # eliminates background, gridlines, and chart border
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  theme(axis.line = element_line(colour = "black", 
                                 size = 1, linetype = "solid"))+
  ylab("TPM")+
  theme(axis.text.x = element_text(face="bold", color="#000000", size=14),
        axis.text.y = element_text(face="bold", color="#000000", size=14),
        axis.title.x = element_blank(),
        #        axis.title.y = element_blank(),
        legend.text = element_text( size = 14)) +
  theme(axis.title=element_text(size=20), plot.title = element_text(size=24), 
        legend.title=element_text(size=16)) + 
  labs(list(title = "Cell Types", y="TPM")) +
  scale_x_discrete(labels=c("AC" = "HCV_Cirr"))+
  #scale_fill_manual(values = c('grey45','red','grey90','pink')) + 
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.75), dotsize = 0.4)
p

pdf("Boxplot_WNT2Fig5B.pdf", height = 4, width = 6, useDingbats=FALSE)
p
dev.off()


########################################
##                                    ##
##       Differential Expression      ##
##                                    ##
########################################

# For each comparison, code was written separately
# I thought this would be easier so you can find which comparison you are interested in to run
# This section runs just like the beginning. It takes in the output from Kallisto, and makes a custom
# Sleuth object to do the comparison.

########################################
##                                    ##
##         FOR DE : ExAH vs HC        ##
##                                    ##
########################################

# Redoing everything but removing unnecesary samples (have to redo the sleuth object anyway)

# Get path to all alignment directories
sample_id <- dir(file.path("XXXX/align_b/"))
sample_id

# Get all paths to final results
kal_dirs <- file.path("XXXX/align_b", sample_id, "abundance.h5")
kal_dirs
kal_dirs <- kal_dirs[-c(22, 43, 54)] # Get rid of 538, 518, and 550, all outliers


# Get metadata, may be able to incorporate more
# The metadata file for all patients
s2c <- read.table(file.path("XXXX/Compiled_Names_v3.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = SRA, AFFECTION_STATUS, Sex_2, Disease_2, Disease_3, Disease_4, Disease_5)
s2c 

# Add paths-to-files to metadata
s2c <- dplyr::mutate(s2c, path = kal_dirs)

# Remove other diseases
s2c <- s2c[!s2c$Disease_4 == "AHL", ]
s2c <- s2c[!s2c$Disease_4 == "EAH", ]
s2c <- s2c[!s2c$Disease_4 == "HCV", ]
s2c <- s2c[!s2c$Disease_4 == "NAFLD", ]
s2c <- s2c[!s2c$Disease_4 == "AC", ]

# Get Gene names from biomart

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

# This will generate a matrix of all Enemble IDs (genes and transcripts) with true gene names

# Step 1 - Initialize the sleuth object
so <- sleuth_prep(s2c, 
                  target_mapping = t2g, 
                  aggregation_column = 'ext_gene',        # Required for aggregation of transcripts
                  gene_mode = TRUE,                       # Tells the function that you want genes and not transcripts
                  extra_bootstrap_summary = TRUE,         # Required for DE analysis
                  read_bootstrap_tpm = TRUE)

# Step 2 - PCA
plot_pca(so, color_by = 'Disease_4')
plot_pca(so, color_by = 'Sex_2')

# What are the drivers of the PCA
plot_loadings(so, pc_input = 1) 

# ALB, HP, FGA, APOB, and ADH1B are driving PC1
plot_bootstrap(so,  target_id = "ADH1B", units = "scaled_reads_per_base", color_by = "Disease_4")

# What are the drivers of the PCA
plot_loadings(so, pc_input = 2) 

# FAM83A, AC010970.1, SERPINA1, MT-CO1, FGA are driving PC2
plot_bootstrap(so,  target_id = "FAM83A", units = "scaled_reads_per_base", color_by = "Sex_2")

# Step 3 - Fittings
# First fit for condition
so <- sleuth_fit(so, ~Disease_4, 'full')

# Second fit to assume most variation is low
so <- sleuth_fit(so, ~1, 'reduced')

# Step 4 - DE Test
so <- sleuth_lrt(so, 'reduced', 'full')

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
write.table(sleuth_table, file='Pitt_ExAHvsHC.txt')

########################################
##                                    ##
##         FOR DE : AHL vs HC         ##
##                                    ##
########################################

# Redoing everything but removing unnecesary samples (have to redo the sleuth object anyway)

# Get path to all alignment directories
sample_id <- dir(file.path("XXXX/align_b/"))
sample_id

# Get all paths to final results
kal_dirs <- file.path("XXXX/align_b", sample_id, "abundance.h5")
kal_dirs
kal_dirs <- kal_dirs[-c(22, 43, 54)] # Get rid of 538, 518, and 550, all outliers


# Get metadata, may be able to incorporate more
# The metadata file for all patients
s2c <- read.table(file.path("XXXX/Compiled_Names_v3.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = SRA, AFFECTION_STATUS, Sex_2, Disease_2, Disease_3, Disease_4, Disease_5)
s2c 

# Add paths-to-files to metadata
s2c <- dplyr::mutate(s2c, path = kal_dirs)

# Remove other diseases
s2c <- s2c[!s2c$Disease_4 == "ExAH", ]
s2c <- s2c[!s2c$Disease_4 == "EAH", ]
s2c <- s2c[!s2c$Disease_4 == "HCV", ]
s2c <- s2c[!s2c$Disease_4 == "NAFLD", ]
s2c <- s2c[!s2c$Disease_4 == "AC", ]

# Get Gene names from biomart

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

# This will generate a matrix of all Enemble IDs (genes and transcripts) with true gene names

# Step 1 - Initialize the sleuth object
so <- sleuth_prep(s2c, 
                  target_mapping = t2g, 
                  aggregation_column = 'ext_gene',        # Required for aggregation of transcripts
                  gene_mode = TRUE,                       # Tells the function that you want genes and not transcripts
                  extra_bootstrap_summary = TRUE,         # Required for DE analysis
                  read_bootstrap_tpm = TRUE)

# Step 2 - PCA
plot_pca(so, color_by = 'Disease_4')
plot_pca(so, color_by = 'Sex_2')

# What are the drivers of the PCA
plot_loadings(so, pc_input = 1) 

# ALB, HP, FGA, APOB, and SERPINA1 are driving PC1
plot_bootstrap(so,  target_id = "ALB", units = "scaled_reads_per_base", color_by = "Disease_4")

# What are the drivers of the PCA
plot_loadings(so, pc_input = 2) 

# FAM83A, SERPINA1, MT-CO1, FGA, APOB are driving PC2
plot_bootstrap(so,  target_id = "FAM83A", units = "scaled_reads_per_base", color_by = "Sex_2")

# Step 3 - Fittings
# First fit for condition
so <- sleuth_fit(so, ~Disease_4, 'full')

# Second fit to assume most variation is low
so <- sleuth_fit(so, ~1, 'reduced')

# Step 4 - DE Test
so <- sleuth_lrt(so, 'reduced', 'full')

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
write.table(sleuth_table, file='Pitt_AHLvsHC.txt')

########################################
##                                    ##
##         FOR DE : EAH vs HC         ##
##                                    ##
########################################

# Redoing everything but removing unnecesary samples (have to redo the sleuth object anyway)

# Get path to all alignment directories
sample_id <- dir(file.path("XXXX/align_b/"))
sample_id

# Get all paths to final results
kal_dirs <- file.path("XXXX/align_b", sample_id, "abundance.h5")
kal_dirs
kal_dirs <- kal_dirs[-c(22, 43, 54)] # Get rid of 538, 518, and 550, all outliers


# Get metadata, may be able to incorporate more
# The metadata file for all patients
s2c <- read.table(file.path("XXXX/Compiled_Names_v3.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = SRA, AFFECTION_STATUS, Sex_2, Disease_2, Disease_3, Disease_4, Disease_5)
s2c 

# Add paths-to-files to metadata
s2c <- dplyr::mutate(s2c, path = kal_dirs)

# Remove other diseases
s2c <- s2c[!s2c$Disease_4 == "ExAH", ]
s2c <- s2c[!s2c$Disease_4 == "AHL", ]
s2c <- s2c[!s2c$Disease_4 == "HCV", ]
s2c <- s2c[!s2c$Disease_4 == "NAFLD", ]
s2c <- s2c[!s2c$Disease_4 == "AC", ]

# Get Gene names from biomart

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

# This will generate a matrix of all Enemble IDs (genes and transcripts) with true gene names

# Step 1 - Initialize the sleuth object
so <- sleuth_prep(s2c, 
                  target_mapping = t2g, 
                  aggregation_column = 'ext_gene',        # Required for aggregation of transcripts
                  gene_mode = TRUE,                       # Tells the function that you want genes and not transcripts
                  extra_bootstrap_summary = TRUE,         # Required for DE analysis
                  read_bootstrap_tpm = TRUE)

# Step 2 - PCA
plot_pca(so, color_by = 'Disease_4')
plot_pca(so, color_by = 'Sex_2')

# What are the drivers of the PCA
plot_loadings(so, pc_input = 1) 

# ALB, HP, FGA, APOB, and SERPINA1 are driving PC1
plot_bootstrap(so,  target_id = "ALB", units = "scaled_reads_per_base", color_by = "Disease_4")

# What are the drivers of the PCA
plot_loadings(so, pc_input = 2) 

# FAM83A, HP, ALB, FGA, MT-CO1 are driving PC2
plot_bootstrap(so,  target_id = "FAM83A", units = "scaled_reads_per_base", color_by = "Sex_2")

# Step 3 - Fittings
# First fit for condition
so <- sleuth_fit(so, ~Disease_4, 'full')

# Second fit to assume most variation is low
so <- sleuth_fit(so, ~1, 'reduced')

# Step 4 - DE Test
so <- sleuth_lrt(so, 'reduced', 'full')

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
write.table(sleuth_table, file='Pitt_EAHvsHC.txt')


########################################
##                                    ##
##         FOR DE : HCV vs HC         ##
##                                    ##
########################################

# Redoing everything but removing unnecesary samples (have to redo the sleuth object anyway)

# Get path to all alignment directories
sample_id <- dir(file.path("XXXX/align_b/"))
sample_id

# Get all paths to final results
kal_dirs <- file.path("XXXX/align_b", sample_id, "abundance.h5")
kal_dirs
kal_dirs <- kal_dirs[-c(22, 43, 54)] # Get rid of 538, 518, and 550, all outliers


# Get metadata, may be able to incorporate more
# The metadata file for all patients
s2c <- read.table(file.path("XXXX/Compiled_Names_v3.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = SRA, AFFECTION_STATUS, Sex_2, Disease_2, Disease_3, Disease_4, Disease_5)
s2c 

# Add paths-to-files to metadata
s2c <- dplyr::mutate(s2c, path = kal_dirs)

# Remove other diseases
s2c <- s2c[!s2c$Disease_4 == "ExAH", ]
s2c <- s2c[!s2c$Disease_4 == "AHL", ]
s2c <- s2c[!s2c$Disease_4 == "EAH", ]
s2c <- s2c[!s2c$Disease_4 == "NAFLD", ]
s2c <- s2c[!s2c$Disease_4 == "AC", ]

# Get Gene names from biomart

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

# This will generate a matrix of all Enemble IDs (genes and transcripts) with true gene names

# Step 1 - Initialize the sleuth object
so <- sleuth_prep(s2c, 
                  target_mapping = t2g, 
                  aggregation_column = 'ext_gene',        # Required for aggregation of transcripts
                  gene_mode = TRUE,                       # Tells the function that you want genes and not transcripts
                  extra_bootstrap_summary = TRUE,         # Required for DE analysis
                  read_bootstrap_tpm = TRUE)

# Step 2 - PCA
plot_pca(so, color_by = 'Disease_4')
plot_pca(so, color_by = 'Sex_2')

# What are the drivers of the PCA
plot_loadings(so, pc_input = 1) 

# ALB, HP, FAM83B, FGA, APOB are driving PC1
plot_bootstrap(so,  target_id = "ALB", units = "scaled_reads_per_base", color_by = "Disease_4")

# What are the drivers of the PCA
plot_loadings(so, pc_input = 2) 

# FAM83A, ALB, HP, MT-CO1, B2M are driving PC2
plot_bootstrap(so,  target_id = "B2M", units = "scaled_reads_per_base", color_by = "Disease_4")

# Step 3 - Fittings
# First fit for condition
so <- sleuth_fit(so, ~Disease_4, 'full')

# Second fit to assume most variation is low
so <- sleuth_fit(so, ~1, 'reduced')

# Step 4 - DE Test
so <- sleuth_lrt(so, 'reduced', 'full')

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
write.table(sleuth_table, file='Pitt_HCVvsHC.txt')


########################################
##                                    ##
##        FOR DE : NAFLD vs HC        ##
##                                    ##
########################################

# Redoing everything but removing unnecesary samples (have to redo the sleuth object anyway)

# Get path to all alignment directories
sample_id <- dir(file.path("XXXX/align_b/"))
sample_id

# Get all paths to final results
kal_dirs <- file.path("XXXX/align_b", sample_id, "abundance.h5")
kal_dirs
kal_dirs <- kal_dirs[-c(22, 43, 54)] # Get rid of 538, 518, and 550, all outliers


# Get metadata, may be able to incorporate more
# The metadata file for all patients
s2c <- read.table(file.path("XXXX/Compiled_Names_v3.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = SRA, AFFECTION_STATUS, Sex_2, Disease_2, Disease_3, Disease_4, Disease_5)
s2c 

# Add paths-to-files to metadata
s2c <- dplyr::mutate(s2c, path = kal_dirs)

# Remove other diseases
s2c <- s2c[!s2c$Disease_4 == "ExAH", ]
s2c <- s2c[!s2c$Disease_4 == "AHL", ]
s2c <- s2c[!s2c$Disease_4 == "EAH", ]
s2c <- s2c[!s2c$Disease_4 == "HCV", ]
s2c <- s2c[!s2c$Disease_4 == "NAFLD", ]

# Get Gene names from biomart

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

# This will generate a matrix of all Enemble IDs (genes and transcripts) with true gene names

# Step 1 - Initialize the sleuth object
so <- sleuth_prep(s2c, 
                  target_mapping = t2g, 
                  aggregation_column = 'ext_gene',        # Required for aggregation of transcripts
                  gene_mode = TRUE,                       # Tells the function that you want genes and not transcripts
                  extra_bootstrap_summary = TRUE,         # Required for DE analysis
                  read_bootstrap_tpm = TRUE)

# Step 2 - PCA
plot_pca(so, color_by = 'Disease_4')
plot_pca(so, color_by = 'Sex_2')

# What are the drivers of the PCA
plot_loadings(so, pc_input = 1) 

# ALB, HP, FGA, APOB, FAM83A  are driving PC1
plot_bootstrap(so,  target_id = "ADH1B", units = "scaled_reads_per_base", color_by = "Disease_4")

# What are the drivers of the PCA
plot_loadings(so, pc_input = 2) 

# FAM83A, ALB, HP, FGA, MT-CO1 are driving PC2
plot_bootstrap(so,  target_id = "ALB", units = "scaled_reads_per_base", color_by = "Disease_4")

# Step 3 - Fittings
# First fit for condition
so <- sleuth_fit(so, ~Disease_4, 'full')

# Second fit to assume most variation is low
so <- sleuth_fit(so, ~1, 'reduced')

# Step 4 - DE Test
so <- sleuth_lrt(so, 'reduced', 'full')

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
write.table(sleuth_table, file='Pitt_HCV_CirrvsHC.txt')


########################################
##                                    ##
##        FOR DE : NAFLD vs HC        ##
##                                    ##
########################################

# Redoing everything but removing unnecesary samples (have to redo the sleuth object anyway)

# Get path to all alignment directories
sample_id <- dir(file.path("XXXX/align_b/"))
sample_id

# Get all paths to final results
kal_dirs <- file.path("XXXX/align_b", sample_id, "abundance.h5")
kal_dirs
kal_dirs <- kal_dirs[-c(22, 43, 54)] # Get rid of 538, 518, and 550, all outliers


# Get metadata, may be able to incorporate more
# The metadata file for all patients
s2c <- read.table(file.path("XXXX/Compiled_Names_v3.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = SRA, AFFECTION_STATUS, Sex_2, Disease_2, Disease_3, Disease_4, Disease_5)
s2c 

# Add paths-to-files to metadata
s2c <- dplyr::mutate(s2c, path = kal_dirs)

# Remove other diseases
s2c <- s2c[!s2c$Disease_4 == "ExAH", ]
s2c <- s2c[!s2c$Disease_4 == "AHL", ]
s2c <- s2c[!s2c$Disease_4 == "EAH", ]
s2c <- s2c[!s2c$Disease_4 == "HCV", ]
s2c <- s2c[!s2c$Disease_4 == "AC", ]

# Get Gene names from biomart

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

# This will generate a matrix of all Enemble IDs (genes and transcripts) with true gene names

# Step 1 - Initialize the sleuth object
so <- sleuth_prep(s2c, 
                  target_mapping = t2g, 
                  aggregation_column = 'ext_gene',        # Required for aggregation of transcripts
                  gene_mode = TRUE,                       # Tells the function that you want genes and not transcripts
                  extra_bootstrap_summary = TRUE,         # Required for DE analysis
                  read_bootstrap_tpm = TRUE)

# Step 2 - PCA
plot_pca(so, color_by = 'Disease_4')
plot_pca(so, color_by = 'Sex_2')

# What are the drivers of the PCA
plot_loadings(so, pc_input = 1) 

# ALB, HP, APOB, FGA, ADH1B  are driving PC1
plot_bootstrap(so,  target_id = "ADH1B", units = "scaled_reads_per_base", color_by = "Disease_4")

# What are the drivers of the PCA
plot_loadings(so, pc_input = 2) 

# FAM83A, MT-CO1, HP, ALB, APOB are driving PC2
plot_bootstrap(so,  target_id = "ALB", units = "scaled_reads_per_base", color_by = "Disease_4")

# Step 3 - Fittings
# First fit for condition
so <- sleuth_fit(so, ~Disease_4, 'full')

# Second fit to assume most variation is low
so <- sleuth_fit(so, ~1, 'reduced')

# Step 4 - DE Test
so <- sleuth_lrt(so, 'reduced', 'full')

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
write.table(sleuth_table, file='Pitt_NAFLDvsHC.txt')

