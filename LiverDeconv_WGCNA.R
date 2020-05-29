
########################################
##                                    ##
##         1. Install Packages        ##
##                                    ##
########################################

install.packages("BiocManager")
BiocManager::install("WGCNA")
BiocManager::install("GO.db")


########################################
##                                    ##
##         2. Load packages           ##
##                                    ##
########################################

# Load the WGCNA package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

setwd("~/XXXXXX/")
getwd()

########################################
##                                    ##
##           3. Load Data             ##
##                                    ##
########################################

# Load geen expression data
Pitt_TPM <- read.delim("XXXX/Pitt_AH_TPM.txt", row.names=1)

# Load clinical Data
traitData <- read.csv("XXXX/Pitt_Clinical_all.csv")
dim(traitData)
names(traitData)

########################################
##                                    ##
##              4. WGCNA              ##
##                                    ##
########################################
# Transpose Data
datExpr0 = as.data.frame(t(Pitt_TPM));

# FIlter bad genes
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Cluster Samples
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Remove outliers - I see 2 possible outlieres (HCV and NAFLD)
# Plot a line to show the cut
abline(h = 150000, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 150000, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# remove columns that hold information we do not need.
allTraits = traitData
#allTraits = traitData[, -c(16:19)];
#allTraits = allTraits[, c(2, 11:36) ];
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.
AllSamples = rownames(datExpr);
traitRows = match(AllSamples, allTraits$SRA);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];
collectGarbage();


# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

save(datExpr, datTraits, file = "PItt_Liver-dataInput.RData")


# Construct network
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# Look at figure for inflection point to determine power
net = blockwiseModules(datExpr, power = 5,                     # Based on graph
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "Pitt_Liver-networkConstruction-auto.RData")

table(net$colors)


#Clinical

lnames = load(file = "C:/Users/thead/Desktop/LiverDeconPaper/Analysis/Pitt_Liver-dataInput.RData");
lnames = load(file = "C:/Users/thead/Desktop/LiverDeconPaper/Analysis/Pitt_Liver-networkConstruction-auto.RData");

# Load data - These are files created above
lnames = load(file = "Pitt_Liver-dataInput.RData");
lnames = load(file = "Pitt_Liver-networkConstruction-auto.RData");

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

pdf("Heatmap_Modules_all.pdf", height = 20, width = 20, useDingbats=FALSE)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = TRUE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


# Define variable weight containing the weight column of datTrait
# REPLACE with Clinical Parameter of interest
weight = as.data.frame(datTraits$AST);
names(weight) = "AST"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

# REPLACE Color module and Clinical Parameter of interest 
module = "yellowgreen"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for AST",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# Output
# REPLACE Color module
names(datExpr)
names(datExpr)[moduleColors=="yellowgreen"]

probes = names(datExpr)

# Create the starting data frame
geneInfo0 = data.frame(substanceBXH = probes,
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
# REPLACE Color module and Clinical Parameter
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.AST));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneInfo_yellowgreen_AST.csv")



# Networks

# NOT WORKING

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")


nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")



# Eigengenes

# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate from the clinical traits
MELD = as.data.frame(datTraits$MELD);
names(MELD) = "MELD"
# Add the Clincal to existing module eigengenes
MET = orderMEs(cbind(MEs, MELD))



AGE = as.data.frame(datTraits$AGE);
names(AGE) = "AGE"
MET = orderMEs(cbind(MET, AGE))

ALB = as.data.frame(datTraits$ALB);
names(ALB) = "ALB"
MET = orderMEs(cbind(MET, ALB))

AST = as.data.frame(datTraits$AST);
names(AST) = "AST"
MET = orderMEs(cbind(MET, AST))

ALT = as.data.frame(datTraits$ALT);
names(ALT) = "ALT"
MET = orderMEs(cbind(MET, ALT))

ALP = as.data.frame(datTraits$ALP);
names(ALP) = "ALP"
MET = orderMEs(cbind(MET, ALP))

Tbil = as.data.frame(datTraits$TBIL);
names(Tbil) = "TBIL"
MET = orderMEs(cbind(MET, Tbil))

PLAT = as.data.frame(datTraits$PLATELET);
names(PLAT) = "PLAT"
MET = orderMEs(cbind(MET, PLAT))

CREAT = as.data.frame(datTraits$CREAT);
names(CREAT) = "CREAT"
MET = orderMEs(cbind(MET, CREAT))

HepOne = as.data.frame(datTraits$NA.Hep1);
names(HepOne) = "HepOne"
MET = orderMEs(cbind(MET, HepOne))

HepTwo = as.data.frame(datTraits$NA.Hep2);
names(HepTwo) = "HepTwo"
MET = orderMEs(cbind(MET, HepTwo))

HepThree = as.data.frame(datTraits$NA.Hep3);
names(HepThree) = "HepThree"
MET = orderMEs(cbind(MET, HepThree))

HepFour = as.data.frame(datTraits$NA.Hep4);
names(HepFour) = "HepFour"
MET = orderMEs(cbind(MET, HepFour))

HepFive = as.data.frame(datTraits$NA.Hep5);
names(HepFive) = "HepFive"
MET = orderMEs(cbind(MET, HepFive))

Chol = as.data.frame(datTraits$NA.Chol);
names(Chol) = "Chol"
MET = orderMEs(cbind(MET, Chol))

PMpOne = as.data.frame(datTraits$Basal.Mp1);
names(PMpOne) = "PMpOne"
MET = orderMEs(cbind(MET, PMpOne))

PMpTwo = as.data.frame(datTraits$Basal.Mp2);
names(PMpTwo) = "PMpTwo"
MET = orderMEs(cbind(MET, PMpTwo))

PMpThree = as.data.frame(datTraits$Basal.Mp3);
names(PMpThree) = "PMpThree"
MET = orderMEs(cbind(MET, PMpThree))

PMpFour = as.data.frame(datTraits$Basal.Mp4);
names(PMpFour) = "PMpFour"
MET = orderMEs(cbind(MET, PMpFour))

PMpFive = as.data.frame(datTraits$Basal.Mp5);
names(PMpFive) = "PMpFive"
MET = orderMEs(cbind(MET, PMpFive))

PNonIn = as.data.frame(datTraits$Basal.NonInf.Mp);
names(PNonIn) = "PNonIn"
MET = orderMEs(cbind(MET, PNonIn))

LMpOne = as.data.frame(datTraits$NA.Mp1);
names(LMpOne) = "LMpOne"
MET = orderMEs(cbind(MET, LMpOne))

LMpThree = as.data.frame(datTraits$NA.Mp3);
names(LMpThree) = "LMpThree"
MET = orderMEs(cbind(MET, LMpThree))

LMpFive = as.data.frame(datTraits$NA.Mp5);
names(LMpFive) = "LMpFive"
MET = orderMEs(cbind(MET, LMpFive))

LNonIn = as.data.frame(datTraits$NA.NonInf.Mp);
names(LNonIn) = "LNonIn"
MET = orderMEs(cbind(MET, LNonIn))




# Plot the relationships among the eigengenes and the trait
sizeGrWindow(10,15);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)

pdf("Heatmap_Pheno_allv2.pdf", height = 20, width = 20, useDingbats=FALSE)
#sizeGrWindow(5,7.5);
#par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
dev.off()


# Plot the dendrogram
sizeGrWindow(10,10);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
