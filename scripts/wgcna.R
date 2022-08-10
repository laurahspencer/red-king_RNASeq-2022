# I used this script to run on Sedna.
# Because many of the code requires examining before moving on to the next step, I ran this interactively. To do so I used this code: 
# ` srun --pty --mem=180GB /bin/bash`  # Initiate an interactive job, specify the max amount of memory on the higher-numbered nodes (nodes29-36 can handle up to 180GB of memory). 
# If needed I can also request a himem node (which you can easily do like 1300GB if not more), via `-p himem -c 24` to specify the himem node with all 24 cores
# `module load R` #load R module
# `R` #start R interactively 

# Install and load packages, and proceed! 
# NOTE: before doing the above I created a new directory on Sedna to save packages, /home/lspencer/R/library,
# then, after logging on to R I followed the directions for "Installing your own packages" here: https://hpcf.umbc.edu/cpu/other-software/how-to-run-r-programs-on-taki/ via
# `.libPaths("/home/lspencer/R/library") 
# I then logged off of of R via `quit()` and added the path to .bashrc via: `nano .bashrc` and added `export R_LIBS_USER=path/to/my_R_Packages` to the end of the file 

# To run this using a slurm scheduler, I wrote the slurm script "wgcna.sh" which simply does this to execute this R code: 
# `module load R`
# `R CMD BATCH /home/lspencer/2022-redking-OA/scripts/wgcna.R /home/lspencer/2022-redking-OA/wgcna/wgcna-output.txt`


# CODE STARTS HERE

# Add all required libraries that are installed with install.packages() here
list.of.packages <- c("tidyverse", "RColorBrewer", "ggpmisc", "ggpubr", "ape", "pheatmap", "plotly", "clipr")
# Add all libraries that are installed using BiocManager here
bioconductor.packages <- c("WGCNA", "DESeq2")

# This commented out code need only be run once per machine (installs packages). Don't re-do it, it can take a while. 
# # Install BiocManager if needed
# if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# 
# # Get names of all required packages that aren't installed
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
# new.bioc.packages <- bioconductor.packages[!(bioconductor.packages %in% installed.packages()[, "Package"])]
# # Install all new packages
# if(length(new.packages)) install.packages(new.packages)
# if(length(new.bioc.packages)) BiocManager::install(new.bioc.packages)

# Load all required libraries
all.packages <- c(list.of.packages, bioconductor.packages)
lapply(all.packages, FUN = function(X) {
  do.call("require", list(X))
})


# Load in datasets
# load data
load(file = "/home/lspencer/2022-redking-OA/wgcna/counts.trans") #counts.trans
#load(file = "/home/lspencer/2022-redking-OA/wgcna/per-gene-summary-stats") #genes.stats
load(file = "/home/lspencer/2022-redking-OA/wgcna/vsd.pH") #vsd.pH

# Read in sample metadata
load(file="/home/lspencer/2022-redking-OA/wgcna/sample.info") #sample.info

# Load annotation objects 
load(file = "/home/lspencer/2022-redking-OA/wgcna/counts.annot.Pcamt")
load(file = "/home/lspencer/2022-redking-OA/wgcna/P.camt.blast.GO")

# load in DEG lists 
load(file = "/home/lspencer/2022-redking-OA/wgcna/diffex.AM")
load(file = "/home/lspencer/2022-redking-OA/wgcna/diffex.AL")
load(file = "/home/lspencer/2022-redking-OA/wgcna/diffex.ML")

# Load in DESeq2 results lists to pull background list of genes 
load(file = "/home/lspencer/2022-redking-OA/wgcna/res.pco2.AM") #res.pco2.AM
load(file = "/home/lspencer/2022-redking-OA/wgcna/res.pco2.AL") #res.pco2.AL
load(file = "/home/lspencer/2022-redking-OA/wgcna/res.pco2.ML") #res.pco2.ML

#Set CPU cores for parallel-related functions
cpucores <- 20
require(parallel)
options("mc.cores"=cpucores)

#Set CPU cores for doParallel-related functions
require(doParallel)
cores <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cores)
enableWGCNAThreads()

# Use WGCNA's `goodSamplesGenes` to remove any genes with very low variability or missing data  
datExpr0 <- counts.trans[,-1:-3]
#datExpr0 <- counts.trans[,4:1000] #for testing purposes use subset of data


gsg = goodSamplesGenes(datExpr0, verbose = 3);
sum(gsg$goodGenes==FALSE)  #were any removed?  Nope. 
ncol(gsg$allOK)

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

# Identify outlier samples 
# Decision March 16th - don't throw out Tank4_Crab1, it's not an outlier on the PCAs (in DESeq2 notebook) 

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)

pdf(file = "/home/lspencer/2022-redking-OA/wgcna/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = 250, col = "red");
dev.off()


# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 250, minSize = 10)  #for bowtie2 aligned
#clust = cutreeStatic(sampleTree, cutHeight = 125, minSize = 10)  #for star aligned
table(clust)

# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# remove columns that hold information we do not need.
(allTraits = counts.trans[, 1:3] %>% rownames_to_column("Sample") %>%
    mutate(pH=Treatment) %>%
    mutate_at("pH", function(x) {
      case_when(
        x == "Ambient" ~ 8.0,
        x == "Moderate" ~ 7.8,
        x == "Low" ~ 7.6)}) %>%
    mutate(H_exp=(10^(-1*pH)) %>% formatC(., format = "e", digits = 2)) %>%
    mutate(pco2=Treatment) %>%
    mutate_at("pco2", function(x) {
      case_when(
        x == "Ambient" ~ 371,
        x == "Moderate" ~ 704,
        x == "Low" ~ 1424)}))

dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the traits.
 
Samples = rownames(datExpr);
traitRows = match(Samples, allTraits$Sample);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];

collectGarbage();


# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(sapply(datTraits,as.numeric), signed = F, colors = colorRampPalette(brewer.pal(8, "Set2"))(14));

# Plot the sample dendrogram and the colors underneath.
pdf(file = "/home/lspencer/2022-redking-OA/wgcna/sampleTree.pdf", width = 12, height = 9);
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram with treatment and tank ids")
dev.off()



options(stringsAsFactors = FALSE);


# Construct a weighted gene network entails the choice of the soft thresholding power β to which co-expression similarity is raised to calculate adjacency.
# The authors of [1] have proposed to choose the soft thresholding power based on the criterion of approximate scale-free topology."

# Choose a set of soft-thresholding powers
powers = c(c(1:25))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")

# Plot the results:
pdf(file = "/home/lspencer/2022-redking-OA/wgcna/soft-thresholding-power.pdf", width = 12, height = 9);
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
dev.off()

# I will choose the power 15, which is the lowest power for which the scale-free topology fit index reaches 0.90.


# I kept getting this error about memory
# > adjacency = adjacency(datExpr, power = softPower, type="signed", corFnc="bicor");
# Error: cannot allocate vector of size 83.4 Gb
# So I required 180GB of memory for this job, and it worked

# Calculate the adjacencies, using the soft thresholding power 15:
softPower = 15; # for bowtie2 aligned to RKC genome 
adjacency = adjacency(datExpr, power = softPower, type="signed", corFnc="bicor", );
save(adjacency, file = "/home/lspencer/2022-redking-OA/wgcna/adjacency")

#To minimize effects of noise and spurious associations, we transform the adjacency into Topological Overlap Matrix, and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency); #this takes a very long time
dissTOM = 1-TOM

save(dissTOM, file = "/home/lspencer/2022-redking-OA/wgcna/dissTom-bowtie-rkc")
#load(file = "/home/lspencer/2022-redking-OA/wgcna/dissTom-bowtie-rkc") #load if needed

#We now use hierarchical clustering to produce a hierarchical clustering tree (dendrogram) of genes. Note that we use the function hclust that provides a much faster hierarchical clustering routine than the standard hclust function.

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
save(geneTree, file = "/home/lspencer/2022-redking-OA/wgcna/geneTree")

# Plot the resulting clustering tree (dendrogram)

pdf(file = "/home/lspencer/2022-redking-OA/wgcna/gene-clustering.pdf", width = 12, height = 9);
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()

# We like large modules, so we set the minimum module size relatively high:
#minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = 30);

dynamicMods.50 = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = 50);

dynamicMods.75 = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = 75);


table(dynamicMods)
save(dynamicMods, file = "/home/lspencer/2022-redking-OA/wgcna/dynamicMods")
save(dynamicMods.50, file = "/home/lspencer/2022-redking-OA/wgcna/dynamicMods.50")
save(dynamicMods.75, file = "/home/lspencer/2022-redking-OA/wgcna/dynamicMods.75")


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram and colors underneath
pdf(file = "/home/lspencer/2022-redking-OA/wgcna/gene-clustering-modules.pdf", width = 12, height = 9);
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene clustering on TOM-based dissimilarity\nDendrogram and module colors")
dev.off()


# The Dynamic Tree Cut may identify modules whose expression profiles are very similar. It may be prudent to merge
# such modules since their genes are highly co-expressed. To quantify co-expression similarity of entire modules, we
# calculate their eigengenes and cluster them on their correlation.
# 
# We choose a height cut of 0.25, corresponding to correlation of 0.75, to merge. 

# Calculate eigengenes
#MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

# Plot the result
pdf(file = "/home/lspencer/2022-redking-OA/wgcna/module-tree.pdf", width = 12, height = 9);
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()

# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

pdf(file = "/home/lspencer/2022-redking-OA/wgcna/gene-clustering-modules-merged.pdf", width = 12, height = 9);
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene clustering on TOM-based dissimilarity\nDendrogram and module colors")
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

# Save module colors and labels for use in subsequent parts
#save(MEs, moduleLabels, moduleColors, geneTree, file = "/home/lspencer/2022-redking-OA/wgcna/networkConstruction-stepByStep.RData")

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# Calculate correlation between eigengenes and treatment 
#moduleTraitCor = cor(MEs, datTraits[,c("pH"), drop=F], use = "p");
moduleTraitCor = cor(MEs, datTraits[,c("pco2"), drop=F], use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# Correlation among module eigengenes 
#cor(MEs0, use="complete.obs") %>% corrplot::corrplot(tl.cex=.75, diag = F)

# library(PerformanceAnalytics)
# chart.Correlation(MEs0[,], histogram=F, pch=19)

# Will display correlations and their p-values
# textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
#                            signif(moduleTraitPvalue, 1), ")", sep = "");
textMatrix =  paste(signif(moduleTraitCor, 2), " (",
                    signif(moduleTraitPvalue, 2), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

pdf(file = "/home/lspencer/2022-redking-OA/wgcna/module-pCO2-corr.pdf", width = 4.5, height = 6);
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               #               xLabels = names(datTraits[,c("pH"), drop=F]),
               xLabels = names(datTraits[,c("pco2"), drop=F]),
               yLabels = names(MEs) %>% str_remove(pattern = "ME"),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               cex.main = 0.85,
               cex.lab.y = 0.75,
               cex.lab.x= 0.85,
               zlim = c(-1,1),
               main = paste("Module-pCO2 relationships"))
dev.off()


# Define variable containing the pH column
# pH = as.data.frame(datTraits$pH);
# names(pH) = "pH"
pco2 = as.data.frame(datTraits$pco2);
names(pco2) = "pco2"

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

#geneTraitSignificance = as.data.frame(cor(datExpr, pH, use = "p"));
geneTraitSignificance = as.data.frame(cor(datExpr, pco2, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

# names(geneTraitSignificance) = paste("GS.", names(pH), sep="");
# names(GSPvalue) = paste("p.GS.", names(pH), sep="");
names(geneTraitSignificance) = paste("GS.", names(pco2), sep="");
names(GSPvalue) = paste("p.GS.", names(pco2), sep="");

# Need to have every gene represented in the annotation file - here I merge the best blast (has all genes) & uniprot blast (some genes) results to have all genes, and to also be able to use either blast result. 
annot = counts.annot.Pcamt %>% 
  filter(geneID %in% names(datExpr))
dim(annot)
names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$geneID)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.

# Create the starting data frame
geneInfo0 = data.frame(id = probes,
                       SPID = annot$SPID[probes2annot],
                       protein_names = annot$protein_names[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance for pH
#modOrder = order(-abs(cor(MEs, pH, use = "p")));
modOrder = order(-abs(cor(MEs, pco2, use = "p")));

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
#geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.pH));
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.pco2));

geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "/home/lspencer/2022-redking-OA/wgcna/geneInfo.csv")
save(geneInfo, file="/home/lspencer/2022-redking-OA/wgcna/geneInfo")


# departure from WGCNA tutorials 

## Generate line plots of genes for each module of interest

# Generate per-gene summary stats 
# first extract the normalized counts
gene.counts.norm <- assay(vsd.pH) %>% as.matrix() %>% as.data.frame() %>% rownames_to_column(var = "gene") %>% 
  pivot_longer(-gene, values_to = "counts.vsd", names_to = "Sample") %>%
  left_join(sample.info[,c("Sample", "Treatment")]) %>% 
  mutate_at(c("gene", "Sample"), factor) %>%
  mutate(Treatment=factor(Treatment, levels=c("Ambient", "Moderate", "Low"))) %>% 
  filter(Sample %in% rownames(datExpr)) # <-- important, keep only samples that were retained in WGCNA analysis 

# calculate mean, sd, cv and annotate 
genes.stats <- gene.counts.norm %>% group_by(Treatment, gene) %>% 
  summarize(mean=mean(counts.vsd), median=median(counts.vsd), sd=sd(counts.vsd), cv=sd(counts.vsd)/mean(counts.vsd)) %>% 
  mutate(DEG=if_else(gene %in% c(rownames(diffex.AM),rownames(diffex.AL), rownames(diffex.ML)), "DEG", "Non-DEG")) %>%
  mutate(change = case_when(
    gene %in% c(rownames(subset(res.pco2.AM, padj < 0.05 & log2FoldChange < 0)), rownames(subset(res.pco2.AL, padj < 0.05 & log2FoldChange < 0))) ~ "Up-regulated",
    gene %in% c(rownames(subset(res.pco2.AM, padj < 0.05 & log2FoldChange > 0)), rownames(subset(res.pco2.AL, padj < 0.05 & log2FoldChange > 0))) ~ "Down-regulated")) %>% 
  mutate(change=replace_na(change,"Non-DEG")) %>% 
  mutate(DEG=factor(DEG, levels=c("Non-DEG", "DEG"))) %>% 
  mutate(change=factor(change, levels=c("Non-DEG", "Down-regulated", "Up-regulated"))) %>%
  ungroup() %>% 
  left_join(counts.annot.Pcamt %>% dplyr::select(geneID, SPID, protein_names, evalue), by=c("gene"="geneID"))

# Generate line plots of significant modules showing only genes significantly correlated to module 
# List of modules (their color IDs) that are significantly associated with pH 
#modules <- moduleTraitPvalue %>% as.data.frame() %>% filter(pH<0.05) %>% rownames() %>% str_remove(pattern = "ME")
modules <- cbind(moduleTraitCor %>% as.data.frame(), moduleTraitPvalue %>% as.data.frame()) %>% dplyr::rename(corr=1, pvalue=2) %>% filter(pvalue<0.05) %>% arrange(corr) %>%
  rownames() %>% str_remove(pattern = "ME")
modules.p <- paste("p.MM.", modules, sep = "")
modules.cor <- noquote(paste("MM.", modules, sep = ""))

save(modules, file="/home/lspencer/2022-redking-OA/wgcna/modules")
save(modules.p, file="/home/lspencer/2022-redking-OA/wgcna/modules.p")
save(modules.cor, file="/home/lspencer/2022-redking-OA/wgcna/modules.cor")

# create empty list to store summary statistics for each modules' sign. genes 
module.stats <- vector("list", length(modules))
names(module.stats) <- modules

# create empty list to store line plots for each modules' sign. genes 
module.lineplots <- vector("list", length(modules))
names(module.lineplots) <- modules

for (i in 1:length(modules)) {
  a <- geneInfo %>%
    filter(moduleColor==modules[i]) %>% 
    #  filter(p.GS.pH<0.05) %>% filter(!!as.symbol(modules.p[i])<0.05) %>% 
    filter(p.GS.pco2<0.05) %>% filter(!!as.symbol(modules.p[i])<0.05) %>% 
    #  dplyr::select(id, moduleColor, GS.pH, p.GS.pH, !!as.symbol(modules.p[i]), !!as.symbol(modules.cor[i])) %>% 
    dplyr::select(id, moduleColor, GS.pco2, p.GS.pco2, !!as.symbol(modules.p[i]), !!as.symbol(modules.cor[i])) %>% 
    left_join(genes.stats, by = c("id"="gene")) %>% 
    #group_by(Treatment, moduleColor) %>% mutate(mean_z = scale(mean)) %>% ungroup() %>% #use this line of code to z-transform 
    mutate(pH=Treatment) %>%
    mutate_at("pH", function(x) {
      case_when(
        x == "Ambient" ~ 8.0,
        x == "Moderate" ~ 7.8,
        x == "Low" ~ 7.6)}) %>%
    mutate(H_exp=(10^(-1*pH)) %>% formatC(., format = "e", digits = 2)) %>%
    mutate(pco2=Treatment) %>%
    mutate_at("pco2", function(x) {
      case_when(
        x == "Ambient" ~ 371,
        x == "Moderate" ~ 704,
        x == "Low" ~ 1424)})
  
  #formula <- y ~ x
  #colors <- c("brown","pink","cyan2", "red","cyan4","green","darkgreen")
  module.lineplots[[i]] <-  
    #ggplot(a, aes(x=pH, y=median)) + 
    ggplot(a, aes(x=pco2, y=median)) + 
    geom_line(aes(group=id), size=.01, color=modules[i]) +
    theme_minimal(base_size = 8) + facet_wrap(~moduleColor + change) + 
    geom_smooth(method="loess",se = FALSE, color="blue", size=1) + #  method=lm, formula = formula,
    #scale_x_continuous(breaks=c(7.6, 7.8, 8.0)) + 
    theme(plot.title = element_text(size=8), axis.text.x = element_text(hjust = 0.2), 
          axis.title.x = element_text(size=6)) +
    #xlab("pH Treatment") #+
    xlab("pCO2 concentration") #+
  # stat_poly_eq(aes(label = ..eq.label..),
  #              label.x.npc = 0.5, label.y.npc = 0.99,
  #              formula = formula, parse = TRUE, size = 3)
  
  module.stats[[i]] <- a %>% dplyr::select(-Treatment, -mean, -sd, -cv, -pH) %>% distinct()
}
#print(module.lineplots)
save(module.lineplots, file = "/home/lspencer/2022-redking-OA/wgcna/module.lineplots")

#Save module lineplots pdf
fig <- ggarrange(plotlist=module.lineplots,
                 ncol=2, nrow=3,
                 common.legend = TRUE, legend="top")
pdf(file = "/home/lspencer/2022-redking-OA/wgcna/module-lineplots.pdf", width = 8.5, height = 11, onefile = TRUE);
fig
dev.off()

#Summary of modules significantly associated with pCO2 
moduleTraitCor %>% nrow()
modules %>% length()
cbind(moduleTraitCor %>% as.data.frame(), moduleTraitPvalue %>% as.data.frame()) %>% dplyr::rename(corr=1, pvalue=2) %>% filter(pvalue<0.05) %>% arrange(corr)

# For pH-associated modules, look at correlation between module membership and gene significane for each gene 
par(mfrow = c(4,2));

for (i in 1:length(modules)) {
  module = modules[i]
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  
  #sizeGrWindow(7, 7);
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for pH",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
}


# Look at overlap among WGCNA and DESeq2 results 

## What percentage of DEGs are captured by pH-associated WGCNA modules? Particularly interested in genes I will use for enrichment analysis, so I will filter genes for: 
# - Only genes whose abundances are significantly correlated with pH  (i.e. Gene Significance (GS) p-value, p.GS, <0.05)
# - NO ---> Only genes whose module membership (MM) is significant (p.MM < 0.05). 


# First how many DEGs do I have total? 
n.DEGs <- genes.stats %>% filter(DEG=="DEG") %>% dplyr::select(gene) %>% distinct() %>% nrow()

# merge module association info with DEG info 
geneInfo %>%
  left_join(genes.stats %>% dplyr::select(gene, DEG) %>% distinct(), by = c("id"="gene")) %>% #join with table indicating DEG
  mutate(moduleColor = as.factor(moduleColor)) %>%
  #group_by(moduleColor, DEG) %>% summarise(n=n(), sign=sum(p.GS.pH<0.05), not.sign=sum(p.GS.pH>0.05)) %>% 
  group_by(moduleColor, DEG) %>% summarise(n=n(), sign=sum(p.GS.pco2<0.05), not.sign=sum(p.GS.pco2>0.05)) %>% 
  filter(moduleColor %in% modules) %>% 
  filter(DEG=="DEG") %>% 
  ungroup() %>% dplyr::select(n, sign) %>% colSums()/n.DEGs*100 

# ANSWER: 
# 81.4% of the DEGs were assigned to a pH-associated module. 
# 73.9% of the DEGs were both assigned to a pH-associated module, AND those genes are sign. associated with pH (p.GS.pH<0.05)


## What percentage of DEGs are also identified as correlated with pH via WGCNA (and assigned to any module)

geneInfo %>%
  left_join(genes.stats %>% dplyr::select(gene, DEG) %>% distinct(), by = c("id"="gene")) %>% #join with table indicating DEG
  mutate(moduleColor = as.factor(moduleColor)) %>%
  #  group_by(moduleColor, DEG) %>% summarise(n=n(), sign=sum(p.GS.pH<0.05), not.sign=sum(p.GS.pH>0.05)) %>% 
  group_by(moduleColor, DEG) %>% summarise(n=n(), sign=sum(p.GS.pco2<0.05), not.sign=sum(p.GS.pco2>0.05)) %>% 
  filter(DEG=="DEG") %>% 
  ungroup() %>% dplyr::select(n, sign, not.sign) %>% colSums()/n.DEGs*100 
# 85% of DEGs are deemed to be significantly correlated with pH


## How many genes are correlated with pH (any module?)

geneInfo %>%
  mutate(moduleColor = as.factor(moduleColor)) %>%
  #  group_by(moduleColor) %>% summarise(n=n(), sign=sum(p.GS.pH<0.05), not.sign=sum(p.GS.pH>0.05)) %>% 
  group_by(moduleColor) %>% summarise(n=n(), sign=sum(p.GS.pco2<0.05), not.sign=sum(p.GS.pco2>0.05)) %>% 
  dplyr::select(moduleColor, n, sign, not.sign) %>% 
  mutate(perc.sig=sign/n) %>% 
  dplyr::select(sign, not.sign) %>% colSums()

# # Venn diagram of DEGs and genes in WGCNA pH-associated modules 
# 
# display_venn <- function(x, ...){
#   library(VennDiagram)
#   grid.newpage()
#   venn_object <- venn.diagram(x, filename = NULL, ...)
#   grid.draw(venn_object)
# }
# 
# module.genes <- vector("list", length(modules)+1)
# names(module.genes) <- c(modules, "all")
# for (i in 1:length(modules)) {
#   module.genes[[i]] <- geneInfo %>%
#     # filter for genes assigned to pH-associated modules 
#     filter(moduleColor==modules[i]) %>% 
#     ## further filter for genes that are significantly associated with pH and significantly assigned to their respective module
#     #filter(p.GS.pH<0.05) %>% 
#     #filter(!!as.symbol(modules.p[i])<0.05) %>%
#     ## alternative approach to just select "hub genes"   
#     #    filter(abs(GS.pH)>0.2) %>% filter(!!as.symbol(modules.cor[i])>0.8) %>% 
#     rownames()
# }
# module.genes[[length(module.genes)]] <- Reduce(c, module.genes)
# 
# x <- list(
#   AvsM = row.names(diffex.AM), 
#   AvsL = row.names(diffex.AL), 
#   MvsL = row.names(diffex.ML),
#   WGCNA = module.genes[["all"]]
# )
# 
# display_venn(
#   x,
#   category.names = c("DEGs\nAmbient vs.\nModerate" , "DEGs\nAmbient vs.\nLow" , "DEGs\nModerate vs.\nLow", "WGCNA\npH-associated\nmodules"),
#   # Circles
#   lwd = 2,
#   lty = 'blank',
#   fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
#   # Numbers
#   cex = .9,
#   fontface = "italic",
#   # Set names
#   cat.cex = .9,
#   cat.fontface = "bold",
#   cat.default.pos = "outer")#,
# #cat.dist = c(0.06, 0.06, 0.06))
# 
# # EFFECT OF VARIOUS FILTERING SETTINGS: 
# # When I include all genes assigned to pH-associated modules, the % of DEGs included in WGCNA modules is 80.1%
# # When I filter for only genes that are significantly associated with pH (p.GS.pH<0.05), that drops to 73.3%
# # When I filter for genEs that are sign. associated with pH (p.GS.pH<0.05) AND module membership is sign. (p.MM<0.05) that drops to 72.2%
# # Because I am using all genes assigned to a pH-associated module in the GO Enrichment analysis and IPA, I shall use the least filtered set. 
# 
# 
# 
# Now look for overlap among DEGs and each module 
# 
# 
# par(mfrow = c(4,2));
# for (i in 1:length(modules)) {
#   y <- list(
#     AvsM = row.names(diffex.AM), 
#     AvsL = row.names(diffex.AL), 
#     MvsL = row.names(diffex.ML),
#     WGCNA = module.genes[[i]]
#   )
#   
#   display_venn(
#     y,
#     category.names = c("DEGs\nAmbient vs.\nModerate" , "DEGs\nAmbient vs.\nLow" , "DEGs\nModerate vs.\nLow", 
#                        paste("WGCNA module: ", modules[i], "\npH-associated\nmodules", sep="")),
#     # Circles
#     lwd = 2,
#     lty = 'blank',
#     fill = c("#999999", "#E69F00", "#56B4E9", modules[i]),
#     # Numbers
#     cex = .9,
#     fontface = "italic",
#     # Set names
#     cat.cex = .9,
#     cat.fontface = "bold",
#     cat.default.pos = "outer")
# }


# How many genes assigned to each pH-associated module?

geneInfo %>%  
  filter(moduleColor %in% modules) %>%
  dplyr::count(moduleColor)

# How many total genes assigned to any pH-associated module? 
geneInfo %>%  
  filter(moduleColor %in% modules) %>%
  dplyr::count(moduleColor) #%>% dplyr::select(n) %>% colSums()

# How many total genes assigned to any pH-associated module AND have high Gene Significance (significantly associated with pH)? 
geneInfo %>%  
  filter(moduleColor %in% modules) %>%
  #  filter(p.GS.pH<0.05) %>%
  filter(p.GS.pco2<0.05) %>%
  dplyr::count(moduleColor) # %>% dplyr::select(n) %>% colSums()


# Generate heatmaps of pH-associated modules

# create empty list to store line plots for each modules' sign. genes 
counts.modules <- vector("list", length(modules))
names(counts.modules) <- modules

#filter(abs(GS.pH)>0.2) %>% filter(!!as.symbol(modules.cor[i])>0.8)

for (i in 1:length(modules)) {
  counts.modules[[i]] <- 
    (assay(vsd.pH))[(geneInfo %>% 
                       filter(moduleColor==modules[i]) %>% 
                       filter(p.GS.pco2<0.05))$id,] %>% 
    #                       filter(p.GS.pH<0.05))$id,] %>% 
    as.data.frame() %>% add_column(module=modules[i])
}

s <- sample.info %>% 
  mutate(Treatment=factor(Treatment, levels=c("Ambient", "Moderate", "Low"))) %>% 
  arrange(factor(Treatment, levels=c("Ambient", "Moderate", "Low"))) %>% 
  remove_rownames() %>% column_to_rownames("Sample")

# define order of modules from most sign. positive correlation <--> most sign. negative correlation 
modules.sorted <- moduleTraitCor %>% as.data.frame() %>% 
  filter(rownames(.) %in% 
           #           rownames(moduleTraitPvalue %>% as.data.frame() %>% filter(pH<0.05))) %>% 
           rownames(moduleTraitPvalue %>% as.data.frame() %>% filter(pco2<0.05))) %>% 
  #  arrange(pH) %>% rownames() %>% str_remove(pattern = "ME")
  arrange(pco2) %>% rownames() %>% str_remove(pattern = "ME")

counts.modules.matrix<- bind_rows(counts.modules[modules.sorted]) %>% 
  dplyr::select(module, rownames(s)) %>%
  mutate(module=factor(module, levels=modules.sorted))

pdf(file = "/home/lspencer/2022-redking-OA/wgcna/wgcna-modules-heatmap.pdf", width = 12, height = 9);
pheatmap(counts.modules.matrix[,colnames(counts.modules.matrix) !="module"], 
         cluster_cols = FALSE, cluster_rows=F, annotation_names_col = FALSE, 
         show_rownames=FALSE, na.rm=TRUE, scale="row", 
         main = "Expression by pH-associated module\n(significant genes only)", fontsize = 8, 
         # annotation_colors = list(Treatment=c(Ambient="#2c7bb6", Moderate="#fdae61", Low="#d7191c"),
         #                          module=c(cyan="cyan", pink="pink", brown="brown", 
         #                                   lightcyan1="lightcyan1", green="green", darkgreen="darkgreen")),
         # annotation_colors = list(Treatment=c(Ambient="#2c7bb6", Moderate="#fdae61", Low="#d7191c"),
         #                           module=c(magenta="magenta", orangered4="orangered4", darkgreen="darkgreen",
         #                                    plum1="plum1", salmon="salmon",
         #                                    darkgrey="darkgrey")),
         annotation_col=s[,"Treatment", drop=F],
         annotation_row = counts.modules.matrix[,"module", drop=FALSE])
dev.off()


## Prep for GO_MWU analysis 

#### Extract module membership values for each gene for each pH-associated module 

# Generate a gene annotation file 
geneInfo.mwu <- geneInfo %>% left_join(
  (P.camt.blast.GO)[,c("gene_id","gene_ontology_i_ds")], by = c("id"="gene_id")) %>%  # could also opt to  %>% filter(evalue < 1e-10)
  filter(!is.na(gene_ontology_i_ds)) #filter for genes that are annotated with GO information 
write_delim(geneInfo.mwu[,c("id", "gene_ontology_i_ds")], "/home/lspencer/2022-redking-OA/wgcna/WGCNA-genes_for-GOMWU.tab",delim = "\t")

# Write out separate tables for each module that contain two columns: 1) gene id, and 2) module membership, i.e. significance. For each module, that column contains MM for genes assigned to that module, and 0 for genes not assigned to that module. 
geneInfo.mwu.mods <- vector("list", length(modules))
names(geneInfo.mwu.mods) <- modules

for (i in 1:length(modules)) {
  geneInfo.mwu.mods[[i]] <- geneInfo.mwu %>%
    mutate(sign = ifelse(moduleColor == modules[[i]],!!as.symbol(modules.cor[i]), 0)) %>%
    dplyr::select(id,sign) %>% remove_rownames()
  write.csv(geneInfo.mwu.mods[[i]], paste("/home/lspencer/2022-redking-OA/wgcna/WGCNA-module_",modules[[i]],".csv", sep = ""), quote = F,row.names = F)
}


# # Helpful R code from Ariana: 
# iris %>%
#   group_by(Species) %>%
#   summarise(across(starts_with("Sepal"), list(mean = mean, sd = sd), .names = "{.col}.{.fn}"))
