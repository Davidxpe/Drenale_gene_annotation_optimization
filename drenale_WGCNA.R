#library
library(tximeta)
library(tximport)
library(DESeq2)
library(WGCNA);            
library(tidyr);         
library(dynamicTreeCut);    
library(flashClust);       
library(lattice);            
library(survival);          
library(Hmisc);            
library(dplyr)
################################################################################

#db
dir <- file.path(here::here("Drenale_counts_TPM/quant"))
samples = read.delim("Drenale_counts_TPM/metadata", header = TRUE, sep = "\t")
row.names(samples) <- samples$run
files <- file.path(dir, paste0(samples$sample, "_quant"), "quant.sf")
txi <- tximport(files, type="salmon", tx2gene=tx2gene, ignoreTxVersion=TRUE)
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ 1)

##inputdata
counts <- counts(ddsTxi, normalized = FALSE, replaced = FALSE)
counts <- as.data.frame(counts)
counts$drm2_tail <- NULL
#TPM
tpm<-txi$abundance
colnames(tpm) <- colnames(counts)
tpm <- as.data.frame(tpm)
tpm$drm2_tail <- NULL
#VST
vst <- vst(ddsTxi)
vst <- assay(vst)
vst <- as.data.frame(vst)
vst$drm2_tail <- NULL
#edgeR
counts <- txi$counts
dge <- DGEList(counts)
dge <- calcNormFactors(dge)
dge = as.data.frame(dge)

#WGCNA
drenalenorm <- read.delim("VST")
traits <- read.delim("traits_tissue")

#variable_date
row.names(drenalenorm) <- sub("\\.t[1-9]$", "", row.names(drenalenorm))
SubGeneNames=row.names(drenalenorm);
typeof(SubGeneNames);
SubGeneNames[0:10];
drenale <- as.data.frame(t(drenalenorm))

#define_power
powers = c(c(1:10), seq(from =10, to=30, by=1)) 
sft = pickSoftThreshold(drenale,
                        dataIsExpr = TRUE,
                        powerVector = powers,
                        corFnc = cor,                             
                        corOptions = list(use = 'p'),             
                        verbose = 5,
                        moreNetworkConcepts = TRUE,              
                        networkType = "signed");  
#Scale-free topology
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, Signed R^2",
     type="n", main = paste("Scale independence for\ D. renale tissues\ (Matrix E)"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red");  
abline(h=0.93,col="blue");  
abline(h=0.95,col="brown");  

#generate_matrix
adjacency = adjacency(drenale,type = "signed",   
                      power = softPower);
TOM = TOMsimilarityFromExpr(adjacency,                         
                            TOMType = "signed", 
                            power = softPower);
dissTOM = 1-TOM

#genetree

geneTree = hclust(as.dist(dissTOM), method="average")
plot(geneTree, xlab="Gene clusters", 
     main="Gene clustering on TOM−based dissimilarity for D. renale Tissues (Matrix E)", 
     cex.main = 1,
     cex.axis = 1,
     cex=0.3)
#set minimum module
minModuleSize = 20;
dynamicMods = cutreeDynamic(dendro= geneTree, 
                            distM = dissTOM,        
                            deepSplit=2,           
                            pamRespectsDendro= FALSE,      
                            minClusterSize = minModuleSize)
#identify by colors
dynamicColors = labels2colors(dynamicMods)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.05,
                    addGuide = TRUE, guideHang = 0.1,
                    cex.main = 1,
                    main = "Gene dendrogram and module colors for D. renale tissues \n Matrix E")

#merge similar modules(PCA)
MEList= moduleEigengenes(drenale, colors= dynamicColors,
                         excludeGrey = TRUE,
                         softPower = softPower)
MEs = MEList$eigengenes

#eigenges dissimilarity
MEDiss = 1-cor(MEs)          
#cluster it
METree = hclust(as.dist(MEDiss),method= "average")
MEDissThres = 0.10
plot(METree, main= "Clustering of module eigengenes for D.renale Tissues (Matrix E)", 
     cex.main = 1.2,
     xlab= "",      
     ylab = "",
     cex = 0.5,     
     sub= "");   
##merging function
merge = mergeCloseModules(drenale, dynamicColors, 
                          cutHeight= MEDissThres, 
                          verbose =3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
#Reasign the new variables to the old one 
moduleColors = mergedColors
#Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));         
moduleLabels = match(mergedColors, colorOrder)-1     
MEs = mergedMEs

#Write modules no merged in the path location
length(dynamicColors)
module_colors = setdiff(unique(dynamicColors), "grey")  
length(module_colors)
for (color in module_colors){
  module = SubGeneNames[which(dynamicColors==color)]
  write.table(module, paste("module_",color, ".txt",sep=""), 
              sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}

#Write MEs merged in the path location
length(mergedColors)
module_colors = setdiff(unique(mergedColors), "grey")     
length(module_colors)
#module_colors
for (color in module_colors){
  module=SubGeneNames[which(mergedColors==color)]
  write.table(module, paste("merged_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}

#correlate to traits
nGenes = ncol(drenale)
nSamples = nrow(drenale)
#recalculate MEs with color labels.   
MEs0 = moduleEigengenes(drenale, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

#link the external values to the correlation matrix
moduleTraitCor2 = cor(MEs, traits, use= "p")
#calculates Student asymptotic p-value for given correlations.
moduleTraitPvalue2 = corPvalueStudent(moduleTraitCor2, nSamples)   
textMatrix2= paste(signif(moduleTraitCor2, 2), " (", signif(moduleTraitPvalue2, 1), ")", sep= "")
dim(textMatrix2) = dim(moduleTraitCor2)

#heatmap
par(mar= c(3.5, 10, 2, 1))  
labeledHeatmap(Matrix= moduleTraitCor2,
               xLabels= names(traits),
               yLabels = names(MEs),
               yLabelsPosition = "left",
               yColorWidth = 0.1,
               cex.lab.y = 0.8,
               cex.lab.x = 1,
               ySymbols= names(MEs),
               colors= blueWhiteRed(50), 
               textMatrix= (textMatrix2),
               setStdMargins= FALSE,
               cex.text= 0.6,
               main= paste("Module-trait"))