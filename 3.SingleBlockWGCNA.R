#!/usr/bin/env Rscript

# generate WGCNA network using single-block calculations
# network output is multiple *.MODULE files named by color along with 
# relevant plotsinto a single directory
# args: <TOM.RData>     RData file with dissTOM object 
#       <data.txt>      Original expression matrix
#       <network.RData> Where to save network RData
#       <outdir>        Custom directory name?
#                       Caution: This script will overwrite any directories
#                                with this provided name or, if not provided,
#                                "MODULES"
# NOTE: User should modify parameters as necessary to fit their data
# author:   Kevin Chau
# date:     2018 02 02

library(WGCNA)

# argument parsing
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
    stop("Must supply at least the TOM and data!\n", call.=FALSE)
} else if (length(args) == 2) {
    args[3] = "network.RData"
    args[4] = "MODULES"
} else if (length(args) == 3) {
    args[4] = "MODULES"
}

# recommended parameters to modify
minModuleSize <- 50         # minimum module size for module identification
MEDissThresh <- 0.15        # cutoff for module eigengene clustering

# load data
load(args[1])       # TOM.RData should include at least dissTOM
datExpr <- t(read.table(args[2], header=TRUE))

# network construction
tree <- hclust(as.dist(dissTOM), method='average')
dynamicMods <- cutreeDynamic(dendro=tree, distM=dissTOM,
                             deepSplit=2, minClusterSize=minModuleSize
                             )
dynamicColors <- labels2colors(dynamicMods)
png("prelimDendrogram.png", width=1920, height=1080)
plotDendroAndColors(tree, dynamicColors, "", 
                    dendroLabels=FALSE, hang=0.03,
                    addGuide=TRUE, guideHang=0.05,
                    mar=c(3,10,7,3), mgp=c(6,0,0),
                    main="Preliminary Dendrogram and Module Colors",
                    cex.main=6, cex.axis=4, cex.lab=5
                    )
dev.off()

# module eigengene clustering and merging
MEList <- moduleEigengenes(datExpr, colors=dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method='average')
png("moduleEigengeneClustering.png", width=1920, height=1080)
par(mar=c(3,11,7,3), mgp=c(6,0,0))
plot(METree, main="Clustring of Module Eigengenes",
     xlab='', ylab="Height", sub='',
     cex.main=6, cex.axis=4, cex.lab=6, cex=3
     )
abline(h=MEDissThresh, col='red')
dev.off()
merge <- mergeCloseModules(datExpr, dynamicColors, 
                           cutHeight=MEDissThresh,
                           verbose=3
                           )
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

# final dendrogram and modules
png("mergedDendro.png", width=1920, height=1080)
par(mgp=c(6,0,0))
plotDendroAndColors(tree, cbind(dynamicColors, mergedColors),
                    c("DynamicTreeCut", "MergedDynamic"),
                    dendroLabels=F, hang=0.03,
                    addGuide=T, guideHang=0.05,
                    mar=c(3,25,7,3),
                    cex.main=6, cex.axis=4, cex.lab=5,
                    cex.colorLabels=4
                    )
dev.off()
moduleColors <- mergedColors
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder) - 1
MEs <- mergedMEs;
save(moduleLabels, moduleColors, tree,
     MEs, mergedMEs, file=args[3]
     )

# output gene lists into modules
dir.create(args[4], showWarnings=F)
assignedGenes <- colnames(datExpr)
names(assignedGenes) <- moduleColors
for (module in unique(names(assignedGenes))) {
    filename <- sprintf("%s/moduleMembers_%s.MODULE", args[4], module)
    moduleGenes <- assignedGenes[names(assignedGenes)==module]
    write.table(moduleGenes, file=filename, sep='\n', quote=F, 
                row.names=F, col.names=F
                )
}
