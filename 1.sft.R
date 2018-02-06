#!/usr/bin/env Rscript

# Calculate soft thresholding statistics and generate relevant plots
# args: <infile.txt>    White-space delimited file, rows=genes, col=samples 
#                       with header
#       <sftimg.png>    Image output file
#       <outfile.RData> Output RData file
# author:   Kevin Chau
# date:     2018 02 02

library(WGCNA)
allowWGCNAThreads()
enableWGCNAThreads()


# argument parsing
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
    stop("Missing arguments.\n", call.=FALSE)
} else if (length(args) == 1) {
    args[2] <- "softThresholding.RData"
    args[3] <- "sft.RData"
}

# load data
datExpr <- t(read.table(args[1], header=TRUE))

# soft-thresholding calculations
powers <- c(c(1:15), seq(from=16, to=30, by=2))
sft <- pickSoftThreshold(datExpr, powerVector=powers, verbose=5)

# output as image file
png(args[2], width=1920, height=1080) 
par(mfrow=c(1, 2), mar=c(10,10,8,2))
cex1 <- 4
cex2 <- 4
plot(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     xlab="", 
     ylab="", 
     type="n", main=paste("Scale Independence"), 
     cex.lab=cex1, cex.axis=cex2, cex.main=cex1,
     mgp=c(0,3,0)
     )
title(xlab="Soft Threshold (power)", 
      ylab="Scale Free Topology Model Fit, signed R^2", 
      mgp=c(7, 1, 0), cex.lab=cex1
      )
text(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     labels=powers, cex=cex2, col="red"
     )

# plot line at height=0.8 (scale-free topology network should level out 
# after this point)
abline(h=0.8, col="red")

# mean connectivity plot
plot(sft$fitIndices[,1], sft$fitIndices[,5], 
     xlab="", ylab="", type="n", 
     main = paste("Mean connectivity"),
     cex.lab=cex1, cex.axis=cex2, cex.main=cex1,
     mgp=c(0,3,0)
     ) 
title(xlab="Soft Threshold (power)", 
      ylab="Mean Connectivity", 
      mgp=c(7, 1, 0), cex.lab=cex1
      )
text(sft$fitIndices[,1], sft$fitIndices[,5], 
     labels=powers, cex=cex2, col="red"
     )
dev.off()

