#!/usr/bin/env Rscript

# generate adjacency matrix and topological overlap matrix
# creates signed adj mat and TOM (also tom dissimilarity matrix)
# args: <infile.txt>    White-space delimited file, rows=genes, col=samples 
#                       with header
#       <sft.RData>     Soft thresholding data
#       <out.RData>     Output RData file to store adj, TOM, and dissTOM
# author:   Kevin Chau
# date:     2018 02 02

library(WGCNA)
allowWGCNAThreads()
enableWGCNAThreads()


# argument parsing
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
    stop("Missing arguments.\n", call.=FALSE)
} else if (length(args) == 2) {
    args[3] <- "TOM.RData"
}

# load data
load(args[2])       # sft.RData
datExpr <- t(read.table(args[1], header=TRUE))

# matrix calculations
adj <- adjacency(datExpr, type='signed', power=sft$powerEstimate)
TOM <- TOMsimilarity(adj)
dissTOM <- as.matrix(1-TOM)
save(adj, TOM, dissTOM, file=args[3])
