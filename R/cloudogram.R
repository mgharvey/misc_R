setwd("/Users/michaelharvey/Documents/Simulations/cloudogram")
getwd()

# Script: cloudogram.r 
# By: Michael G. Harvey 
# Date: 17 April 2014 

# A script to overlay multiple phylogenetic trees to get a feel for the tree depths and topologies.
#
# MAJOR CAVEATS: This script does NOT reorder the tips, so if the taxa do not appear in the same
# order in the trees they will be out of order on the plot. Also, the heights of the trees are 
# nearly on the same scale, but there seems to be some very small variation across trees in the y 
# scale. Really this script should just be used for eyeballing patterns at this point.

require(ape)
require(plyr)

# Read files 

trees <- read.tree("trees.tre") # read in treefile, example available in ../data/

par(mfrow=c(1,1), new=F)

# Get tree depths 

depthlist <- rep(0, 100)
for (i in 1:length(trees)) {
	tree <- trees[[i]] # pull out one tree
	treevector <- data.frame(as.matrix(dist.nodes(tree))[, length(tree$tip)+1], rownames=c(rbind(tree$tip.label), c((length(tree$tip)+1):max(length(dist.nodes(tree)[,1]))))) # vector of node depths
	treedepth = max(treevector[,1]) # tree depth (deepest node depth)
	depthlist[i] <- treedepth # add to list
} 
maxdepth <- max(depthlist)

# Plot trees 

for (i in 1:length(trees)) {
	tree <- trees[[i]] # pull out one tree
	code <- ifelse(is.monophyletic(tree, c(1,2)), ifelse(is.monophyletic(tree, c(3,4)), ifelse(is.monophyletic(tree, c(5,6)), ifelse(is.monophyletic(tree, c(7,8)), T, F), F), F), F) # check for monophyly of all 4 species 
	#code <- ifelse(depthlist[i] > 10.5, F, T)
	col <- ifelse(code == T, "black", "red") # color based on some criterion (in this case monophyly of each species)
	par(fig=c(0.1,0.8,0.1,0.1+((depthlist[i]/maxdepth)*0.8)), mar=c(1,1,0.1,1))
	plot(tree, type="cladogram", direction = "downwards", y.lim=c(0, depthlist[i]), edge.width=0.5, show.tip.label=F, edge.color=col)
	par(new=T) # keep using same plot window
}
