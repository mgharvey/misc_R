setwd("/Users/michaelharvey/Documents/Simulations/cloudogram")
getwd()

##################################
### cloudogram.r ################
### By Michael G. Harvey #########
### 17 April 2014 ################
##################################

library(ape)
library(plyr)

### Read files ###

trees <- read.tree("trees.tre") # read in treefile

par(mfrow=c(1,1), new=F)

### Get tree depths ###

depthlist <- rep(0, 100)
for (i in 1:length(trees)) {
	tree <- trees[[i]] # pull out one tree
	treevector <- data.frame(as.matrix(dist.nodes(tree))[, length(tree$tip)+1], rownames=c(rbind(tree$tip.label), c((length(tree$tip)+1):max(length(dist.nodes(tree)[,1]))))) # vector of node depths
	treedepth = max(treevector[,1]) # tree depth (deepest node depth)
	depthlist[i] <- treedepth # add to list
} 
maxdepth <- max(depthlist)

### Plot trees ###

for (i in 1:length(trees)) {
	tree <- trees[[i]] # pull out one tree
	code <- ifelse(is.monophyletic(tree, c(1,2)), ifelse(is.monophyletic(tree, c(3,4)), ifelse(is.monophyletic(tree, c(5,6)), ifelse(is.monophyletic(tree, c(7,8)), T, F), F), F), F) # check for monophyly of all 4 species 
	col <- ifelse(code == T, "black", "red") # color based on some criterion (in this case monophyly of each species)
	par(mai=c(0.5,0.5,(maxdepth+0.5)-depthlist[i],0.5))
	plot(tree, type="cladogram", direction = "downwards", y.lim=c(0, depthlist[i]), edge.width=0.5, show.tip.label=F, edge.color=col)
	par(new=T) # keep using same plot window
}
