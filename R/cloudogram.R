setwd("/Users/michaelharvey/Documents/Simulations/cloudogram")
getwd()

##################################
### cloudogram.r ################
### By Michael G. Harvey #########
### 16 April 2014 ################
##################################

library(ape)
library(plyr)

### Read files ###

trees <- read.tree("trees.tre") # read in treefile
codes <- c(rep(1, 80), rep(2, 20)) # make a vector of some code

par(mfrow=c(1,1), new=F)

### Get tree depths ###

depthlist <- rep(0, 100)
for (i in 1:length(trees)) {
	tree <- trees[[i]] # pull out one tree
	treevector <- data.frame(as.matrix(dist.nodes(tree))[, length(tree$tip)+1], rownames=c(rbind(tree$tip.label), c((length(tree$tip)+1):max(length(dist.nodes(tree)[,1])))))
	treedepth = max(treevector[,1])
	depthlist[i] <- treedepth
} 
maxdepth <- max(depthlist)

### Plot trees ###

for (i in 1:length(trees)) {
	tree <- trees[[i]] # pull out one tree
	col <- ifelse(codes[i] == 1, "red", "black") # color based on some criterion
	par(mai=c(0.5,0.5,(maxdepth+0.5)-depthlist[i],0.5))
	plot(tree, type="cladogram", direction = "downwards", y.lim=c(0, depthlist[i]), edge.width=0.5, show.tip.label=F, edge.color=col)
	par(new=T)
}
