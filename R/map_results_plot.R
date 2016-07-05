setwd("~/Desktop/map_plot_test/")
getwd()

# Script: map_results_plot.r 
# By: Michael G. Harvey 
# Date: 5 July 2016 

# A script to take a table of NCBI Blastn hits and plot them along lines representing 
# chromosomes or scaffolds present in the reference. Hits appear as blue dots above the lines.
# There is also the option to flag certain hits of interest (e.g. outlier loci) by plotting a
# set of red dots below the line representing the reference. The file representing the reference
# should be a tab-delimited text file with three columns (scaffold/chromosome name, RefSeq ID,
# and the size of the scaffold/chromosome in Mbp). The file of hits of interest should just be
# a column listing the names of the loci of interest (identical to how they appear in the Blastn
# results). The Blastn hits should be in the "hit table (.csv)" NCBI format. Some of the code 
# below for parsing the input files may have to be tweaked depending on the vagaries of formatting
# in the various input files.

# Read input files, examples available in ../data/

scaffolds <- read.table("Reference_genome_chromosomes.txt", sep="\t", header=T, row.names=1) # Reference info
hits <- read.table("top_blast_hit_table.txt", sep=",", header=F) # NCBI Blastn hit table
outliers <- read.table("list_of_loci_of_interest.txt", sep="\t", header=F) # List of loci of interest
scaffolds <- scaffolds[order(-scaffolds$Size),] # Order scaffolds from largest to smallest

# Process input, get information for plotting

outliers.hits.rows <- subset(hits, matrix(unlist(strsplit(as.character(hits[,1]), split="\\|")), ncol=2, byrow=T)[,1] %in% outliers$V1) # pull out subset of Blastn results for loci of interest
z.f.chroms <- scaffolds[2:nrow(scaffolds),] # remove first scaffold/chromosome (often "unknown")
sep <- (8/length(z.f.chroms[,1])) # determine distance between reference fragments
scaf.lengths <- 5*(z.f.chroms$Size/max(z.f.chroms$Size)) # determine the size for each fragment
seg.pos <- vector() # vector to record Y positions for axis labels
par(mar=c(5,8,1,2)) # set up plot margins

# Plot

plot(NA, xlim=c(0,5), ylim=c(0,8), xlab="Length (Mbp)", ylab="", axes=F) # empty plot
for(i in 1:nrow(z.f.chroms)) {
	segments(0, 8-(sep*i), scaf.lengths[i], 8-(sep*i), col="gray", lwd=4) # plot fragment
	seg.pos <- c(seg.pos, 8-(sep*i)) # record fragment Y positions for axis labels
	# Add in all hits
	hit.rows <- hits[which(gsub("\\|", "", substr(hits[,2], 18, 31)) == as.character(z.f.chroms$RefSeq[i])),]
	hit.pos <- ((hit.rows[10]-((hit.rows[10]-hit.rows[9])/2))/max(z.f.chroms$Size*1000000))*5
	points(hit.pos[,1], rep(8-(sep*i)+0.05, length(hit.pos[,1])), col="blue", pch=19, cex=0.6)
	# Add in loci of interest
	outliers.rows <- outliers.hits.rows[which(gsub("\\|", "", substr(outliers.hits.rows[,2], 18, 31)) == as.character(z.f.chroms$RefSeq[i])),]
	outliers.pos <- ((outliers.rows[10]-((outliers.rows[10]-outliers.rows[9])/2))/max(z.f.chroms$Size*1000000))*5
	points(outliers.pos[,1], rep(8-(sep*i)-0.05, length(outliers.pos[,1])), col="red", pch=19, cex=0.6)
}
interval <- 20/(max(z.f.chroms$Size)/5) # interval for x axis ticks
axis(1, at=c(0, interval, interval*2, interval*3, interval*4, interval*5, interval*6, interval*7), labels=c(0, 20, 40, 60, 80, 100, 120, 140)) # x axis ticks/labels may need to be tweaked depending on reference
axis(2, at=seg.pos, labels=rownames(z.f.chroms), las=2, col="white") # y axis ticks/labels
