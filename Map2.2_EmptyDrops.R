##### NOTE ####
#
# This script is specifically designed for datasets where an excessive number of cells were 
# called using cellranger's pipeline. DO NOT USE on datasets where cell ranger identified 
# fewer than 20,000 cells!!!
#
##############

require(methods)
require(Matrix)
require("DropletUtils")
require(dplyr)
require(Seurat)
script_dir = "/cluster/home/tandrews/scripts/LiverMap2.0"

args <- as.numeric(as.character(commandArgs(trailingOnly=TRUE)))
# folder of raw data
# name of project
# max cells
# min cells
folder <- args[1]
name <- args[2];
if (length(args) < 4) {
	n_min_cells <- 100
	if (length(args) < 3){
		n_max_cells <- 20000
	} else {
		n_max_cells <- args[3];
	}
} else {
	n_max_cells <- args[3];
	n_min_cells <- args[4];
}

FDR=0.01

# Read the raw matrix
mydata <- Read10X(data.dir = paste(folder, sep="/"))
# Remove rows & columns that are completely zero
mydata <- mydata[,Matrix::colSums(mydata) > 0]
mydata <- mydata[Matrix::rowSums(mydata) > 0,]
print(paste("Stats out:", dim(outdata), "input dimensions of", name), collapse=" ")

# Rank droplet barcodes by total UMIs
br.out <- barcodeRanks(mydata)

# Get total UMIs that correspond to the min & max thresholds.
plausible_cell_threshold <- max(br.out$total[br.out$rank > n_max_cells]);
mandatory_cell_threshold <- max(br.out$total[br.out$rank < n_min_cells]);

# Knee and Inflection point plot
#plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
#o <- order(br.out$rank)
#lines(br.out$rank[o], br.out$fitted[o], col="red")

#abline(h=br.out$knee, col="dodgerblue", lty=2)
#abline(h=br.out$inflection, col="forestgreen", lty=2)
#legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
#    legend=c("knee", "inflection"))


set.seed(100)
# Run EmptyDrops
e.out <- emptyDrops(mydata, lower=max(plausible_cell_threshold, min(br.out$knee, br.out$inflection)), niters=100000, ignore=3, retain=mandatory_cell_threshold)

# Clean up results
e.out <- e.out[!is.na(e.out$PValue),] # Remove NAs (too few UMIs to estimate PValue)
e.out$q.value <- e.out$PValue;
e.out$q.value[e.out$Limited] <- 0; # Set those with p < 1/n_simulations to p = 0 to ensure they are kept after multiple testing correction
e.out$q.value <- p.adjust(e.out$q.value, method="fdr"); # apply FDR multiple testing correction.

# Get number of detected genes/droplet
mydata <- mydata[,match(rownames(e.out),colnames(mydata))]
e.out$ngenes <- Matrix::colSums(mydata>0);

# Significance threshold
is.cell <- e.out$q.value <= FDR
sum(is.cell, na.rm=TRUE)

# Plot of selected genes
png(paste(name, "emptyDrops.png", sep="_"), width=6, height=6, units="in", res=150)
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
    xlab="Total UMI count", ylab="-Log Probability")
dev.off()

# Subset the matrix to the selected droplets and save it to a file.
outdata <- mydata[,is.cell]
print(paste("Stats out:", dim(outdata), "output dimensions of", name), collapse=" ")
saveRDS(outdata, paste(name, "emptyDrops_table.rds", sep="_"));
