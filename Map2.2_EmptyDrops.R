##### NOTE ####
#
# Arguments:
# folder of the raw data matrix from 10X
# name of project to use for output files
# maximum plausible number of cells (default: 20,000)
# minimum number of cells (default: 100)
# FDR threshold (default: 0.01  (1%))
##############



require(methods)
require(Matrix)
require("DropletUtils")
require(dplyr)
require(Seurat)
script_dir = "/cluster/home/tandrews/scripts/LiverMap2.0"
source(paste(script_dir, "My_R_Scripts.R", sep="/"));

args <- as.character(commandArgs(trailingOnly=TRUE))
folder <- args[1]
name <- args[2];
print(args);
#maximum number of plausible cells
if (length(args) < 3) { 
	MAX_CELLS <- 20000 
} else {
	MAX_CELLS <- args[3]
}
#mimum number of cells
if (length(args) < 4) {
	MIN_CELLS <- 100
} else {
	MIN_CELLS <- args[4]
}
#FDR for calling a cell
if (length(args) < 5) {
	FDR=0.01
} else {
	FDR <- args[5]
}
# completely exclude droplets below this threshold
if (length(args) < 6) {
	TRIM=10
} else {
	TRIM <- args[6]
}


# Read the raw matrix
if (!grepl("raw", folder)) {
	folder <- paste(folder, "raw_gene_bc_matrices/GRCh38", sep="/");
}
mydata <- Read10X(data.dir = paste(folder, sep="/"))
# Remove rows & columns that are completely zero
mydata <- mydata[,Matrix::colSums(mydata) > TRIM]
mydata <- mydata[Matrix::rowSums(mydata) > 0,]
print(paste("Stats out:", dim(mydata), "input",c("gene","cells"), "of", name))

# Rank droplet barcodes by total UMIs
br.out <- barcodeRanks(mydata)

# Get total UMIs that correspond to the min & max thresholds.
plausible_cell_threshold <- max(br.out$total[br.out$rank > MAX_CELLS]);
mandatory_cell_threshold <- max(br.out$total[br.out$rank < MIN_CELLS]);

# Knee and Inflection point plot
#plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
#o <- order(br.out$rank)
#lines(br.out$rank[o], br.out$fitted[o], col="red")

#abline(h=br.out$knee, col="dodgerblue", lty=2)
#abline(h=br.out$inflection, col="forestgreen", lty=2)
#legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
#    legend=c("knee", "inflection"))

# My calling threshold
n_umi_sorted <- br.out$total[order(br.out$total, decreasing = T)]
rank_sorted <- 1:length(n_umi_sorted);

slope <- diff(log(n_umi_sorted))/diff(log(rank_sorted))
smooth_slope <- smooth.spline(slope, spar=0.5)
inflection <- which(smooth_slope$y == min(smooth_slope$y[MIN_CELLS:MAX_CELLS]))
my_inflection <- n_umi_sorted[inflection];




if (br.out$knee > MIN_CELLS*10 || br.out$inflection > MIN_CELLS*10) {
	is_empty_threshold = max(plausible_cell_threshold, 
				br.out$total[br.out$knee], 
				br.out$total[br.out$inflection], my_inflection)
} else {
	is_empty_threshold = plausible_cell_threshold
}

set.seed(100)
# Run EmptyDrops
e.out <- emptyDrops(mydata, lower=is_empty_threshold, niters=100000, ignore=TRIM, retain=mandatory_cell_threshold)

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
    xlab="Total UMI count", ylab="-Log Probability", pch=18)
dev.off()

# Subset the matrix to the selected droplets and save it to a file.
outdata <- mydata[,is.cell]
print(paste("Stats out:", dim(outdata), "output",c("gene","cells"), "of", name))
saveRDS(outdata, paste(name, "emptyDrops_table.rds", sep="_"));


my_background_correction <- function(sample_mat, is.cell) {
	my_background <- Matrix::rowSums(sample_mat[,!is.cell]);
	zeros <- my_background==0;
	# don't want zeros in the background.
	# set them to 1s
	norm_factor <- (sum(my_background)+sum(zeros))/sum(my_background)
	my_background <- my_background*norm_factor;
	my_background[zeros] <- 1;
	my_background <- my_background/sum(my_background);

	correct_cell <- function(c) {
		# Ignore zeros since we assume all genes expressed in the background to some
		# extent, thus the background distribution is far more distributed than the cell
		# due to sparsity. 
		# Assuming a cell is strictly made up of contamination + cell then this means 
		# the zeros cannot tell us anything about the contamination.
		# Expectation if all background
		#expected <- my_background;
		#expected[c==0] <- 0;
		#expected <- expected/sum(expected)*sum(c);

		# if not much background theen there should be many big difference
		# however there will also be differences due to chance, and the differences
		# due to chance are proportional to the total.

		# But we also expect those with many counts to have more non-background. 
		# Because of the scaling always expect about half to be -vs and half to be +vs
		# if same as background.
		#diff <- c-expected;

		# top most expressed genes in the background should be the guide to fitting the 
		# amount. 
		expected <- my_background;
		expected[c==0] <- 0;
		top_genes <- quantile(expected[expected>0], 0.80)

		top_bg_genes <- expected[expected >= top_genes]
		c_bg_genes <- c[names(c) %in% names(top_bg_genes)]
		bg_factor <-sum(c_bg_genes)/sum(top_bg_genes)
		expected <- my_background*bg_factor;

		diff <- c-expected;
		diff[diff < 0] <- 0
		return(diff)
	}
	corrected_mat <- apply(sample_mat, 2, correct_cell);
	return(corrected_mat);
}


saveRDS(my_background_correction(mydata, is.cell), paste(name, "emptyDrops_correctedtable.rds", sep="_"));
