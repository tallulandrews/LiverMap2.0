require("Seurat")
source("/cluster/home/tandrews/scripts/LiverMap2.0/My_R_Scripts.R")
require("DropletUtils")

loc_file <- read.delim("/cluster/home/tandrews/scripts/LiverMap2.0/TCR_BCR_clonetype_files.csv", sep=",", header=T, stringsAsFactors=FALSE)

MAX_CELLS=25000
MIN_CELLS=100

i = 1;

name=loc_file[i,1]
# Read and TCRs & BCRs
tcrs <- read.delim(loc_file[i,"TCR_file"], sep=",", header=T)
bcrs <- read.delim(loc_file[i,"BCR_file"], sep=",", header=T)

# Read raw mat and tidy
filt <- Read10X(loc_file[i,"X5pr_filtdir"])
raw_mat <- Read10X(loc_file[i,"X5pr_rawdir"])
raw_mat <- raw_mat[,Matrix::colSums(raw_mat) > 0]
raw_mat <- raw_mat[Matrix::rowSums(raw_mat) > 0,]
print(paste("Stats out:", dim(raw_mat), "input", c("genes","cells"),"of", name))
colnames(raw_mat) <- paste(colnames(raw_mat), "1", sep="-")

# Get total UMIs that correspond to the min & max thresholds.
br.out <- barcodeRanks(raw_mat)
umi_per_barcode <- Matrix::colSums(raw_mat)



# get lower threshold for cells
plausible_cell_threshold <- max(br.out$total[br.out$rank > MAX_CELLS]);
t_cell_threshold <- quantile(umi_per_barcode[colnames(raw_mat) %in% tcrs[,1]], 0.25)
b_cell_threshold <- quantile(umi_per_barcode[colnames(raw_mat) %in% bcrs[,1]], 0.25)
lower_threshold = min(plausible_cell_threshold, t_cell_threshold, b_cell_threshold)
# get higher threshold for cells
mandatory_cell_threshold <- max(br.out$total[br.out$rank < MIN_CELLS]);
t_cell_threshold <- min(umi_per_barcode[colnames(raw_mat) %in% tcrs[tcrs$high_confidence=="True" & tcrs$is_cell == "True",1]])
b_cell_threshold <- min(umi_per_barcode[colnames(raw_mat) %in% bcrs[bcrs$high_confidence=="True" & bcrs$is_cell=="True",1]])
higher_threshold = min(mandatory_cell_threshold)

# Run EmptyDrops to try to clean up while allowing more permissive thershold to keep tcrs and bcrs

set.seed(100)
# Run EmptyDrops
e.out <- emptyDrops(raw_mat, lower=lower_threshold, niters=100000, ignore=3, retain=higher_threshold)

# Clean up results
e.out <- e.out[!is.na(e.out$PValue),] # Remove NAs (too few UMIs to estimate PValue)
e.out$q.value <- e.out$PValue;
e.out$q.value[e.out$Limited] <- 0; # Set those with p < 1/n_simulations to p = 0 to ensure they are kept after multiple testing correction
e.out$q.value <- p.adjust(e.out$q.value, method="fdr"); # apply FDR multiple testing correction.

keep_cells <- rownames(e.out)[e.out$q.value < 0.01]
filt_mat <- raw_mat[,match(keep_cells, colnames(raw_mat))]
tcr_meta <- tcrs[tcrs$barcode %in% keep_cells,]
bcr_meta <- bcrs[bcrs$barcode %in% keep_cells,]

# Can't turn tcrs into a metadata table b/c more than one TCR/cell...

saveRDS(list(mat=filt_mat, tcr=tcr_meta, bcr=bcr_meta), file=paste(name,"TCR_BCR_EmptyDrops.rds", sep="_"))


