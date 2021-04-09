# Combine all mean files into one.

dir <- "/cluster/projects/macparland/TA/LiverMap2.0/RawData/Analysis/EmptyDrops/Subcluster/ManualAnnotation"


files <- Sys.glob("*_cellphonedb_means.txt")



all_mats <- list()
all_comparisons <- c()
all_interactors <- c()
for (f in files) {
	mat <- read.table(f, header=T, sep="\t")
	ids <- paste(mat[,1], mat[,2], sep="_")
	mat[,1] <- ids;
#	print(head(mat[,1:5])) # Debug
	tag <- sub("_cellphonedb_means.txt", "", f)
	all_mats[[tag]] <- mat
	
	all_interactors <- unique(rbind(all_interactors, mat[,1:11]))
	all_comparisons <- unique(c(all_comparisons, colnames(mat)[12:ncol(mat)]))
#	print(dim(mat)) # Debug
}


full_matrix <- c()
n_mat <- c();
sum_matrix <- c();
full_matrix_cell_types <- c()
full_matrix_sample <- c();

all_comparisons <- sort(all_comparisons)

# Count number ofsample in which each pair of cell-types are "interacting" 
# threshold of expression for a pair of proteins to be interacting
# threshold of number of pairs for a cell-type to be interating
expr_threshold = 0.5;
n_interaction_threshold = 5;
n_present <- rep(0, length(all_comparisons))
n_sample_inter <- rep(0, length(all_comparisons))

for (i in 1:length(all_mats)) {
	this_mat <- all_mats[[i]]
	this_mat <- this_mat[match( all_interactors[,1], this_mat[,1]),]
	rownames(this_mat) <- all_interactors[,1]
	tmp <- as.matrix(this_mat[,colnames(this_mat) %in% all_comparisons])
	tmp <- tmp[,match(all_comparisons, colnames(tmp))]
	colnames(tmp) <- paste(names(all_mats)[i], all_comparisons, sep="_")
	this_mat <- tmp;
	full_matrix_cell_types <- c(full_matrix_cell_types, all_comparisons)
	full_matrix_sample <- c(full_matrix_sample, rep(names(all_mats)[i], length(all_comparisons)))

	if (is.null(nrow(full_matrix))) {
		full_matrix <- this_mat
		this_mat[is.na(this_mat)] <- 0;
		sum_matrix <- this_mat;
		n_mat <- !is.na(this_mat) + 0;
	} else {
		full_matrix <- cbind(full_matrix, this_mat)
		n_mat <- n_mat + !is.na(this_mat);
		this_mat[is.na(this_mat)] <- 0;
		sum_matrix <- sum_matrix + this_mat;
	}

	# Count T/F interactions
	present <- !( colSums(is.na(tmp)) == nrow(tmp))
	n_present <- n_present + present
	n_interpairs <- colSums(tmp > expr_threshold, na.rm=T)
	n_sample_inter <- n_sample_inter + (n_interpairs > n_interaction_threshold)
}

mean_across_samples_mat <- sum_matrix/n_mat

overall_score <- colSums(mean_across_samples_mat)
cell_types <- sub("C37_", "", names(overall_score))
cell_types <- gsub("\\.qHSC", "_qHSC", cell_types)
cell_types <- gsub("\\.aHSC", "_aHSC", cell_types)

cell_types <- strsplit(cell_types, "\\.")

cell_type_2_id <- factor(unique(unlist(cell_types)))

heatmap_dat <- matrix(0, nrow=length(cell_type_2_id), ncol=length(cell_type_2_id))
heatmap_nsamp_dat <- matrix(0, nrow=length(cell_type_2_id), ncol=length(cell_type_2_id))
for(r in 1:length(cell_types)) {
	i <- which(levels(cell_type_2_id) == cell_types[[r]][1])
	j <- which(levels(cell_type_2_id) == cell_types[[r]][2])
	heatmap_dat[i,j] <- overall_score[r]
	heatmap_nsamp_dat[i,j] <- n_sample_inter[r]/n_present[r]
}

colnames(heatmap_dat) <- rownames(heatmap_dat) <- levels(cell_type_2_id)
colnames(heatmap_nsamp_dat) <- rownames(heatmap_nsamp_dat) <- levels(cell_type_2_id)

require(gplots)

png("Mean_pair_expr_prelim.png", width=10, height=10, units="in", res=300)
heatmap.2(heatmap_dat, margins=c(8, 8), trace="none")
dev.off()
png("Mean_pair_replicated_prelim.png", width=10, height=10, units="in", res=300)
heatmap.2(heatmap_nsamp_dat, margins=c(8, 8), trace="none")
dev.off()
