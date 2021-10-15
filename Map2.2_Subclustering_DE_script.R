#require(edgeR)
source("~/scripts/LiverMap2.0/My_R_Scripts.R")
require(matrixStats)

files <- Sys.glob("*subcluster_pseudobulks.rds")
f_names <- unlist(lapply(strsplit(files, "_"), function(x){x[[1]]}));


for (i in 4:length(files)) {

bulks <- readRDS(files[i])
proj_name <- f_names[i]

pseudobulk_mat <- bulks$pseudobulk_mat
metadata <- bulks$metadata

metadata$cluster <- unlist(lapply( strsplit( colnames(pseudobulk_mat), "_"), function(x){unlist(x[1])}))

# Issue = Background / Batch effects!
# Solution = scale by each donor.
# Then use glm (not NB anymore).

# Scale & normalize pseudobulks!
pseudobulk_mat <- pseudobulk_mat[!grepl("^MT-", rownames(pseudobulk_mat)),]
pseudobulk_norm <- pseudobulk_mat #t(t(pseudobulk_mat)/colSums(pseudobulk_mat) * min(colSums(pseudobulk_mat)))
pseudobulk_scaled <- pseudobulk_norm

#avgs <- group_rowmeans(pseudobulk_norm, metadata$donor)
#vars <- group_rowvars(pseudobulk_norm, metadata$donor)
#vars <- sqrt(vars)
#vars[vars == 0] <- 0.0001
#for(d in levels(metadata$sample)) {
#	d_col <- metadata$sample == d
#	d_z_corrected <- (pseudobulk_norm[,d_col]-avgs[,d])/vars[,d]
#	pseudobulk_scaled[,d_col] <- d_z_corrected;
#}


out_list <- list()

for(type in unique(metadata$cluster)) {
	mat <- pseudobulk_mat[,metadata$cluster == type]
	non_mat <- pseudobulk_mat[,metadata$cluster != type]
	meta <- metadata[metadata$cluster == type,]
	non_meta <- metadata[metadata$cluster != type,]

	# Straight Means
	if (sum(meta$assay_type == "3pr") < 2) {
		this_3pr_mean <- rep(NA, nrow(mat));
		this_3pr_var <- rep(NA, nrow(mat));
	} else {
		this_3pr_mean <- rowMeans(mat[,meta$assay_type == "3pr"])
		this_3pr_var <- rowVars(mat[,meta$assay_type == "3pr"])/ sqrt(sum(meta$assay_type == "3pr"))
	}
	if (sum(meta$assay_type == "5pr") < 2) {
		this_5pr_mean <- rep(NA, nrow(mat));
		this_5pr_var <- rep(NA, nrow(mat));
	} else {
		this_5pr_mean <- rowMeans(mat[,meta$assay_type == "5pr"])
		this_5pr_var <- rowVars(mat[,meta$assay_type == "5pr"])/sqrt(sum(meta$assay_type == "5pr"))
	}
	if (sum(non_meta$assay_type == "3pr") < 2) {
		other_3pr_mean <- rep(NA, nrow(mat));
		other_3pr_var <- rep(NA, nrow(mat));
	} else {
		other_3pr_mean <- rowMeans(non_mat[,non_meta$assay_type == "3pr"])
		other_3pr_var <- rowVars(non_mat[,non_meta$assay_type == "3pr"])/ sqrt(sum(non_meta$assay_type == "3pr"))
	}
	if (sum(non_meta$assay_type == "5pr") < 2) {
		other_5pr_mean <- rep(NA, nrow(mat));
		other_5pr_var <- rep(NA, nrow(mat));
	} else {
		other_5pr_mean <- rowMeans(non_mat[,non_meta$assay_type == "5pr"])
		other_5pr_var <- rowVars(non_mat[,non_meta$assay_type == "5pr"])/sqrt(sum(non_meta$assay_type == "5pr"))
	}

	out <- cbind(this_3pr_mean, this_5pr_mean, other_3pr_mean, other_5pr_mean)
	score <- (out[,1]-out[,3]) + (out[,2]-out[,4])
	out <- cbind(out, score)
	colnames(out) <- c(paste(type, "3pr", sep="_"), paste(type, "5pr", sep="_"), "Other_3pr", "Other_5pr", "Sum_Diff_Both")
	out <- out[order( score ), ]
	out_list[[type]] <- out;	

	#assay <- metadata$assay_type
	#cluster <- metadata$cluster

	# gene-specific DE testing:
	#x <- mat["IL10RA",]
	#cluster_de <- function(x) {
	#	test = summary(glm(x~assay+cluster))
	#}
}		

saveRDS(out_list, paste(proj_name, "subcluster_pseudobulkDE_results.rds", sep="_"));
}
