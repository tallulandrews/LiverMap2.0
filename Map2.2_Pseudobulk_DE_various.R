args <- commandArgs(trailingOnly=TRUE)
# Seurat Object
# projname
# cluster_col_name (Use_clusters)

require(Seurat)
obj <- readRDS(args[1])
source("/cluster/home/tandrews/scripts/LiverMap2.0/Map2.2_add_metadata.R")
meta_data <- get_metadata(obj@meta.data$donor);

obj@meta.data <- cbind(obj@meta.data, meta_data);

if (length(args) >= 4) {
	obj@meta.data$Use_clusters <- obj@meta.data[,args[4]]
}
if (length(args) >= 5) {
	obj@meta.data$donor <- obj@meta.data[,args[5]]
}

exclude1 <- is.na(obj@meta.data$Use_clusters) | obj@meta.data$Use_clusters == "";
exclude2 <- is.na(obj@meta.data$donor) | obj@meta.data$donor == "";
obj <- obj[,!(exclude1 | exclude2)];

## Pseudobulk ##
source("~/scripts/LiverMap2.0/My_R_Scripts.R")

bulks <- get_pseudobulk(obj@assays$RNA@counts, 
			factor(obj@meta.data$Use_clusters), 
			factor(obj@meta.data$donor))
labs <- strsplit(colnames(bulks), "_")
c_labs <- sapply(labs, function(x){unlist(x[[1]])})
d_labs <- sapply(labs, function(x){unlist(x[[2]])})

# DE straight pseudobulks


## Rescale within each donor ##
c_freqs <- unlist(table(factor(obj@meta.data$Use_clusters)));
c_freqs <- c_freqs/sum(c_freqs);
norm_bulks <- bulks
for (d in unique(d_labs)) {
	dmat <- bulks[,d_labs==d]
	if(ncol(dmat) != length(c_freqs)) {
		print(paste("Error:", d, "only has", ncol(dmat), "cluster means when expect", length(c_freqs)));
	}
	norm <- t(log2(t(dmat)/colSums(dmat)*10000+1)) # clusters as columncs
	weighted <- t(t(norm)*c_freqs) 
	scale_factor <- rowMeans(weighted);
	mean_centered <- norm-scale_factor;
	norm_bulks[,d_labs==d] <- mean_centered;
}

# DE sample-specific mean_centered pseudobulks




# DE list:aaiii 
# Each Cluster : Cluster vs all others, M vs F and Elderly vs Young in each cluster - conditioning on each other. Ac


if(!identical(c_labs, d_labs)) {
require("edgeR")
edger_obj <- DGEList(bulks, samples=data.frame(cluster=c_labs, donor=d_labs), group=c_labs)
edger_obj <- calcNormFactors(edger_obj)
design <- model.matrix(~cluster+donor, data=edger_obj$samples)
design <- design[,colSums(design) > 0];

coef_names <- colnames(design);

edger_file <- paste(args[3], "_pseudobulk_edger.rds", sep="")
if (!file.exists(edger_file)) {
        edger_obj <- estimateDisp(edger_obj, design)
        fit <- glmQLFit(edger_obj, design)
        saveRDS(fit, edger_file);
} else {
        fit <- readRDS(edger_file)
}

contrast_vec <-rep(0, length(coef_names));
this_cluster <- paste("cluster", args[2], sep="");
contrast_vec[grepl("cluster", coef_names)] <- -1;
contrast_vec[1] <- -1
if (this_cluster %in% coef_names) {
        contrast_vec[coef_names==this_cluster]<- 1;
} else {
        contrast_vec[1] <- 1;
}

this_vs_all <- glmQLFTest(fit, contrast=contrast_vec)
out <- topTags(this_vs_all, nrow(obj));
all_DE <- list(seur_wilcox=out_seur_wilcox, seur_mast=out_seur_MAST, edger=out)
saveRDS(all_DE, paste(args[3], args[2], "_pseudobulk_DE.rds", sep="_"))
}

#######   END ########


### Diffcyt? ###
#require(diffcyt)

## make input:
#input_list <- list()
#experiment_info <- c()
#for (d in unique(obj@meta.data$donor)) {
#	for (c in unique(obj@meta.data$Use_clusters)) {
#		this <- obj@assays$RNA@counts[,obj@meta.data$donor ==d & obj@meta.data$Use_cluster ==c]
#		sample <- paste(d, c, sep="_")
#		if (sum(obj@meta.data$donor ==d & obj@meta.data$Use_cluster ==c) > 0) {
#			input_list[[sample]] <- as.matrix(t(this));
#			experiment_info <- rbind(experiment_info, c(sample, d, c))
#		}
#	}
#}
#colnames(experiment_info) <- c("sample_id", "patient_id", "cell_type")
#
##input data - list of matrices for each 'sample'
##experiment_info - dataframe of metadata for samples
#
##marker_info - dataframe of metadata for markers
#
##preprocess
#d_se <- prepareData(input_list, data.frame(experiment_info), marker_info=data.frame(gene=rownames(obj), marker_name=rownames(obj)))
##d_se <- transformData(d_se) <- don't do this b/c designed for CyTOF
##d_se <- generateClusters(d_se, seed_clustering = 123) <- FlowSOM clustering
## Calculate cluster cell counts
#d_counts <- calcCounts(d_se)
#
#d_se@elementMetadata$cluster_id <- matched_cell_cell_type
#
## Calculate cluster medians
#d_medians <- calcMedians(d_se)
#
## Create design matrix
## note: selecting columns containing group IDs and patient IDs (for an 
## unpaired dataset, only group IDs would be included)
#design <- createDesignMatrix(
#  data.frame(experiment_info), cols_design = c("cell_type", "patient_id")
#)
## Create contrast matrix
#coef_names <- colnames(design)
#contrast_vec <-rep(0, length(coef_names));
#this_cluster <- paste("cell_type", args[2], sep="");
#contrast_vec[grepl("cell_type", coef_names)] <- -1;
#contrast_vec[1] <- -1
#if (this_cluster %in% coef_names) {
#        contrast_vec[coef_names==this_cluster]<- 1;
#} else {
#        contrast_vec[1] <- 1;
#}
#contrast <- createContrast(contrast_vec)
#
#
## check
#nrow(contrast) == ncol(design)
#
## Test for differential abundance (DA) of clusters
#res_DA <- testDA_edgeR(d_counts, design, contrast)
#
## display table of results for top DA clusters
#topTable(res_DA, format_vals = TRUE)
#
#
