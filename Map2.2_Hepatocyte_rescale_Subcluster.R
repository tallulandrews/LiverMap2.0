dir ="/cluster/projects/macparland/TA/LiverMap2.0/RawData/Analysis/EmptyDrops/Subcluster"
source("~/scripts/LiverMap2.0/My_R_Scripts.R")

obj <- readRDS("../Merged_EmptyOnly_obj_Map2.2_ImportedClusters.rds")

require(Seurat)

all_scaled <- c();
my_cell_ids <- c();
for(s in unique(obj@meta.data$sample)) {
	norm_mat <- obj@assays$RNA@data[,obj@meta.data$sample == s]
	type <- obj@meta.data$Core_clusters[obj@meta.data$sample == s]
	pseudobulks <- group_rowmeans(norm_mat, type)
	weighted_mean <- rowMeans(pseudobulks[,1:10])
	stddevs <- apply(norm_mat, 1, sd)

	scaled <- (norm_mat - weighted_mean)/stddevs
	my_cell_ids <- c(my_cell_ids, obj@meta.data$cell_ID[obj@meta.data$sample == s])
	all_scaled <- cbind(all_scaled, as.matrix(scaled));
}

identical(my_cell_ids, obj@meta.data$cell_ID)

obj@assays$RNA@scale.data <- all_scaled

saveRDS(obj, "Merged_EmptyOnly_obj_Map2.2_ImportedClusters_RescaledPseudobulk.rds")

obj_hepatocyte <- obj[,obj@meta.data$Coarse_clusters %in% c(0,1,4,5,9,19)]
obj_hepatocyte2 <- obj[,obj@meta.data$Coarse_clusters %in% c(2,17,18)]
rm(obj)

# copied from: Map2.2_IndiScaled_Subcluster.R
subset = obj_hepatocyte2

gene_filter <- group_rowmeans(subset@assays$RNA@counts, subset@meta.data$donor);
gene_filter <- apply(gene_filter, 1, median)
gene_filter <- gene_filter > 0;

# Subcluster
subset <- FindVariableFeatures(subset[gene_filter,], selection.method="vst", nfeatures=2500)
VariableFeatures(subset) <- unique(c(key_marker_sets[[set_i]], VariableFeatures(subset)));
subset <- RunPCA(subset, features=VariableFeatures(object=subset))


set.seed(2009)
require(harmony)
subset <- RunHarmony(subset, "sample", plot_convergence=FALSE)

npcs <- 10

res <- seq(from=0.3, to=1.5, by=0.2)
nkNN <- seq(from=30, to=60, by=10)

for(res_param in res) {
for(nkNN_param in nkNN){
        subset <- FindNeighbors(subset, reduction="harmony", dims = 1:npcs, k.param=nkNN_param)
        subset <- FindClusters(subset, reduction="harmony", resolution = res_param, k.param=nkNN_param)
        name <- paste("knn_",nkNN_param,"_res_", res_param, sep="");

        subset@meta.data[[name]] <- subset@meta.data$seurat_clusters;

}}

saveRDS(subset, "Hepatocyte2_Rescaled_allgenes_Subcluster.rds")
