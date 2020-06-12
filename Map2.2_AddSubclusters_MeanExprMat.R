
require(Seurat)

all <- readRDS("All_genes_Merged_obj_v2_with_Analysis.rds")
all@meta.data$Use_clusters <- all@meta.data$Fine_clusters
all@meta.data$Sub_clusters <- paste("All", all@meta.data$Fine_clusters,sep="_")
sub <- Sys.glob("Subcluster/*analysis.rds")
all$cell_ID <- paste(all$sample, all$cell_barcode, sep="_")

for (i in 1:length(sub)) {
	sub_obj <- readRDS(sub[i])
	sub_obj$cell_ID <- paste(sub_obj$sample, sub_obj$cell_barcode, sep="_")
	sub_name <- strsplit(gsub("Subcluster/","", sub[i],), "_")[[1]][[1]]
	sub_obj@meta.data$Unique_cluster_lab <- paste(sub_name, sub_obj@meta.data$Use_clusters, sep="_")
	all@meta.data[match(sub_obj@meta.data$cell_ID, all@meta.data$cell_ID),"Sub_clusters"] <- sub_obj@meta.data$Unique_cluster_lab
}

saveRDS(all, "All_genes_Merged_obj_v2_with_Analysis2.rds")

source("~/scripts/LiverMap2.0/My_R_Scripts.R")

subcluster_mat <- get_rel_expression(all@assays$RNA@data, all@meta.data$Sub_clusters, all@meta.data$sample);
generalcluster_mat <- get_rel_expression(all@assays$RNA@data, all@meta.data$Use_clusters, all@meta.data$sample);

huge_mat <- cbind(generalcluster_mat, subcluster_mat);
write.table(huge_mat, "All_cluster_subcluser_mean_expr_mat.csv", sep=",", col.names=TRUE, row.names=TRUE)
