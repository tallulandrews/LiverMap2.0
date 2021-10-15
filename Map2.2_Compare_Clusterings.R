require("Seurat")

proj_name="Map2.2_merged_integrated"

liver.integrated <- readRDS(paste(proj_name, "harmony_plus_clusters.rds", sep="_"))

# Compare all these clusterings
require(igraph)
require(gplots)
clust_table <- liver.integrated@meta.data[, grepl("^knn_", colnames(liver.integrated@meta.data))]

out <- scClustViz::CalcAllSCV(liver.integrated, clust_table, assayType="RNA", assaySlot="scale.data", 
