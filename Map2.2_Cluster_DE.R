require("Seurat")
source("~/scripts/LiverMap2.0/My_R_Scripts.R")

# Read in Data
obj <- readRDS("Merged_EmptyOnly_obj_Map2.2_ImportedClusters_ManualAnno_dimreduce.rds")

data_3pr <- obj[,obj@meta.data$assay_type == "3pr"]
#data_5pr <- obj[,obj@meta.data$assay_type == "5pr"]


# 20 clusters
# Manual Anno cell types

# Seurat Wilcox.test
for (i in unique(obj@meta.data$Coarse_clusters)) {
	res <- FindMarkers(data_3pr, group.by="Coarse_clusters", ident.1=i)
	write.table(res, file=paste("DE", i, "Coarse_clusters", "wilcox.csv", sep="_"), sep=",", row.names=T, col.names=T)
}

for (i in unique(obj@meta.data$Coarse_Manual_Anno)) {
	res <- FindMarkers(data_3pr, group.by="Coarse_Manual_Anno", ident.1=i)
	write.table(res, file=paste("DE", i, "Coarse_Manual", "wilcox.csv", sep="_"), sep=",", row.names=T, col.names=T)
}


# Pseudobulk edgeR test
#pseudo_3pr <- get_pseudobulk(data_3pr@assays$RNA@counts, data_3pr@meta.data$Coarse_clusters, data_3pr@meta.data$sample)
#pseudo_5pr <- get_pseudobulk(data_5pr, factor(data_5pr@meta.data$Coarse_clusters), factor(data_5pr@meta.data$sample))


