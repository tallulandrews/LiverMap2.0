require("Seurat")

liver.integrated <- readRDS("basic_integration_analysis.rds")

set.seed(7742)
#liver.integrated <- RunPCA(liver.integrated, features = VariableFeatures(object = liver.integrated))
#ElbowPlot(liver.integrated)
npcs <- 10

res <- seq(from=0.3, to=2, by=0.2)
nkNN <- seq(from=30, to=90, by=10)

for(res_param in res) {
for(nkNN_param in nkNN){
liver.integrated <- FindNeighbors(liver.integrated, dims = 1:npcs, k.param=nkNN_param)
liver.integrated <- FindClusters(liver.integrated, resolution = res_param, k.param=nkNN_param)
name <- paste("knn_",nkNN_param,"_res_", res_param, sep="");

liver.integrated@meta.data[[name]] <- liver.integrated@meta.data$seurat_clusters;

}}

require(igraph)
require(gplots)
clust_table <- liver.integrated@meta.data[, grepl("^knn_", names(liver.integrated@meta.data))]
clust_table <- as.matrix(apply(clust_table,2,as.numeric))
require("proxy")
clust_dists <- proxy::dist(clust_table, method=function(x,y){igraph::compare(x,y,method="vi")}, by_rows=FALSE)
clust_similr1 <- proxy::simil(clust_table, method=function(x,y){igraph::compare(x,y,method="nmi")}, by_rows=FALSE)
clust_similr2 <- proxy::simil(clust_table, method=function(x,y){igraph::compare(x,y,method="adjusted.rand")}, by_rows=FALSE)


require("apcluster")
set.seed(18371)

res1 <- apcluster(-1*as.matrix(clust_dists))
res2 <- apcluster(as.matrix(clust_similr1))
res3 <- apcluster(as.matrix(clust_similr1))

valid_clusterings <- res1@exemplars[which(res1@exemplars %in% res2@exemplars & res1@exemplars %in% res3@exemplars)]

#manually select which exemplar to use
liver.integrated@meta.data$Coarse_clusters <- liver.integrated@meta.data$knn_70_res_0.3
liver.integrated@meta.data$Fine_clusters <- liver.integrated@meta.data$knn_90_res_1.5

png("compare_clusterings_heatmap.png", width=6, height=6, units="in", res=300)
apcluster::heatmap(res1, -1*as.matrix(clust_dists))
dev.off()

saveRDS(liver.integrated, "Integrate_with_robust_clusters.rds")

#create a heatmap where: cell = average similarity of this clustering to all other clusterings
# distance between clusterings is measured using igraph::compare(c1, c2, method="vi")

require(clustree)
clustree(liver.integrated@meta.data, prefix="knn_70_res_")
clustree(liver.integrated@meta.data, prefix="knn_60_res_")

saveRDS(liver.integrated, file="alt_clusterings.rds")

