require("Seurat")
source("~/scripts/LiverMap2.0/My_R_Scripts.R")

args <- commandArgs(trailingOnly=TRUE);

liver.integrated <- readRDS(args[1])
proj_name <- args[2]

# require all genes to be detected in at least one cell in a majority of donors
gene_filter <- group_rowmeans(liver.integrated@assays$RNA@counts, liver.integrated@meta.data$donor);
gene_filter <- apply(gene_filter, 1, median)
gene_filter <- gene_filter > 0;

liver.integrated <- liver.integrated[gene_filter,]

set.seed(1039)

# Dimensionality Reduction
### Probably need to re-run harmony on just these cell for best results....####
liver.integrated <- NormalizeData(liver.integrated)
liver.integrated <- ScaleData(liver.integrated);
liver.integrated <- FindVariableFeatures(liver.integrated, selection.method = "vst", nfeatures = 2500)
liver.integrated <- RunPCA(liver.integrated, features = VariableFeatures(object = liver.integrated))
ElbowPlot(liver.integrated)



require("harmony")
set.seed(10131)
liver.integrated <- RunHarmony(liver.integrated, "donor", plot_convergence = TRUE)


npcs <- as.numeric(15)

# Cluster with many different parameters
res <- seq(from=0.3, to=1.5, by=0.2)
nkNN <- seq(from=30, to=60, by=10)

for(res_param in res) {
for(nkNN_param in nkNN){
liver.integrated <- FindNeighbors(liver.integrated, reduction="harmony", dims = 1:npcs, k.param=nkNN_param)
liver.integrated <- FindClusters(liver.integrated, reduction="harmony", resolution = res_param, k.param=nkNN_param)
name <- paste("knn_",nkNN_param,"_res_", res_param, sep="");

liver.integrated@meta.data[[name]] <- liver.integrated@meta.data$seurat_clusters;

}}

# Compare all these clusterings
require(igraph)
require(gplots)
clust_table <- liver.integrated@meta.data[, grepl("^knn_", names(liver.integrated@meta.data))]
clust_table <- as.matrix(apply(clust_table,2,as.numeric))
require("proxy")
clust_dists <- proxy::dist(clust_table, method=function(x,y){igraph::compare(x,y,method="vi")}, by_rows=FALSE)


# Find robust exemplar clustering(s)
require("apcluster")
require("gplots")
set.seed(18371)

res1 <- apcluster(-1*as.matrix(clust_dists), p=-2.5)

ex_clust_sz <- unlist(lapply(res1@clusters, length))

top <-which( ex_clust_sz == max(ex_clust_sz));
ex_clust_sz[top] <- 0
top2 <-which( ex_clust_sz == max(ex_clust_sz));
top <- names(res1@exemplars)[top]
top2 <- names(res1@exemplars)[top2]

use_lvl = top;


#manually select which exemplar to use
liver.integrated@meta.data$Use_clusters <- liver.integrated@meta.data[[use_lvl]]

#apcluster::heatmap(res1, -1*as.matrix(clust_dists))

png(paste(proj_name, "compare_clusterings_heatmap.png", sep="_"), width=6, height=6, units="in", res=300)
lab <- matrix("", ncol=ncol(clust_table), nrow=ncol(clust_table))
lab[colnames(clust_table)==use_lvl, colnames(clust_table)==use_lvl] <- "1"
lab[colnames(clust_table)==top2, colnames(clust_table)==top2] <- "2"
heatmap.2(as.matrix(clust_dists), trace="none", distfun=function(x){return(as.dist(clust_dists))}, cellnote=lab)
dev.off()


# Visualize the Chosen clusterings
nkNN <- as.numeric(40)

set.seed(1092)

liver.integrated <- RunTSNE(liver.integrated, reduction="harmony",  dims = 1:npcs)
liver.integrated <- RunUMAP(liver.integrated, reduction="harmony",  dims = 1:npcs, parallel=FALSE, n.neighbour=nkNN)
png(paste(proj_name, "_umap.png", sep="_"), width=6, height=6, units="in", res=50)
DimPlot(liver.integrated, reduction = "umap", group.by="Use_clusters")
dev.off()
png(paste(proj_name, "_tsne.png", sep="_"), width=6, height=6, units="in", res=50)
DimPlot(liver.integrated, reduction = "tsne", group.by="Use_clusters")
dev.off()
png(paste(proj_name, "_umap_donor.png", sep="_"), width=6, height=6, units="in", res=50)
DimPlot(liver.integrated, reduction = "umap", group.by="orig.ident")
dev.off()
DimPlot(liver.integrated, reduction = "tsne", group.by="orig.ident")

barplot(table(liver.integrated@meta.data$orig.ident, liver.integrated@meta.data$Use_clusters), col=rainbow(20))

tab <- table(liver.integrated@meta.data$orig.ident, liver.integrated@meta.data$Use_clusters)

png(paste(proj_name, "_umap_autoanno.png", sep="_"), width=7.5, height=6, units="in", res=50)
DimPlot(liver.integrated, reduction = "umap", group.by="consistent_labs")
dev.off()

saveRDS(liver.integrated, paste(proj_name, "_analysis.rds", sep=""))

## Pseudobulk
source("~/scripts/LiverMap2.0/My_R_Scripts.R")
pseudobulks <- get_pseudobulk(liver.integrated@assays$RNA@data, liver.integrated@meta.data$Use_clusters, liver.integrated@meta.data$sample);
saveRDS(pseudobulks, paste(proj_name, "_pseudobulks.rds", sep=""))

