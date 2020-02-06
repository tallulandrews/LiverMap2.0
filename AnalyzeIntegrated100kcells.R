require("Seurat")

liver.integrated <- readRDS("basicintegration.rds")

set.seed(7742)
require("sctransform")
liver.integrated <- ScaleData(liver.integrated);
liver.integrated <- RunPCA(liver.integrated, features = VariableFeatures(object = liver.integrated))
ElbowPlot(liver.integrated)

set.seed(8291)
npcs <- 10
res <- 1.2
nkNN <- 20

liver.integrated <- FindNeighbors(liver.integrated, dims = 1:npcs, k.param=nkNN)
liver.integrated <- FindClusters(liver.integrated, resolution = res, k.param=nkNN)
liver.integrated <- RunTSNE(liver.integrated, dims = 1:npcs)
liver.integrated <- RunUMAP(liver.integrated, dims = 1:npcs, parallel=FALSE, n.neighbour=nkNN)
png("Highres_integrated_umap.png", width=6, height=6, units="in", res=50)
DimPlot(liver.integrated, reduction = "umap")
dev.off()
png("Highres_integrated_tsne.png", width=6, height=6, units="in", res=50)
DimPlot(liver.integrated, reduction = "tsne")
dev.off()
png("Donor_integrated_umap.png", width=6, height=6, units="in", res=50)
DimPlot(liver.integrated, reduction = "umap", group.by="orig.ident")
dev.off()
png("Donor_integrated_tsne.png", width=6, height=6, units="in", res=50)
DimPlot(liver.integrated, reduction = "tsne", group.by="orig.ident")
dev.off()

png("Integrated_clusters_by_donor.png", width=6, height=6, units="in", res=50)
barplot(table(liver.integrated@meta.data$orig.ident, liver.integrated@meta.data$seurat_clusters), col=rainbow(20))
dev.off()

#saveRDS(liver.integrated, "basic_integration_analysis.rds");

# Automated annotation? - doesn't work as too few genes left, but does work on original datasets. 
# automated annotation at individual dataset level then add it here.

source("Setup_autoannotation.R")
#liver.integrated <- readRDS("basic_integration_analysis.rds")

#test <- Use_markers_for_anno(liver.integrated@assays$RNA$data, liver.integrated@meta.data$seurat_clusters)

all_anno <- readRDS("All20_automatedannotation.rds");

#cell_barcodes <- unlist(lapply(strsplit(rownames(liver.integrated@meta.data), "_"), function(x){return(x[[1]])}))
#liver.integrated@meta.data$cell_barcode <- cell_barcodes

liver.integrated@meta.data$scmap_anno <- rep("unknown", ncol(liver.integrated));
liver.integrated@meta.data$scmap_anno2 <- rep("unknown", ncol(liver.integrated));
for (donor in unique(liver.integrated@meta.data$orig.ident)) {
	cell_ids <- liver.integrated@meta.data$cell_barcode[liver.integrated@meta.data$donor == donor]
	anno <- all_anno[[donor]];
	anno <- anno[anno$cell_barcode %in% cell_ids,]
	anno <- anno[match(cell_ids, anno$cell_barcode),]
	liver.integrated@meta.data$scmap_anno[liver.integrated@meta.data$orig.ident == donor] <- as.character(anno$scmap_cluster_anno$lm1)
	liver.integrated@meta.data$scmap_anno2[liver.integrated@meta.data$orig.ident == donor] <- as.character(anno$scmap_cell_anno)
}


png("AutoAnno_integrated_umap.png", width=8, height=6, units="in", res=50)
DimPlot(liver.integrated, reduction = "umap", group.by="scmap_anno")
dev.off()
png("AutoAnno_integrated_tsne.png", width=8, height=6, units="in", res=50)
DimPlot(liver.integrated, reduction = "tsne", group.by="scmap_anno")
dev.off()
png("AutoAnno2_integrated_umap.png", width=8, height=6, units="in", res=50)
DimPlot(liver.integrated, reduction = "umap", group.by="scmap_anno2")
dev.off()
png("AutoAnno2_integrated_tsne.png", width=8, height=6, units="in", res=50)
DimPlot(liver.integrated, reduction = "tsne", group.by="scmap_anno2")
dev.off()

saveRDS(liver.integrated, "basic_integration_analysis.rds");
