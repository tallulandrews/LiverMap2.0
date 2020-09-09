#args <- commandArgs(trailingOnly=TRUE)
# seurat object RDS file.
# plot file name prefix

require("Seurat")
#seur_obj <- readRDS(as.character(args[1]));
#prefix <- as.character(args[2]);

png(paste(prefix, "perMT.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::FeaturePlot(seur_obj, "percent.mt", reduction="umap")
dev.off()
png(paste(prefix, "nFeature.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::FeaturePlot(seur_obj, "nFeature_RNA", reduction="umap")
dev.off()
png(paste(prefix, "CCphase.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::DimPlot(seur_obj, group.by="Phase", reduction="umap")
dev.off()
png(paste(prefix, "sample.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::DimPlot(seur_obj, group.by="sample", reduction="umap")
dev.off()

png(paste(prefix, "assay.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::DimPlot(seur_obj, group.by="assay_type", reduction="umap")
dev.off()

# More Complex Plots

source("/cluster/home/tandrews/scripts/LiverMap2.0/Colour_Scheme.R")
source("/cluster/home/tandrews/scripts/LiverMap2.0/My_R_Scripts.R")

png(paste(prefix, "marker_anno.png", sep="_"), width=6, height=6, units="in", res=150)
Type_DimPlot(seur_obj, type_col="marker_labs", reduction="umap", cluster_col="seurat_clusters")
dev.off()
png(paste(prefix, "scmap_anno.png", sep="_"), width=6, height=6, units="in", res=150)
Type_DimPlot(seur_obj, type_col="consistent_labs", reduction="umap", cluster_col="seurat_clusters")
dev.off()

png(paste(prefix, "Seurclusters.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::DimPlot(seur_obj, group.by="seurat_clusters")
dev.off()


png(paste(prefix, "Sample_by_cluster.png", sep="_"), width=6, height=6, units="in", res=150)
perc <- table(seur_obj@meta.data$sample, seur_obj@meta.data$seurat_clusters)
perc <- t(t(perc)/colSums(perc))
barplot(perc, col=get_seurat_colours(seur_obj, "sample"), main="Sample by Cluster")
dev.off()

png(paste(prefix, "CCphase_by_cluster.png", sep="_"), width=6, height=6, units="in", res=150)
perc <- table(seur_obj@meta.data$Phase, seur_obj@meta.data$seurat_clusters)
perc <- t(t(perc)/colSums(perc))
barplot(perc, col=get_seurat_colours(seur_obj, "Phase"), main="CC by Cluster")
dev.off()

png(paste(prefix, "scmap_by_cluster.png", sep="_"), width=6, height=6, units="in", res=150)
perc <- table(map_cell_types(seur_obj@meta.data$consistent_labs), seur_obj@meta.data$seurat_clusters)
perc <- t(t(perc)/colSums(perc))
type_cols <- Cell_type_colours[match(rownames(perc),Cell_type_colours[,1]),2]
barplot(perc, col=type_cols, main="Scmap AutoAnno by Cluster")
dev.off()

png(paste(prefix, "marker_by_cluster.png", sep="_"), width=6, height=6, units="in", res=150)
perc <- table(map_cell_types(seur_obj@meta.data$marker_labs), seur_obj@meta.data$seurat_clusters)
perc <- t(t(perc)/colSums(perc))
type_cols <- Cell_type_colours[match(rownames(perc),Cell_type_colours[,1]),2]
barplot(perc, col=type_cols, main="Marker AutoAnno by Cluster")
dev.off()
