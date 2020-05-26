require("Seurat")
obj <- readRDS("Merged_all_genes_plus_metadata.rds")
obj$cell_id <- paste(obj$orig.ident, obj$cell_barcode, sep="_")
obj_analyzed <- readRDS("../integration_harmony_plus_analysis.rds")
obj_analyzed$cell_id <- paste(obj_analyzed$orig.ident, obj_analyzed$cell_barcode, sep="_")
reorder <- match(obj$cell_id, obj_analyzed$cell_id)



obj@reductions <- obj_analyzed@reductions

for( reduc in names(obj@reductions)) {
	obj@reductions[[reduc]]@cell.embeddings <- obj@reductions[[reduc]]@cell.embeddings[reorder,]
}

obj@meta.data <- cbind(obj@meta.data, obj_analyzed@meta.data[reorder,])

VariableFeatures(obj) <- VariableFeatures(obj_analyzed);

saveRDS(obj, file="Merged_all_genes_plus_analysis.rds")

png("160k_cells_umap_w_cellcycle.png", width=6, height=6, units="in", res=300)
DimPlot(obj, reduction="umap", group.by="CCPhase")
dev.off()
