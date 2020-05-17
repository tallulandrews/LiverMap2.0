#args <- commandArgs(trailingOnly=TRUE)
# seurat object RDS file.
# plot file name prefix

require("Seurat")
#seur_obj <- readRDS(as.character(args[1]));
#prefix <- as.character(args[2]);

png(paste(prefix, "perMT.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::FeaturePlot(seur_obj, "percent.mt")
dev.off()
png(paste(prefix, "nFeature.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::FeaturePlot(seur_obj, "nFeature_RNA")
dev.off()
png(paste(prefix, "CCphase.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::DimPlot(seur_obj, group.by="Phase")
dev.off()
png(paste(prefix, "sample.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::DimPlot(seur_obj, group.by="sample")
dev.off()
