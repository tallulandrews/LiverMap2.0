# Auto-annotation
liver.integrated <- readRDS("All_merged_universal_genes_harmony_integrated_v2.rds");
prefix <- "All_harmony_integrated";

all_doublets <- readRDS("All20_doubletDetection2.rds")

liver.integrated@meta.data$DoubletScore <- rep(-1, ncol(liver.integrated));
liver.integrated@meta.data$is.Doublet <- rep("unknown", ncol(liver.integrated));

require(stringr)
for (sample in unique(liver.integrated@meta.data$orig.ident)) {
        cell_ids <- liver.integrated@meta.data$cell_barcode[liver.integrated@meta.data$sample == sample]

	sample2 <- gsub("_NPC", "-NPC", sample)
        doublets <- all_doublets[[sample2]];
        doublets <- doublets[doublets$cell_barcode %in% cell_ids,]
        doublets <- doublets[match(cell_ids, doublets$cell_barcode),]

	doublets$DoubletScore <- as.numeric(doublets[,ncol(doublets)-1])
	doublets$is.Doublet <- doublets[,ncol(doublets)-1]

        liver.integrated@meta.data$DoubletScore[liver.integrated@meta.data$orig.ident == sample] <- as.numeric(doublets$DoubletScore)
        liver.integrated@meta.data$is.Doublet[liver.integrated@meta.data$orig.ident == sample] <- as.character(doublets$is.Doublet)
}

require("Seurat")
png(paste(prefix, "DoubletScore.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::FeaturePlot(liver.integrated, "DoubletScore")
dev.off()

png(paste(prefix, "isDoublet.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::DimPlot(liver.integrated, group.by="is.Doublet")
dev.off()
