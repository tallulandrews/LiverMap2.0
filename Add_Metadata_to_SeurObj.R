# Auto-annotation
liver.integrated <- readRDS("all_gene_merged_obj.rds");

all_CC <- readRDS("All20_cellcycle.rds")
#all_doublets <- readRDS("All20_doubletDetection.rds")

liver.integrated@meta.data$CCPhase <- rep("unknown", ncol(liver.integrated));
liver.integrated@meta.data$DoubletScore <- rep("unknown", ncol(liver.integrated));
for (donor in unique(liver.integrated@meta.data$orig.ident)) {
        cell_ids <- liver.integrated@meta.data$cell_barcode[liver.integrated@meta.data$donor == donor]

        CC <- all_CC[[donor]];
        CC <- CC[CC$cell_barcode %in% cell_ids,]
        CC <- CC[match(cell_ids, CC$cell_barcode),]

#        doublets <- all_doublets[[donor]];
#        doublets <- doublets[doublets$cell_barcode %in% cell_ids,]
#        doublets <- doublets[match(cell_ids, doublets$cell_barcode),]

        liver.integrated@meta.data$CCPhase[liver.integrated@meta.data$orig.ident == donor] <- as.character(CC$Phase)

#        liver.integrated@meta.data$DoubletScore[liver.integrated@meta.data$orig.ident == donor] <- as.character(doublets$DoubletScore)
}

saveRDS(liver.integrated, "Merged_all_genes_plus_metadata.rds")

