All_meta <- readRDS("/cluster/projects/macparland/TA/LiverMap2.0/RawData/Analysis/EmptyDrops/Merged_EmptyOnly_obj_Map2.2_ImportedClusters_ManualAnno_MetaDataOnly.rds")
subcluster <- Sys.glob("/cluster/projects/macparland/TA/LiverMap2.0/RawData/Analysis/EmptyDrops/Subcluster/ManualAnnotation/*fullmeta*.rds")

cell_ID_celltype <- All_meta[,c("cell_ID", "Coarse_clusters", "Coarse_Manual_Anno")]

cell_ID_celltype[,"Coarse_Manual_Anno"] <- as.character(cell_ID_celltype[,"Coarse_Manual_Anno"])

# Hepatocyte1 - using original coarse cluster
cell_ID_celltype[cell_ID_celltype[,2] == "0","Coarse_Manual_Anno"] <- "PortalHep"
cell_ID_celltype[cell_ID_celltype[,2] == "1","Coarse_Manual_Anno"] <- "PortalHep"
cell_ID_celltype[cell_ID_celltype[,2] == "4","Coarse_Manual_Anno"] <- "InterzonalHep"
cell_ID_celltype[cell_ID_celltype[,2] == "5","Coarse_Manual_Anno"] <- "CentralHep"
cell_ID_celltype[cell_ID_celltype[,2] == "9","Coarse_Manual_Anno"] <- "CentralHep"
cell_ID_celltype[cell_ID_celltype[,2] == "19","Coarse_Manual_Anno"] <- "ProlifInterzonalHep"

for( f in subcluster) {
	meta.data <- readRDS(f)
	tag <- strsplit(f, "_")[[1]][1]
	tag <- unlist(strsplit(tag, "/")); tag <- tag[length(tag)]
	subset <- cell_ID_celltype[,1] %in% meta.data$cell_ID
	meta.data <- meta.data[match(cell_ID_celltype[,1], meta.data$cell_ID),]
	cell_ID_celltype[subset,"Coarse_Manual_Anno"] <- paste(tag, meta.data[subset, "Subcluster_Manual"], sep="_")
}

cell_ID_celltype[,"Coarse_Manual_Anno"]<- gsub("\\+", "_", cell_ID_celltype[,"Coarse_Manual_Anno"])
	
# Set-up for CellPhonedb

full_counts <- readRDS("/cluster/projects/macparland/TA/LiverMap2.0/RawData/Analysis/EmptyDrops/Merged_EmptyOnly_obj_Map2.2_ImportedClusters_NormCounts.rds")

for (s in unique( All_meta$sample)) {
	cells <- All_meta$cell_ID[All_meta$sample == s]
	counts <- full_counts[,All_meta$cell_ID %in% cells]
	write.table(as.matrix(counts), paste(s, "counts.txt", sep="_"), col.names=T, row.names=T, sep="\t")
	meta <- cell_ID_celltype[cell_ID_celltype[,1] %in% cells, c("cell_ID", "Coarse_Manual_Anno")]
	colnames(meta) <-  c("Cell", "cell_type")
	write.table(meta, paste(s, "meta.txt", sep="_"), col.names=T, row.names=F, quote=F, sep="\t")
}

# Version 2 of Annotation

cell_ID_celltype[cell_ID_celltype[,2] == "19","Coarse_Manual_Anno"] <- "ProlifInterzonalHep"
for (s in unique( All_meta$sample)) {
	cells <- All_meta$cell_ID[All_meta$sample == s]
	meta <- cell_ID_celltype[cell_ID_celltype[,1] %in% cells, c("cell_ID", "Coarse_Manual_Anno")]
	colnames(meta) <-  c("Cell", "cell_type")
	write.table(meta, paste(s, "meta_general.txt", sep="_"), col.names=T, row.names=F, quote=F, sep="\t")
}

