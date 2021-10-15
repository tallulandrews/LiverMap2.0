# Criteria :
#	- 3' RNAseq
#	- Not Flush
#	- Not biased towards specific cell-types / clusters: 
#	- Downsample keeping some of all cell-types

set.seed(101)

dir="/cluster/projects/macparland/TA/LiverMap2.0/RawData/Analysis/EmptyDrops/"

all_obj <- readRDS("Merged_EmptyOnly_obj_Map2.2_ImportedClusters_ManualAnno.rds")

all_obj <- all_obj[,all_obj@meta.data$assay_type == "3pr"]
all_obj <- all_obj[,! (all_obj@meta.data$sample %in% c("C51_Flush"))]

# Final Manual Annotation:

subcluster <- Sys.glob("/cluster/projects/macparland/TA/LiverMap2.0/RawData/Analysis/EmptyDrops/Subcluster/ManualAnnotation/*fullmeta*.rds")


cell_ID_celltype <- all_obj@meta.data[,c("cell_ID","sample", "Coarse_clusters", "Coarse_Manual_Anno")]

cell_ID_celltype[,"Fine_Manual_Anno"] <- as.character(cell_ID_celltype[,"Coarse_Manual_Anno"])

# Hepatocyte1 - using original coarse cluster
cell_ID_celltype[cell_ID_celltype[,2] == "0","Fine_Manual_Anno"] <- "PortalHep"
cell_ID_celltype[cell_ID_celltype[,2] == "1","Fine_Manual_Anno"] <- "PortalHep"
cell_ID_celltype[cell_ID_celltype[,2] == "4","Fine_Manual_Anno"] <- "InterzonalHep"
cell_ID_celltype[cell_ID_celltype[,2] == "5","Fine_Manual_Anno"] <- "CentralHep"
cell_ID_celltype[cell_ID_celltype[,2] == "9","Fine_Manual_Anno"] <- "CentralHep"
cell_ID_celltype[cell_ID_celltype[,2] == "19","Fine_Manual_Anno"] <- "ProlifInterzonalHep"

for( f in subcluster) {
        meta.data <- readRDS(f)
        tag <- strsplit(f, "_")[[1]][1]
	tag <- strsplit(tag, "/")[[1]]; tag <- tag[length(tag)]
        subset <- cell_ID_celltype[,1] %in% meta.data$cell_ID
        meta.data <- meta.data[match(cell_ID_celltype[,1], meta.data$cell_ID),]
        cell_ID_celltype[subset,"Fine_Manual_Anno"] <- paste(tag, meta.data[subset, "Subcluster_Manual"], sep="_")
}

#### Fine Annotation #####
 table(cell_ID_celltype$sample, cell_ID_celltype$Fine_Manual_Anno)

# Filter - Remove Doublets / Contamination
cell_ID_celltype <- cell_ID_celltype[! cell_ID_celltype$Fine_Manual_Anno %in% c("AntiB_Doublet", "Endo_HepContam", "Hepatocyte2_Female", "Macrophage_Debris", "Macrophage_Doublet", "NKT_Hepato", "NKT_HepDoublets", "Stellate_Contam"), ]

# Filter - Merge rare cell-types under general annotation
tmp <- table(cell_ID_celltype$sample, cell_ID_celltype$Fine_Manual_Anno)
# Ideally at least 10 cells in each included sample for each cell-type
apply(tmp, 2, function(x){sum(x > 5)}) < 10

# Should I exclude all proliferating cell types? Can't exclude just a few or may bias results. 
	# I think excluding, CC genes from markers / profiles will be better

cell_ID_celltype <- cell_ID_celltype[! cell_ID_celltype$Fine_Manual_Anno %in% c("NKT_Flush", "Endo_VasEndo"),] 
cell_ID_celltype[ cell_ID_celltype$Fine_Manual_Anno %in% 
			c("Stellate_Fibro-aHSC", "Stellate_Fibro-qHSC", "Stellate_aHSC", "Stellate_AP1+", "Stellate_qHSC"), 
			"Fine_Manual_Anno"] <- "Stellate_general"
cell_ID_celltype[ cell_ID_celltype$Fine_Manual_Anno %in% 
			c("NKT_Erythroblasts", "Hepatocyte2_Erythroblasts", "Hepatocyte2_Erythroid"), "Fine_Manual_Anno"] <- "Erythroid"
cell_ID_celltype[ cell_ID_celltype$Fine_Manual_Anno %in% 
			c("AntiB_IgK+IgA+", "AntiB_IgK+IgG+", "AntiB_IgL+IgA+", "AntiB_IgL+IgG+", "AntiB_Naive"), 
			"Fine_Manual_Anno"] <- "AntiB_general"

tmp <- table(cell_ID_celltype$sample, cell_ID_celltype$Fine_Manual_Anno)

# Filter Samples with most cell-types
keep <- apply(tmp, 1, quantile, probs=0.1) > 2

cell_ID_celltype <- cell_ID_celltype[cell_ID_celltype$sample %in% names(keep)[keep],]

tmp <- table(cell_ID_celltype$sample, cell_ID_celltype$Fine_Manual_Anno)
# Downsample #

max_cells_per_type = 1000

for (i in which(colSums(tmp) > max_cells_per_type)) {
	cells <- which(cell_ID_celltype$Fine_Manual_Anno == colnames(tmp)[i])
	keep <- sample(cells, max_cells_per_type)
	cell_ID_celltype <- cell_ID_celltype[-cells[!cells %in% keep],]
}

tmp <- table(cell_ID_celltype$sample, cell_ID_celltype$Fine_Manual_Anno)

subset <- all_obj[,all_obj$cell_ID %in% cell_ID_celltype$cell_ID]
identical(subset@meta.data$cell_ID, cell_ID_celltype$cell_ID)

saveRDS(cell_ID_celltype, file="for_MuSiC_subset_cells.rds")

require(Biobase)
for_MuSiC <- ExpressionSet(as.matrix(subset@assays$RNA@counts), phenoData=AnnotatedDataFrame(cell_ID_celltype))

#Key genes 
markers <- readRDS("~/Annotation_Package/markers.rds")
hvgs <- VariableFeatures(all_obj)

for_MuSiC@featureData <-AnnotatedDataFrame(data.frame(is.hvg = rownames(for_MuSiC) %in% hvgs, 
				is.marker = rownames(for_MuSiC) %in% markers$best[,2]))

saveRDS(for_MuSiC, file="for_MuSiC_ExpressionSet.rds")

