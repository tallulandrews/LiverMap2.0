### Currently only supports Human Data
# Automatically generates various plots to help with manual annotation of clusters.
# Includes: 
#  Dotplot / heatmaps of key marker genes
#  TCR and BCR constant chain expression
#  Correlation of hepatocytes with Halpern spatial stratification reference
#  Correlation of T-cells with Zheng sorted reference



args <- commandArgs(trailingOnly=TRUE)
# RDS file
# Cluster column.
# Coarse Annotation column (from automatic?)
# file_prefix

require("Seurat")

obj <- readRDS(args[1])

clusters <- obj@meta.data[,args[2]]

# For subsetting & specific analyses.
autoanno <- obj@meta.data[,args[3]]

prefix <- args[4]

coarse_anno <- as.character(autoanno)
coarse_anno[grep("Hep", coarse_anno)] <- "Hep"
coarse_anno[grep("NK*cell", coarse_anno)] <- "Lympho"
coarse_anno[grep("T*cell", coarse_anno)] <- "Lympho"
coarse_anno[grep("B*cell", coarse_anno)] <- "Lympho"

tab <- table(clusters, coarse_anno); tab <- tab/rowSums(tab);
Lympho_clusters <- rownames(tab)[tab[,"Lympho"] > 0.5]
Hepato_clusters <- rownames(tab)[tab[,"Hep"] > 0.5]

source("~/scripts/LiverMap2.0/My_R_Scripts.R")

#norm_mat <- obj@assays[[Assays(obj)[1]]]@data
#norm_mat <- obj@assays[[1]]@data
gene_cluster_means <- group_rowmeans(obj@assays[[1]]@data, clusters)


#scale_mat <- obj@assays[[1]]@scale.data
scale_cluster_means <- group_rowmeans(obj@assays[[1]]@scale.data, clusters)


# Key Markers - Immune cells 
# Source: https://www.nature.com/articles/ncomms14049
# Zheng et al. (2017) Massively parallel digita transcriptional profiling of single cells

zheng_immune_markers <- rbind(
		c("CD3D", "Tcells"),
		c("CD8A", "cytotoxic-T"),
		c("CD8B", ""),
		c("NKG7", "NKcells"),
		c("CD56", "NKcells"),
		c("GNLY", "NKcells"),
		c("FCER1A", "dendritic"),
		c("CLEC9A", "dendritic"),
		c("CLEC4C", ""),
		c("CD16", "CD16monocyte"),
		c("S100A8", "myeloid"),
		c("S100A9", "myeloid"),
		c("CD79A", "Bcells"),
		c("CD19", "Bcells"),
		c("CCR10", "Tmem"),
		c("TNFRSF18", "Treg"),
		c("CD25", "Treg"),
		c("ID3", "naive"),
		c("PF4", "megakaryocyte"),
		c("PTCRA", ""),
		c("CD45", ""),
		c("CD45RA", ""),
		c("CD45RO", ""),
		c("LGALS3", ""),
		c("SIGLEC7", ""),
		c("GZMK", "")
		)

Immune_Profiles <- readRDS("/cluster/projects/macparland/TA/ExternalData/Zheng_Immune/Zheng_profiles.rds");


####### Halpern Spatial Stuff
Halpern_data <- read.delim("/cluster/projects/macparland/TA/ExternalData/HalpernSpatialLocHep/NIHMS70855-supplement-Supplementary_Table_3.csv", sep=",", header=T)

# Remap gene IDs
gene_ids <- Halpern_data[,1]
Hepatocyte_spatial_profiles <- Halpern_data[,2:10]
exclude <- Halpern_data$q.values > 0.05 | is.na(Halpern_data$q.values);
Hepatocyte_spatial_profiles <- Hepatocyte_spatial_profiles[!exclude,];
gene_ids <- gene_ids[!exclude];
Hepatocyte_spatial_profiles <- t(apply(Hepatocyte_spatial_profiles, 1, scale));
colnames(Hepatocyte_spatial_profiles) <- paste("Layer", 1:9, sep="")


id_list <- strsplit( as.character(gene_ids), ";")

n_ids <- unlist(lapply(id_list, length))
max_n_ids <- max(unlist(lapply(id_list, length)))

source("~/R-Scripts/Ensembl_Stuff.R")

map_id <- function(ids) {
	h_ids <- General_Map(ids, in.org="Mmus", out.org="Hsap", in.name="symbol", out.name="symbol")
	h_ids <- h_ids[h_ids != ""]
	if (length(h_ids) == 0) {return("")}
	if (length(h_ids) == 1) {return(h_ids)}
	tab <- table(h_ids)
	h_ids <- names(tab)[which(tab ==max(tab))]
	return(h_ids[1]);
}

h_genes_simple <- General_Map(unlist(id_list[n_ids==1]), in.org="Mmus", out.org="Hsap", in.name="symbol", out.name="symbol")
h_genes_complex <- sapply(id_list[n_ids > 1], map_id)

h_genes <- as.character(gene_ids);
h_genes[n_ids==1] <- h_genes_simple;
h_genes[n_ids>1] <- h_genes_complex;

# deal with duplicated genes
q.values <- Halpern_data$q.values[!exclude]
dups <- unique(h_genes[duplicated(h_genes)])
dups <- dups[dups != ""]

for (g in dups) {
	qs <- q.values[h_genes==g];
	best <- which(qs == min(qs));
	best <- best[1]
	keep <- which(h_genes==g)[best]
	h_genes[h_genes == g] <- ""
	h_genes[keep] <- g
}

Hepatocyte_spatial_profiles <- Hepatocyte_spatial_profiles[h_genes != "",]
h_genes <- h_genes[h_genes != ""]
rownames(Hepatocyte_spatial_profiles) <- h_genes
saveRDS(Hepatocyte_spatial_profiles, file="/cluster/projects/macparland/TA/ExternalData/HalpernSpatialLocHep/Human_ortho_sig_profiles.rds");

Hepatocyte_spatial_profiles <- readRDS("/cluster/projects/macparland/TA/ExternalData/HalpernSpatialLocHep/Human_ortho_sig_profiles.rds");
####### Halpern Done

# Read in Marker Gene Tables
immune <- read.table("/cluster/projects/macparland/TA/ExternalData/My_Markers/Immune_markers.csv", sep =",", header=T)
liver <- read.table("/cluster/projects/macparland/TA/ExternalData/My_Markers/Liver_markers.csv", sep=",", header=T)
soupX <- read.table("/cluster/projects/macparland/TA/ExternalData/My_Markers/SoupXGenesets_markers.csv", sep=",", header=T)


# Read in Map 1 marker genes

source("~/scripts/LiverMap2.0/Setup_autoannotation.R")







# Hepatocyte Annotation
sync_profiles <- function(to_anno_means, ref_profiles) {
	common_genes <- sort(rownames(to_anno_means)[rownames(to_anno_means) %in% rownames(ref_profiles)])
	to_anno_means <- to_anno_means[match(common_genes, rownames(to_anno_means)),]
	ref_profiles <- ref_profiles[match(common_genes, rownames(ref_profiles)),]

	tmp <- colnames(to_anno_means)
	to_anno_means <- t(apply(to_anno_means, 1, scale));
	colnames(to_anno_means) <- tmp

	tmp <- colnames(ref_profiles)
	ref_profiles <- t(apply(ref_profiles, 1, scale));
	colnames(ref_profiles) <- tmp
	return(list(anno=to_anno_means, ref=ref_profiles));
}


#hep_profiles <- scale_cluster_means[,Hepato_clusters]
#common_genes <- sort(rownames(hep_profiles)[rownames(hep_profiles) %in% rownames(Hepatocyte_spatial_profiles)])

#hep_profiles <- hep_profiles[match(common_genes, rownames(hep_profiles)),]
#hep_ref <- Hepatocyte_spatial_profiles[match(common_genes, rownames(Hepatocyte_spatial_profiles)),]
#tmp <- colnames(hep_ref)
#hep_ref <- t(apply(hep_ref, 1, scale))
#colnames(hep_ref) <- tmp;

#cluster_profiles <- t(apply(hep_profiles, 1, scale))
#colnames(cluster_profiles) <- colnames(hep_profiles)

synced <- sync_profiles(scale_cluster_means[,Hepato_clusters], Hepatocyte_spatial_profiles);
cluster_profiles <- synced$anno
hep_ref <- synced$ref;

require(proxy)
cors <- simil(t(cluster_profiles), t(hep_ref))
cos <- simil(t(cluster_profiles), t(hep_ref), method="cosine")

require(gplots)
require(RColorBrewer)
heat_col <- rev(c(rev(brewer.pal(4, "Reds")), "white", brewer.pal(4,"Blues")))

png(paste(prefix, "Hep_vs_Halpern_cosine.png", sep="_"), width=8, height=8, units="in", res=300)
heatmap.2(cos, trace="none", Rowv=TRUE, Colv=FALSE, dendrogram=c("none"), symm=TRUE, scale="none", symbreaks=TRUE, col=heat_col)
dev.off()
png(paste(prefix, "Hep_vs_Halpern_cors.png", sep="_"), width=8, height=8, units="in", res=300)
heatmap.2(cors, trace="none", Rowv=TRUE, Colv=FALSE, dendrogram=c("none"), symm=TRUE, scale="none", symbreaks=TRUE, col=heat_col)
dev.off()







# Immune Annotation

immune_ref <- Immune_Profiles$profiles[Immune_Profiles$specificity > 5,]

zheng_unique_ness <- apply(Immune_Profiles$profiles, 1, function(x) {
				fst <- max(x); 
				snd <- max(x[x < max(x)])
				(fst - snd)/snd
				})
zheng_unique_ness[is.na(zheng_unique_ness)] <- 0

#tmp <- unlist(apply(Immune_Profiles$profiles[Immune_Profiles$specificity > 5,], 1, function(x){if (max(x) > 1) {which(x==max(x))} else{return(0)}}))
tmp <- unlist(apply(Immune_Profiles$profiles[zheng_unique_ness>0.5,], 1, function(x){if (max(x) > 1) {which(x==max(x))} else{return(0)}}))

zheng_empiric_markers <- rbind(
		c("MZT2A","t_memory"), c("IL7A","t_memory"),
		c("HINT1","t_memory"), c("NDUFB9","t_memory"),
		c("CORO1B","t_memory"),	c("TRADD","t_memory"),
		c("MIF","t_memory"), c("CD8C","t_cyto_naive"),
		c("SNHG8","t_cyto_naive"), c("HSPB1","t_cyto_naive"),
		c("FBL","t_cyto_naive"), c("S100B","t_cyto_naive"),
		c("SELL","t_naive"), c("MAL","t_naive"),
		c("LEF1","t_naive"), c("CCR7","t_naive"),
		c("CALM3","t_naive"), c("PIK3IP1","t_naive"),
		c("PBXIP1","t_reg"), c("AQP3","t_reg"),
		c("CD44","t_reg"), c("UCP2","t_reg"),
		c("IL10RA","t_reg"), c("S1PR4","t_reg"),
		c("GIMAP5","t_helper"), c("CMPK1","t_helper"),
		c("SSR2","t_helper"), c("EIF3F","t_helper"),
		c("EEF2","t_helper"), c("CD3G","t_cyto"),
		c("DUSP2","t_cyto"), c("CXCR4","t_cyto"),
		c("CD79B","bcells"), c("CD79A","bcells"),
		c("BANK1","bcells"), c("BLK","bcells"),
		c("SPIB","bcells"), c("VPREB3","bcells"),
		c("FCER2","bcells"), c("TCL1A","bcells"),
		c("FCGR3A","nkcells"), c("XCL2","nkcells"),
		c("SPON2","nkcells"), c("FGFBP2","nkcells"),
		c("CLIC3","nkcells"), c("PRF1","nkcells"),
		c("GZMB","nkcells"), c("CCL4","nkcells"),
		c("S100A9","monocyte"), c("S100A8","monocyte"),
		c("MNDA","monocyte"), c("CD14","monocyte"),
		c("LST1","monocyte"), c("CPVL","monocyte"),
		c("CFP","monocyte"), c("FCN1","monocyte"),
		c("LYZ","monocyte"), c("SERPINA1","monocyte"),
		c("CD68","monocyte"), c("LGALS2","monocyte"),
		c("TYMP","monocyte"), c("CST3","monocyte")
		)

if (length(Lympho_clusters) > 1) {
	synced <- sync_profiles(scale_cluster_means[,Lympho_clusters], immune_ref)
	require(proxy)
	cors <- simil(t(synced$anno), t(synced$ref))
	cos <- simil(t(synced$anno), t(synced$ref), method="cosine")

	require(gplots)
	require(RColorBrewer)
	heat_col <- rev(c(rev(brewer.pal(4, "Reds")), "white", brewer.pal(4,"Blues")))

	png(paste(prefix, "Immune_vs_Zheng_cosine.png", sep="_"), width=8, height=8, units="in", res=300)
	heatmap.2(cos, trace="none", Rowv=TRUE, Colv=FALSE, dendrogram=c("none"), symm=TRUE, scale="none", symbreaks=TRUE, col=heat_col)
	dev.off()
	png(paste(prefix, "Immune_vs_Zheng_cors.png", sep="_"), width=8, height=8, units="in", res=300)
	heatmap.2(cors, trace="none", Rowv=TRUE, Colv=FALSE, dendrogram=c("none"), symm=TRUE, scale="none", symbreaks=TRUE, col=heat_col)
	dev.off()
}

#Immune + scamp

require("scmap")
require("SingleCellExperiment")
tmp <- colnames(immune_ref);
immune_scmap <- t(apply(immune_ref, 1, scale))
colnames(immune_scmap) <- tmp;

immune_cells <- obj[,clusters %in% Lympho_clusters]
immune_cells_sce <- as.SingleCellExperiment(immune_cells)
immune_cells_sce <- immune_cells_sce[rownames(immune_cells_sce) %in% rownames(immune_cells@assays[[1]]@scale.data),]

immune_cells_sce <- immune_cells_sce[match(rownames(immune_cells@assays[[1]]@scale.data), rownames(immune_cells_sce)),]

assays(immune_cells_sce)[["logcounts"]] <- immune_cells@assays[[1]]@scale.data;
rowData(immune_cells_sce)$feature_symbol <- rownames(immune_cells@assays[[1]]@scale.data);

cell_level <- scmapCluster(projection=immune_cells_sce, index_list = list(zheng=immune_scmap), threshold=0.1)

cell_immune_lab <- cell_level$combined_labs;
names(cell_immune_lab) <- colnames(immune_cells_sce)

obj@meta.data$Immune_scmap <- cell_immune_lab[match(rownames(obj@meta.data), names(cell_immune_lab))]
obj@meta.data$Immune_scmap[is.na(obj@meta.data$Immune_scmap)] <- "Non-Lympho"

png(paste(prefix, "scmap_immune.png", sep="_"), width=9, height=7, units="in", res=300)
DimPlot(obj, reduction="umap", group.by="Immune_scmap", label=TRUE, pt.size=0.5)
dev.off()



#Immune marker gene plots.

zheng_immune_markers=zheng_immune_markers[zheng_immune_markers[,1] %in% rownames(obj),]
zheng_empiric_markers=zheng_empiric_markers[zheng_empiric_markers[,1] %in% rownames(obj),]

png(paste(prefix, "zheng_markers.png", sep="_"), width=8, height=8, units="in", res=300)
Seurat::DotPlot(obj, features=zheng_immune_markers[,1], group.by=args[2])
dev.off()
png(paste(prefix, "zheng_emp_markers.png", sep="_") width=8, height=8, units="in", res=300)
Seurat::DotPlot(obj,features=zheng_empiric_markers[,1], group.by=args[2])
dev.off()

# Key Markers
# From Spatial - Central vs Portal Hepatocytes
Central <- c("CYP1A2", "CYP2E1", "CYP3A4", "GLUL", "DCXR", "FTL", "GPX2", "GSTA1")
Portal <- c("CYP2A7", "FABP1", "HAL", "AGT", "ALDOB", "SDS")
# Cross reference Map1 with marker tables
Stellate <- c("ACTA2", "COL1A1", "RBP1", "TAGLN", "ADAMTSL2", "GEM", "LOXL1", "LUM")
pLSEC <- c("GJA5", "SPARCL1", "CLEC14A", "PLVAP", "EGR3")
cvLSEC <- c("FCN2", "CLEC1B", "CLEC4G", "PVALB", "S100A13")
Eryth <- c("HBB", "HBA1", "HBA2")
AntiBcell <- c("IGKC", "JCHAIN", "IGHA1", "IGLC1", "IGLC2", "IGLC3")
Bcell <- c("CD22", "CD37", "CD79B", "FCRL1", "LTB", "DERL3", "IGHG4")
NonInfMac <- c("VCAM1", "TTYH3", "TIMD4", "SLC40A1", "RAB31", "MARCO", "HMOX1", "C1QC")
InfMac <- c("VCAN", "S100A8", "MNDA", "LYZ", "FCN1", "CXCL8")
CD3abTcell <- c("CD8A", "CD8B", "CD3D", "CD3G", "TRAC", "IL32", "TRBC1", "TRBC2")
gdTcell <- c("CSTW", "IL7R", "GZMB", "GZMH", "TBX21", "HOPX", "PRF1", "S100B", "TRDC", "TRGC1", "TRGC2")
NKcell <- c("NCR1", "NCR2", "NCR3", "NKG7")
endo <- c("VWF", "RAMP3", "LIFR", "MYL9", "ACTA2", "CLDN4", "ANXA4", "COL1A2")


best_markers <- as.character(Central, Portal, Stellate, pLSEC, cvLSEC, endo, Eryth, AntiBcell, Bcell, NonInfMac, InfMac, CD3abTcell, gdTcell, NKcell)
immune_phenotype <- c(soupX[,2],"TNF", "JCHAIN", "IGKV1-12", "IGKV4-1", "IGLV3-1", "IGLV6-57", "IGLL5", "IGLC7", rownames(obj)[grep("CXC",rownames(obj))], rownames(obj)[grep("^IL",rownames(obj))], rownames(obj)[grep("^HLA",rownames(obj))])


png(paste(prefix, "SoupX_markers.png", sep="_"), width=8, height=8, units="in", res=300)
Seurat::DotPlot(obj, features=soupX[,2], group.by=args[2])
dev.off()
png(paste(prefix, "Immune_pheno.png", sep="_"), width=8, height=8, units="in", res=300)
Seurat::DotPlot(immune_cells, features=immune_pheontype, group.by=args[2])
dev.off()

png(paste(prefix, "best_markers.png", sep="_"), width=8, height=8, units="in", res=300)
Seurat::DotPlot(obj, features=best_markers, group.by=args[2])
dev.off()
### General Marker Genes.


