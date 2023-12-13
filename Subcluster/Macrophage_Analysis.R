source("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/IndividualVariability.R")
source("C:/Users/tandrews/Documents/UHNSonya/scripts/LiverMap2.0/Colour_Scheme.R")
source("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/SubColour_Scheme.R")
source("C:/Users/tandrews/Documents/UHNSonya/scripts/LiverMap2.0/My_R_Scripts.R")

require(fgsea)
require(Seurat)
require(ggplot2)

immune_path <- gmtPathways("../../ExternalData/MSigDb_immune_signatures_c7.all.v7.1.symbols.gmt")
Hallmark_path <- gmtPathways("../../ExternalData/MSigDbHalmarkPathways.gmt")
MSigAll <- gmtPathways("../../ExternalData/MSigDb_curated_c2.all.v7.1.symbols.gmt")
reactome <- gmtPathways("../../ExternalData/ReactomePathways.gmt")
BaderMSig <- gmtPathways("../../ExternalData/BaderLab25Aug2020/Human_MSigdb_August_01_2020_symbol.gmt.txt")
BaderReact <- gmtPathways("../../ExternalData/BaderLab25Aug2020/Human_Reactome_August_01_2020_symbol.gmt.txt")

flow_cyto_genes <- c("CD68", "MRC1", "CD14", "CD274", "PTPRC", "CD86", "CD163", "IL2RA", "FOXP3",
				"NCAM1", "NCR1", "CD3E", "PTPRC", "CD8A", "IL7R", "ICOS", "KLRD1", "IL2RA",
				"CD4", "PDCD1", "SELL", "HAVCR2", "CCR3", "CD3E", "IL7R", "CD8A", "PTPRC", 
				"LAG3", "CTLA4", "CD27", "CD4")


label_genes <- function(g, name) { names(g)<- rep(name, length(g)); return(g)}

#Removed Genes: "JUND", "FOS", "NFKBIA", "ACTG1", "CD14", 

Macrophage_genes_dot <- c(
		label_genes(c("MARCO", "CD5L", "SLC40A1", "FTL", "CD163", "SEPP1", "C1QC", "C1QB", "C1QA", "CTSB", "HMOX1", "VCAM1"), "NonInfammatory"), 
		label_genes(c("HLA-DRA", "HLA-DPA1", "HLA-DQB1", "HLA-DPB1", "CD74"), "MHCII"), "VSIG4", "NINJ1", "IL18", 
		label_genes(c("LYZ", "S100A4", "S100A6", "S100A8", "S100A9", "S100A12", "MNDA", "VCAN", "FCN1"), "Infammatory"),
		label_genes(c("ACP5", "PLD3", "PSAP", "CSTB", "LGMN"), "Lysosome"),"FTH1", "CD68", "APOE", "GLUL", 
		label_genes(c("FABP5", "GPNMB", "CD9", "SPP1", "APOC1"), "LAM"),
		"RBP7", "PLTP", label_genes(c("FOLR2", "TIMD4", "LYVE1", "FCER1G", "MS4A7", "TIMP1"),"Resident"), "CD14",
		label_genes(c("CXCL3", "THBS1", "NAMPT", "CXCL2", "CD83", "IL1B", "AREG", "CCL3","PLAUR", "SRGN"), "Activation"),
		"PLAC8", "CD54", label_genes(c("LST1", "IFITM3", "AIF1", "COTL1"), "Synapse"),
		label_genes(c("DNASE1L3", "FCN2", "CCL14", "FCN3", "SPARC", "CLEC1B", "ENG"), "LSECs")
		) #label_genes(c("ALB", "SERPINA1", "HP", "FGA"),"Hepto")

# from: https://www.frontiersin.org/articles/10.3389/fimmu.2019.01084/full
# Key M1/M2 and Classic vs Alternate 

Classical <- c("MS4A4A", "FPR2", "CXCL9", "LCN2", "CD38", "IL12B", "FPR1", "CXCL10", 
			"CFB", "IL1B", "IL1A", "IFIT1", "GBP6", "IFIT2", "GPR31B", "PTGES", 
			"OASL1", "RSAD2", "ISG20", "ZFP811", "ADGB", "FSCN1")
Alternate <- c("ARG1", "MGL2", "TMEM26", "RNASE2A", "MRC1", "EGR2", "FLT1", "CHIL3", "CLEC10A",
			"MATK", "SOCS2", "ITGB3", "OCSTAMP", "PTGS1", "S100a4", "CLEC7A", 
			"PLXDC2", "HBEGF", "CCL24", "EMP1", "OLFM1", "UBE2C")
M1_Only <- c("IRG1", "CCL2", "CSF2", "MARCKSL1", "IL23A", "CMPK2", "RTP4", "SLC7A11", 
			"MAFF", "GPR85", "PDE4B", "TREX1", "IER3", "IFI47", "NFKBIE", "OASL2",
			"GDF15")

# Dendritic Cell marker from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4403526/
#cDC1 <- CD45+("PTPRC") HLA-DR+("HLA-DRA") CD141+("THBD") CD123-("IL3RA")CD11c+("ITGAX") CD14-("CD14")
#cDC2 <- CD45+("PTPRC") HLA-DR+("HLA-DRA") CD1c+("CD1C") CD123-("IL3RA") CD11c+("ITGAX") CD141-("THBD") CD14+ ("CD14")
#pDC <- HLA-DR+("HLA-DRA") CD123+("IL3RA") CD11c-("ITGAX") CD303+("CLEC4C") CD304+ ("NRP1")
#Mac <- CD68+("CD68") CD163+("CD163")
#Monocytes <- CD14hi("CD14") CD16- ("FCGR3A")

# Sonya dendritic cells
# cDC2 <- CLEC10A + MHCII
# pDC <- CLEC4C + SCT
# cDC1 <- XCR1 + BATF3 + CLEC9A

#BioRad
cDC <- c("CD8A", "CD1C", "ITGAX", "ITGAM", "ITGAE", "LY75", "HLA-DRA")
pDC <- c("TLR7", "TLR9", "IL3RA", "LILRA4", "NRP1", "CLEC4C")

DotPlot(obj_all_genes, group.by="Subcluster_Manual", features=c("PTPRC", "HLA-DRA", "IL3RA", "ITGAX", "CD14", "THBD", "CLEC4C", "NRP1", "CD68", "CD163", "FCGR3A")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

obj_all_genes@reductions <- obj@reductions
FeaturePlot(obj_all_genes, features=c(cDC))
png("SupplFigure_DC_markers.png", width=6, height=6, units="in", res=150)
FeaturePlot(obj_all_genes, features=c("CLEC10A", "XCR1", "CLEC9A", "HLA-DRA"))
dev.off()
DotPlot(obj_all_genes, features=c(cDC, pDC), group.by="Subcluster_Manual")


suppl_tab1 <- get_cluster_summary(obj, samples="sample")
suppl_tab2 <- get_cluster_summary(obj, samples="donor")
write.table(suppl_tab2, file="Macrophages_cluster_summary.csv", sep=",")



extract_pathway_dat <- function(res, this_pathway, this_name) {
	genes <- unlist(res$rich[unlist(res$rich[,1]) == this_pathway,8])
	names(genes) <- rep(this_name, length(genes));

	bar_point <- unlist(abs(log10(res$rich[unlist(res$rich[,1]) == this_pathway,"padj"]))); names(bar_point)<- this_name;
	return(list(genes=genes, point=bar_point))
}




## Read in data ##
require(Seurat)
obj <- readRDS("Macrophage_varimax_Subcluster.rds")
cluster_col = "Coarse_clusters"
obj_all_genes <- readRDS("AllGenes/Macrophage_harmony_Subcluster_Allgenes.rds")


manual_cluster_anno <- c(
		"Synapse", "Kupffer", "Kupffer", "Monocyte", "Resident", "Mono-Act", 
		"Monocyte", "MHCII/DC", "Activated", "MHCII/DC", "LAM-like", 
		"Debris", "Monocyte", "LSEC-Doublet", "Kupffer")
# Changed 10 June 2022
#manual_cluster_anno <- c(
#		"InfSynap", "NonInf", "NonInf", "Inflam", "ResNonInf", "InflamActiv", 
#		"Inflam", "MHCII", "Activated", "MHCII", "PhagoNonInf", 
#		"Debris", "Inflam", "Doublet", "NonInf")
# Changed 25 May 2021 so that Cluster 2 (third in list) is "NonInf" rather than "Debris", 
# and Cluster 8 is "Activated" instead of "Repair"


obj@meta.data$Subcluster_Manual <- manual_cluster_anno[obj@meta.data[,cluster_col]]
obj_full <- obj
obj_all_genes@meta.data <- obj@meta.data
obj_all_genes@reductions <- obj@reductions


saveRDS(obj@meta.data, "Macrophage_fullmetadata.rds")
# Annotation Dotplot
png("Macrophage_Supplementary_Dotplot.png", width=16, height=7, units="in", res=300)
DotPlot(obj_all_genes, group.by="Subcluster_Manual", features=c(Macrophage_genes_dot)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# Marker Heatmap
tmp <- obj_all_genes[rownames(obj_all_genes) %in% Macrophage_genes_dot,]
tmp2 <- t(apply(tmp@assays$RNA@scale.data, 1, scale)); colnames(tmp2) <- colnames(tmp)
tmp@assays$RNA@scale.data <- tmp2;
png("Macrophage_Supplementary_heatmap.png", width=7, height=16, units="in", res=300)
DoHeatmap(tmp, group.by="Subcluster_Manual", features=Macrophage_genes_dot, slot="scale.data", assay="RNA") 
dev.off()



obj_all_genes@meta.data$Anno_cluster <- paste(obj_all_genes@meta.data$Subcluster_Manual, " (", obj_all_genes@meta.data$Coarse_clusters, ") ", sep="")
png("Macrophage_Supplementary_Anno_Dotplot.png", width=16, height=7, units="in", res=300)
DotPlot(obj_all_genes, group.by="Anno_cluster", features=Macrophage_genes_dot) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


DotPlot(obj_all_genes, group.by="Subcluster_Manual", features=M1_Only) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

subtype_order <- c("Monocyte", "Mono-Act", "Activated", "Synapse", "MHCII/DC", "LAM-like", "Resident", "Kupffer", "LSEC-Doublet", "Debris")
subtype_cols <- get_seurat_colours(obj, "Subcluster_Manual"); subtype_cols <- subtype_cols[c(4,5,1,6,7,9,10,8,3,2)];
names(subtype_cols) <- subtype_order
obj@meta.data$Subcluster_Manual <- factor(obj@meta.data$Subcluster_Manual, levels=subtype_order)
png("Figure_Macrophage_FinalSubcluster_UMAP.png", width=8*0.75, height=7*0.75, units="in", res=300)
DimPlot(obj, group.by="Subcluster_Manual", label=F, cols=subtype_cols)
dev.off()


### Object for Shiny ###
shiny_obj <- obj_all_genes;
detect_rate <- group_rowmeans(shiny_obj@assays$RNA@counts > 0, shiny_obj@meta.data$Subcluster_Manual)
exclude_genes <- apply(detect_rate, 1, max) < 0.005
shiny_obj@meta.data <- shiny_obj@meta.data[,c("nCount_RNA","nFeature_RNA","percent.mt","Phase", "donor", "sample", "donor_sex", "donor_age", "Coarse_clusters", "Subcluster_Manual")]
shiny_obj@reductions <- obj@reductions
shiny_obj@assays$RNA@scale.data <- matrix(0)
shiny_obj@assays$RNA@var.features <- c(0)
shiny_obj <- shiny_obj[!exclude_genes,]
expr_mat <- shiny_obj@assays$RNA@data
metadata <- shiny_obj@meta.data
metadata$UMAP1 <- shiny_obj@reductions$umap@cell.embeddings[,1]
metadata$UMAP2 <- shiny_obj@reductions$umap@cell.embeddings[,2]
metadata$cell_colour <- subtype_cols[match(shiny_obj@meta.data$Subcluster_Manual, names(subtype_cols))]


saveRDS(shiny_obj, "ShinyApps/Macrophage_Lattice_obj.rds")
saveRDS(list(expr_mat=expr_mat, metadata=metadata), "ShinyApps/Macrophage_shiny.rds")

# ----------------------------------- #





neutrophils <- DotPlot(obj_all_genes, features=c("FCGR3B", "CXCL8", "CXCR2", "CXCR1", "IFITM2", "CSF3R", "FPR1", "S100A11", "BASP1", "NAMPT", "G0S2"), group.by="Subcluster_Manual")
#Remmerie
lipid_ass_mac <- DotPlot(obj_all_genes, features=c("CLEC4F", "TIMD4", "CD5L", "ZEB2", "LDHA", "PGK1", "MIF", "ALDOA", "TPI1", "PGAM1"), group.by="Subcluster_Manual")
#Guilliams
lipid_ass_mac <- DotPlot(obj_all_genes, features=c("GPNMB", "SPP1", "CD63", "CD93", "FABP5", "TREM2", "FAM20C", "EMP1", "F7", "GDF15", "ITGA6"), group.by="Subcluster_Manual")

#Consensus
lipid_ass_mac <- DotPlot(obj_all_genes, features=c("GPNMB", "TREM2",  "CD93", "CD63", "FABP5", "MS4A7", "MMP14",  "PTAFR", "FCGR2B", "SPP1"), group.by="Subcluster_Manual")
png("Macrophage_Supplementary_LAM_Dotplot.png", width=16*0.4, height=7*0.5, units="in", res=300)
print(lipid_ass_mac+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()

## Dendritic cells for Paper ##
 obj_all_genes@reductions$umap <- obj@reductions$umap
png("Macrophage_DC_markers.png", width=6, height=4, units="in", res=300)
FeaturePlot(obj_all_genes, features=c("CLEC10A", "XCR1", "CLEC9A"))
dev.off()


### ------------ TLM -------------###

require(Seurat)
obj <- readRDS("Macrophage_varimax_Subcluster.rds")
cluster_col = "Coarse_clusters"
obj_all_genes <- readRDS("AllGenes/Macrophage_harmony_Subcluster_Allgenes.rds")

manual_cluster_anno <- c(
		"Monocyte", "Kupffer", "Kupffer", "Monocyte", "Kupffer", "Monocyte-Act", 
		"Monocyte", "MHCII/DC", "Activated", "MHCII/DC", "LAM-like", 
		"Debris", "Monocyte", "LSEC-Doublet", "Kupffer")



obj@meta.data$Subcluster_Manual <- manual_cluster_anno[obj@meta.data[,cluster_col]]
obj_all_genes@meta.data <- obj@meta.data

obj_all_genes <- obj_all_genes[,! obj_all_genes@meta.data$Subcluster_Manual %in% c("Debris", "LSEC-Doublet")]

Macrophage_genes_TLM <- c(
		label_genes(c("MARCO", "CD5L", "SLC40A1", "CD163", "C1QC", "C1QB", "C1QA", "HMOX1", "VCAM1"), "Kupffer"), 
		label_genes(c("LYZ", "S100A4", "S100A6", "S100A8", "S100A9", "S100A12", "MNDA", "VCAN", "FCN1"), "Monocyte"),
		label_genes(c("FABP5", "ACP5", "PLD3", "LGMN", "GPNMB", "TREM2", "CD9", "CD68", "SPP1"), "LAM"),
		label_genes(c("FOLR2", "TIMD4", "LYVE1", "FCER1G", "MS4A7", "TIMP1"),"Resident"),
		label_genes(c("CXCL3", "THBS1", "NAMPT", "CXCL2", "CD83", "IL1B", "AREG", "CCL3","PLAUR", "SRGN"), "Activated")
)

png("Macrophage_TLM_Dotplot.png", width=16, height=7, units="in", res=300)
DotPlot(obj_all_genes, group.by="Subcluster_Manual", features=Macrophage_genes_TLM) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



### ------ AUCell vs Guilliams ------ ###
genesets <- readRDS("../../ExternalData/GuilliamsCellPaper/Guilliams_Genesets.rds")

require(AUCell)
cells_rankings <- AUCell_buildRankings(obj_all_genes@assays$RNA@scale.data)
cells_AUC <- AUCell_calcAUC(genesets, cells_rankings, aucMaxRank=5, nCores=1)

tmp <- cells_AUC@assays@data$AUC
tmp[tmp < 0] <- 0

assignment<-apply(tmp, 2, function(x) { 

							out <- which(x==max(x)); 
							if (length(out) > 1) {return("None")} 
							else {return (rownames(cells_AUC@assays@data$AUC)[out])}
							}
			)
obj@meta.data$AUCell <- assignment
DimPlot(obj, split.by="AUCell", group.by="Subcluster_Manual")

obj@meta.data$cDC1 <- cells_AUC@assays@data$AUC["HumancDC1",]
obj@meta.data$cDC2 <- cells_AUC@assays@data$AUC["HumancDC2",]
obj@meta.data$Mono <- cells_AUC@assays@data$AUC["HumanMono",]
obj@meta.data$KC <- cells_AUC@assays@data$AUC["HumanKC",]
obj@meta.data$imLAM <- cells_AUC@assays@data$AUC["HumanLAM",]
obj@meta.data$matLAM <- cells_AUC@assays@data$AUC["HumanMatLAM",]

png("Macrophage_Guilliams_AUCell.png", width=6, height=8, units="in", res=300)
FeaturePlot(obj, feature=c("cDC1", "cDC2", "Mono", "KC", "imLAM", "matLAM"))
dev.off()

## ------------ ##


# GLUL Plot
png("Macrophage_GLUL_umap.png", width=4, height=4, units="in", res=300)
FeaturePlot(obj, features="GLUL")
dev.off()
png("Macrophage_GLUL_dotplot.png", width=4, height=4, units="in", res=300)
DotPlot(obj_all_genes, group.by="Subcluster_Manual", features="GLUL") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

## DE # For Lewis
cell_types <- obj_full@meta.data$Subcluster_Manual
#cell_types[cell_types=="InflamActiv"] <- "Inflam"
#cell_types[cell_types=="InfSynap"] <- "Inflam"
#cell_types[cell_types=="PhagoNonInf"] <- "NonInf"
#cell_types[cell_types=="ResNonInf"] <- "NonInf"

pseudo_detect <- group_rowmeans(obj_all_genes@assays$RNA@counts > 0, cell_types)
pseudo_mean <- group_rowmeans(obj_all_genes@assays$RNA@data, cell_types)
is.clean <- ! colnames(pseudo_detect) %in% c("Debris", "Dublet")

de_synap <- FindMarkers(obj_all_genes, group.by="Subcluster_Manual", ident.1="InfSynap", ident.2="Inflam",  logfc.threshold=-Inf)

for (i in colnames(pseudo_detect[,is.clean])) {
#	de <- FindMarkers(obj_all_genes, group.by="Subcluster_Manual", ident.1=i,  logfc.threshold=-Inf)
#	write.table(de, file=paste("SubCluster_DE_Tables/Macrophage_", i, "_DE.csv", sep=""), sep=",", row.names=T, col.names=T, quote=T)
	de <- read.table(file=paste("SubCluster_DE_Tables/Macrophage_", i, "_DE.csv", sep=""), sep=",", header=TRUE)

	l2fc <- log2(pseudo_mean[,i]/ (apply(pseudo_mean[,is.clean & colnames(pseudo_detect) != i], 1, max)) )
	detect <- pseudo_detect[,i]- apply(pseudo_detect[,is.clean & colnames(pseudo_detect) != i], 1, max)

	write.table(de[rownames(de) %in% names(detect)[detect > 0.1],], paste("ForLewis_All", i, "Markers.csv", sep="_"), sep=",", row.names=T, col.names=T)

}



## Macrophage Endothelial ##
# Significance of doublets.
doublet_rate = 0.05
LSEC_global_freq = 0.06
obs_doublets = sum(obj@meta.data$Subcluster_Manual == "Doublet")
pbinom(obs_doublets-1, size = ncol(obj), prob=doublet_rate*LSEC_global_freq, lower.tail=F)

# Which Doublets? 
obj_endo <- readRDS("Endo_varimax_Subcluster.rds")
obj_endo@meta.data <- readRDS("Endo_fullmetadata.rds")
obj_endo_all_genes <- readRDS("AllGenes/Endo_harmony_Subcluster_Allgenes.rds")
obj_endo_all_genes@meta.data <- obj_endo@meta.data


# SCMAP
set.seed(101)
require(scmap)
require(SingleCellExperiment)
doublet_obj <- obj_all_genes[,obj_all_genes@meta.data$Subcluster_Manual == "LSEC-Doublet"]
doublet_obj_sce <- Seurat::as.SingleCellExperiment(doublet_obj)

# Endothelial Reference
Endo_genes_dot <- c(label_genes(c("FCN2", "FCN3", "STAB1", "STAB2", "CLEC1B", "CLEC4G", "CLEC4M", "CRHBP"), "cvLSEC"),
			label_genes(c("SPARCL1", "TM4SF1", "CLEC14A", "ID1", "IGFBP7", "MGP", "ADRIF", "S100A6", "AQP1"), "ppLSEC"), "VWF", 
			label_genes(c("ENG", "PECAM1",  "DNASE1L3", "TIMP3", "LIFR", "C7"), "Endo"), "RAMP3",
			label_genes(c("RSPO3", "ACKR1", "WNT2"), "cvEndo"), "INMT", "PLAC8",
			label_genes(c("PODXL", "PLVAP", "CD34"), "Arterial"), "GSN", "RBP7",
			label_genes(c("CCL21", "CST3", "ALB", "APOA1", "PDPN"), "Hepato"))
endo_markers <- Endo_genes_dot[ !(names(Endo_genes_dot) %in% c("Hepato"))]
endo_profiles <- group_rowmeans(obj_endo_all_genes@assays$RNA@data, factor(obj_endo_all_genes@meta.data[,"Subcluster_Manual"]))
endo_scmap_index <- endo_profiles[rownames(endo_profiles) %in% endo_markers,]
endo_scmap_index <- endo_scmap_index[, colnames(endo_scmap_index) != "HepContam"]

# Macrophage Reference
macrophage_markers <- Macrophage_genes_dot[! (names(Macrophage_genes_dot) %in% c("Debris", "LSEC-Doublet"))]
macrophage_profiles <- group_rowmeans(obj_all_genes@assays$RNA@data, factor(obj_all_genes@meta.data[,"Subcluster_Manual"]))
macrophage_scmap_index <- macrophage_profiles[rownames(macrophage_profiles) %in% macrophage_markers,]
macrophage_scmap_index <- macrophage_scmap_index[, !(colnames(macrophage_scmap_index) %in% c("Debris", "LSEC-Doublet"))]

rowData(doublet_obj_sce)$feature_symbol <- rownames(doublet_obj_sce)
scmapCluster(doublet_obj_sce, index_list = list(mac=macrophage_scmap_index, endo=endo_scmap_index))

doublet_expr_mac <- doublet_obj@assays$RNA@data[match(rownames(macrophage_scmap_index), rownames(doublet_obj)),]
cor_mac1 <- cor(as.matrix(doublet_expr_mac), macrophage_scmap_index)
cor_mac2 <- cor(as.matrix(doublet_expr_mac), macrophage_scmap_index, method="spearman")
hits1 <- t(apply(cor_mac1, 1, function(x) {c(max(x), which(x == max(x)))}))
hits2 <- t(apply(cor_mac2, 1, function(x) {c(max(x), which(x == max(x)))}))
mac_hits1 <- data.frame(hits1, type=colnames(cor_mac1)[hits1[,2]])
mac_hits2 <- data.frame(hits2, type=colnames(cor_mac1)[hits2[,2]])

doublet_expr_endo <- doublet_obj@assays$RNA@data[match(rownames(endo_scmap_index), rownames(doublet_obj)),]
cor_endo1 <- cor(as.matrix(doublet_expr_endo), endo_scmap_index)
cor_endo2 <- cor(as.matrix(doublet_expr_endo), endo_scmap_index, method="spearman")
hits1 <- t(apply(cor_endo1, 1, function(x) {c(max(x), which(x == max(x)))}))
hits2 <- t(apply(cor_endo2, 1, function(x) {c(max(x), which(x == max(x)))}))
endo_hits1 <- data.frame(hits1, type=colnames(cor_endo1)[hits1[,2]])
endo_hits2 <- data.frame(hits2, type=colnames(cor_endo1)[hits2[,2]])

sum(mac_hits1[,2] == mac_hits2[,2] & endo_hits1[,2] == endo_hits2[,2])
keep <- (mac_hits1[,2] == mac_hits2[,2] & endo_hits1[,2] == endo_hits2[,2])


mac_endo_results <- cbind(mac_hits1[keep,3], endo_hits1[keep,3])
heat_dat <- matrix(0, nrow=length(unique(mac_endo_results[,1])), ncol=length(unique(mac_endo_results[,2])))
rownames(heat_dat) <- sort(unique(mac_endo_results[,1]));
colnames(heat_dat) <- sort(unique(mac_endo_results[,2]));
for(i in rownames(heat_dat)) {
	for(j in colnames(heat_dat)) {
		heat_dat[i,j] <- sum(mac_endo_results[,1]==i & mac_endo_results[,2] == j)
	}
}
require(gplots)
png("Figure_Mac_Endo_Doublets_Scmap_freq.png", width=5, height=5, units="in", res=150)
heatmap.2(heat_dat, cellnote=heat_dat, trace="none", col=colorRampPalette(c("grey60", "black")), notecol="white", mar=c(10,10),cexRow=1.3, cexCol=1.3)
dev.off()


# Macrophage-Endo Communication Barplots #
mac_detect <- group_rowmeans(obj_all_genes@assays$RNA@data > 0, obj@meta.data$Subcluster_Manual)
endo_detect <- group_rowmeans(obj_endo_all_genes@assays$RNA@data > 0, obj_endo_all_genes@meta.data$Subcluster_Manual) 

make_paired_barplot_full <- function(L_Gene, R_Gene) {
	colours <- c("#fdbf6f", "#ff7f00", "#ffff99", subtype_cols["Monocyte"], subtype_cols["Kupffer"])
	names(colours) <- c("cvLSEC", "ppLSEC", "cvEndo", "Monocyte", "Kupffer")
	L_mac <- mac_detect[L_Gene,!colnames(mac_detect) %in% c("Debris", "Doublet", "HepContam")]
	L_endo <- endo_detect[L_Gene,!colnames(endo_detect) %in% c("Debris", "Doublet", "HepContam")]
	R_mac <- mac_detect[R_Gene,!colnames(mac_detect) %in% c("Debris", "Doublet", "HepContam")]
	R_endo <- endo_detect[R_Gene,!colnames(endo_detect) %in% c("Debris", "Doublet", "HepContam")]

	# bar order = Activated, Inflam, InflamActiv, InfSynap, 
	#			MHCII, NonInf, PhagoNonInf, ResNonInf, 
	#			cvEndo, cvLSEC, interLSEC, ppLSEC, VasEndo
	bar.colours <- c(colours["Kupffer"], colours["Kupffer"], colours["Kupffer"], colours["Kupffer"],
				colours["Kupffer"], colours["Kupffer"], colours["Kupffer"], colours["Kupffer"],
				colours["cvEndo"], colours["cvLSEC"], colours["cvEndo"], colours["ppLSEC"], colours["cvEndo"])

	par(mfrow=c(1,2))
	# Left = Ligand
	par(mar=c(4,6,2,1))
	barplot(c(L_mac, L_endo), las=1, horiz=T, main=L_Gene, xlab="Detection Rate", col=bar.colours)

	# Right = Receptor
	par(mar=c(4,1,2,6))
	ats <- barplot(-1*c(R_mac, R_endo), las=1, horiz=T, main=R_Gene, xlab="Detection Rate", , col=bar.colours, names="", xaxt="n")
	axis(1, at=seq(-0, -1*round(max(c(R_mac, R_endo)), digits=1), by=-0.1), labels=-1*seq(-0, -1*round(max(c(R_mac, R_endo)), digits=1), by=-0.1))
	axis(4, at=ats, labels=c(names(R_mac), names(R_endo)), las=2, tick=FALSE)
}


make_paired_barplot_trimmed <- function(L_Gene, R_Gene, add=FALSE, xlab="Detection Rate") {
	colours <- c("#fdbf6f", "#ff7f00", "#ffff99", subtype_cols["Monocyte"], subtype_cols["Kupffer"])
	names(colours) <- c("cvLSEC", "ppLSEC", "cvEndo", "Monocyte", "Kupffer")
	L_mac <- mac_detect[L_Gene,colnames(mac_detect) %in% names(colours)]
	L_endo <- endo_detect[L_Gene,colnames(endo_detect) %in% names(colours)]
	R_mac <- mac_detect[R_Gene,colnames(mac_detect) %in% names(colours)]
	R_endo <- endo_detect[R_Gene,colnames(endo_detect) %in% names(colours)]

	# bar order = Inflam, NonInf, cvEndo, cvLSEC
	bar.colours <- c(colours["Monocyte"], colours["Kupffer"], colours["cvEndo"], colours["cvLSEC"], colours["ppLSEC"])

	if (!add) {
		par(mfrow=c(1,2))
	}
	if (xlab == "") {
		x_mar = 2;
	} else {
		x_mar = 4
	}

	# Left = Ligand
	par(mar=c(x_mar,6,2,1))
	barplot(c(L_mac, L_endo), las=1, horiz=T, main=L_Gene, xlab=xlab, col=bar.colours)

	# Right = Receptor
	par(mar=c(x_mar,1,2,6))
	ats <- barplot(-1*c(R_mac, R_endo), las=1, horiz=T, main=R_Gene, xlab=xlab, , col=bar.colours, names="", xaxt="n")
	by_interval = -1*signif(max(c(R_mac, R_endo))/5, digits=1)
	axis(1, at=seq(from=-0, to=-1*signif(max(c(R_mac, R_endo)), digits=1), by=by_interval), labels=-1*seq(-0, -1*signif(max(c(R_mac, R_endo)), digits=1), by=by_interval))
	axis(4, at=ats, labels=c(names(R_mac), names(R_endo)), las=2, tick=FALSE)
}


png("Figure_Inflam_cvEndo1.png", width=4.5, height=3.5, units="in", res=300)
par(mfrow=c(2,2))
make_paired_barplot_trimmed("SELL", "CD34",add=TRUE, xlab="")
make_paired_barplot_trimmed("SELL", "PODXL",add=TRUE)
dev.off()

png("Figure_Inflam_cvEndo2.png", width=4.5, height=3.5, units="in", res=300)
par(mfrow=c(2,2))
#make_paired_barplot_trimmed("HLA-F", "LILRB2",add=TRUE, xlab="") # not specific to Mac in full
make_paired_barplot_trimmed("SELP", "SELPLG",add=TRUE, xlab="") # not specific to Mac in full
#make_paired_barplot_trimmed("ITGAL", "ICAM1",add=TRUE, xlab="") # not specific to InfMac in full
make_paired_barplot_trimmed("SPN", "ICAM1",add=TRUE) # not specific to InfMac in full
dev.off()



# Random Sample Heatmap

set.seed(1891)
doublets <- obj[,sample(which(obj@meta.data$Subcluster_Manual == "Doublet"), 50)]
cvLSECs <- obj_endo[,sample(which(obj_endo@meta.data$Subcluster_Manual == "cvLSEC"), 20)]
ppLSECs <- obj_endo[,sample(which(obj_endo@meta.data$Subcluster_Manual == "ppLSEC"), 20)]
Endo <- obj_endo[,sample(which(obj_endo@meta.data$Subcluster_Manual %in% c("cvEndo", "VasEndo")), 20)]
InfMacs <- obj[,sample(which(obj@meta.data$Subcluster_Manual == "Inflam"), 20)]
NonInfMacs <- obj[,sample(which(obj@meta.data$Subcluster_Manual == "NonInf"), 20)]

markers <- c(label_genes(c("MARCO", "CD5L", "CD163", "C1QC", "HLA-DRA", "HLA-DPA1"), "NonMac"), 
		label_genes(c("S100A4", "LYZ", "S100A9", "VCAN", "FCN1"), "InfMac"),
		label_genes(c("FCN2", "STAB1", "CLEC1B", "CLEC4M"), "cvLSEC"), 
		label_genes(c("IGFBP7", "S100A6", "AQP1"), "ppLSEC"),
		label_genes(c("VWF", "PODXL", "PLAC8", "PLVAP", "RSPO3"), "Endo")
		)


sets <- list(doublets, cvLSECs, ppLSECs, Endo, InfMacs, NonInfMacs)
allmat <- c();
for (obj in sets) {
	mat <- obj@assays$RNA@scale.data[match(markers, rownames(obj)),]
	allmat <- cbind(allmat, mat)
}

origin <- factor(c(rep("Doublet", 50), rep("cvLSEC",20), rep("ppLSEC",20), rep("Endo", 20), rep("InfMac", 20), rep("NonMac", 20)),
		levels=c("Doublet", "cvLSEC", "ppLSEC", "Endo", "InfMac", "NonMac"))
gene_origin <- factor(names(markers), levels=c("Doublet", "cvLSEC", "ppLSEC", "Endo", "InfMac", "NonMac"))
colours <- c("grey65", "#fdbf6f", "#ff7f00", "#ffff99", "#a6cee3", "#1f78b4")


require(gplots)
require(Seurat)
heatcol=PurpleAndYellow(20)
png("Mac_Endo_Doublets_Examples.png", width=5, height=8, units="in", res=300)
heatmap.2(allmat, trace="none", scale="none", ColSideColors=colours[origin], RowSideColors=colours[gene_origin], 
		Colv=FALSE, Rowv=FALSE, dendrogram="none", col=PurpleAndYellow(20), cexRow=1.3, cexCol=1.3, mar=c(7,7),
		key.title="", key.xlab="Scaled Expression", symkey=F, density.info="none")
dev.off()


png("Mac_Endo_Doublets_Examples_Legend.png", width=4, height=4, units="in", res=150)
plot(1,1, col="white")
legend("left", levels(origin), fill=colours, bty="n")
dev.off()

# Dotplots
Cellphonedb_genes_MAC <- c("CD2", "FPR3", "FPR2", "CD74", "C5AR1", "ADORA3", "CD44", "SELPLG", "PLXNB2", "LILRB2", "NOTCH4",
				"NOTCH1", "NOTCH2", "CD46", "NOTCH3")
Cellphonedb_genes_LSEC <- c("CD58", "HEBP1", "APP", "RPS19", "ENTPD1", "SELE", "SELP", "SELL", "SEMA4G", "HLA-F", "CD1D",
				"DLL1", "DLL4", "JAG2", "JAG1")
png("MAC2LSEC_Cellphonedb_Dotplot_Mac.png", width=6, height=4, units="in", res=300)
DotPlot(obj_all_genes, group.by="Subcluster_Manual", features=Cellphonedb_genes_MAC) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
png("MAC2LSEC_Cellphonedb_Dotplot_LSEC.png", width=6, height=4, units="in", res=300)
DotPlot(obj_endo_all_genes, group.by="Subcluster_Manual", features=Cellphonedb_genes_LSEC) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()





# -------------------- Defunct ----------------- #



## Frequency Plots ##

obj <- obj[,obj$Subcluster_Manual != "Doublet"]
obj <- obj[,obj$Subcluster_Manual != "Debris"]

png("Macrophage_Sex_Proportions.png", width=6, height=8, units="in", res=300)
barplot_freq(obj, cluster_col="Subcluster_Manual", colours=colours_sex, metadata="donor_sex")
dev.off()
png("Macrophage_Age_Proportions.png", width=6, height=8, units="in", res=300)
barplot_freq(obj, cluster_col="Subcluster_Manual", colours=colours_age, metadata="donor_age_group")
dev.off()
png("Macrophage_Rejection_Proportions.png", width=6, height=8, units="in", res=300)
barplot_freq(obj, cluster_col="Subcluster_Manual", colours=colours_reject, metadata="trans.rejected")
dev.off()


png("Macrophage_donor_Frequency.png", width=9, height=6, units="in", res=300)
indivar_freq(obj@meta.data$Subcluster_Manual, obj@meta.data$donor, obj@meta.data)
dev.off()

cluster_colours <- get_seurat_colours(obj, cluster_col)

cluster_colours <- cluster_colours[c(9, 7, 6, 1, 8, 2, 11, 5)]

indi_expr(obj, tag="Macrophage", width=12, bar_col=cluster_colours) 

## ----- For Paper ----- ##
table(obj@meta.data$Subcluster_Manual)
cluster_colours <- get_seurat_colours(obj, "Subcluster_Manual")
names(cluster_colours) <- levels(factor(obj@meta.data$Subcluster_Manual))
cluster_colours <- cluster_colours[!(names(cluster_colours) %in% c("Debris", "Doublet"))]
cleaned_obj <- obj[, !(obj@meta.data$Subcluster_Manual %in% c("Debris", "Doublet"))]
# Freequency #
png("Macrophage_paper_donor_Frequency.png", width=6, height=3, units="in", res=300)
out <- indivar_freq(cleaned_obj@meta.data$Subcluster_Manual, cleaned_obj@meta.data$donor, cleaned_obj@meta.data)
dev.off()
donor_sex <- donor2meta(cleaned_obj@meta.data, "donor_sex")
donor_age <- as.numeric(as.character(donor2meta(cleaned_obj@meta.data, "donor_age")))
donor_bmi <- donor2meta(cleaned_obj@meta.data, "donor_bmi_group")
donor_bmi <- apply()

for (i in 1:ncol(out)){
print(summary(lm(out[,i]~donor_age)))
}


# Similarity #
donor_col="donor"
ns <- table(cleaned_obj@meta.data[,donor_col])
mat <- get_pseudobulk_means(cleaned_obj@assays$RNA@data, factor(cleaned_obj@meta.data[,"Subcluster_Manual"]), cleaned_obj@meta.data[,donor_col])

mat <- mat[rownames(mat) %in% VariableFeatures(cleaned_obj),]
id <- factor(unlist(lapply(strsplit(colnames(mat), "_"), function(x){x[[1]]})))
d <- as.dist(1-cor(mat, method="pearson"))

scores <- c();
stderr <- c();
for (i in unique(id)) {
       scores <- c(scores, mean(as.matrix(1-d)[id==i, id==i]))
       stderr <- c(stderr, sd(as.matrix(1-d)[id==i, id==i])/sqrt(sum(id==i)^2)) # Fix 25 May 2021
}
png("Macrophage_similarity_paper.png", width=5.5, height=4, units="in", res=300)
par(mar=c(4,6,1,1))
bar_loc <- barplot(rev(scores), name=rev(unique(id)), xlab="Cross Donor Correlation (average)", xlim=c(0,1), las=1, col=rev(cluster_colours), horiz=T)
arrows(rev(scores),bar_loc,  rev(scores+2*stderr), bar_loc,angle=90, len=0)
#legend("bottomright", lty=1, c("95% CI"), bty="n")
dev.off()



#### DO DE - MM - NEBULA ####

for(assay_type %in% c("5pr", "3pr")) {

require(nebula)
tofit <- obj_all_genes@assays$RNA@counts[,obj@meta.data$assay_type==assay_type & !is.na(obj@meta.data$donor_bmi)];
tofit <- tofit[rowSums(tofit > 0) > 20,]
metadata <- obj@meta.data[!is.na(obj@meta.data$donor_bmi) & obj@meta.data$assay_type==assay_type,]
predictors <- model.matrix(~metadata$donor_age+metadata$donor_sex+metadata$donor_bmi)
colnames(predictors) <- c("(Intercept)", "Age", "SexM", "BMI")

test <- nebula(tofit, metadata$sample, pred=predictors, offset=colSums(tofit))
output <- test$summary; 
output$q_Age <- p.adjust(output$p_Age)
output$q_Sex <- p.adjust(output$p_Sex)
output$q_BMI <- p.adjust(output$p_BMI)

write.table(output, paste("Macrophage_Nebula", assay_type, "DE.csv", sep="_"), sep=",", row.names=FALSE, col.names=TRUE, quote=TRUE)
}

# Subclusters
for (subtype in unique(obj@meta.data$Subcluster_Manual)) {

	assay_type = "3pr"

	require(nebula)
	tofit <- obj_all_genes@assays$RNA@counts[, obj@meta.data$Subcluster_Manual == subtype & obj@meta.data$assay_type==assay_type & !is.na(obj@meta.data$donor_bmi)];
	tofit <- tofit[rowSums(tofit > 0) > 20,]
	metadata <- obj@meta.data[obj@meta.data$Subcluster_Manual == subtype & !is.na(obj@meta.data$donor_bmi) & obj@meta.data$assay_type==assay_type,]
	predictors <- model.matrix(~metadata$donor_age+metadata$donor_sex+metadata$donor_bmi)
	colnames(predictors) <- c("(Intercept)", "Age", "SexM", "BMI")
	
	test <- nebula(tofit, metadata$sample, pred=predictors, offset=colSums(tofit))
	output <- test$summary; 
	output$q_Age <- p.adjust(output$p_Age)
	output$q_Sex <- p.adjust(output$p_Sex)
	output$q_BMI <- p.adjust(output$p_BMI)

	write.table(output, paste("Macrophage_Nebula", assay_type, "Subcluster", subtype, "DE.csv", sep="_"), sep=",", row.names=FALSE, col.names=TRUE, quote=TRUE)

}




#### DO DE - MM ####
# Nebula
require(nebula)
obj <- obj[,obj$Subcluster_Manual != "Doublet"]
obj <- obj[,obj$Subcluster_Manual != "Debris"]

obj_3pr <- obj[,obj@meta.data$assay_type == "3pr"]

predictors = model.matrix(~obj_3pr@meta.data$donor_sex+obj_3pr@meta.data$donor_age)
low_expr_threshold = 0.005 # (counts per cell < 0.5%)

de <- nebula(obj_3pr@assays$RNA@counts, id = obj_3pr@meta.data$sample, pred=predictors, cpc=low_expr_threshold)

de_Female <- de$summary[p.adjust(de$summary[,8]) < 0.05 & de$summary[,2] < 0,]
de_Male <- de$summary[p.adjust(de$summary[,8]) < 0.05 & de$summary[,2] > 0,]
de_Old <- de$summary[p.adjust(de$summary[,9]) < 0.05 & de$summary[,3] > 0,]
de_Young <- de$summary[p.adjust(de$summary[,9]) < 0.05 & de$summary[,3] < 0,]

de_Female <- de$summary[de$summary[,8] < 0.05 & de$summary[,2] < 0 & de$convergence == 1,]
de_Male <- de$summary[de$summary[,8] < 0.05 & de$summary[,2] > 0 & de$convergence == 1,]
de_Old <- de$summary[de$summary[,9] < 0.05 & de$summary[,3] > 0 & de$convergence == 1,]
de_Young <- de$summary[de$summary[,9] < 0.05 & de$summary[,3] < 0 & de$convergence == 1,]


scores_sex <- de$summary[,2]; names(scores_sex) <- de$summary$gene;
scores_age <- de$summary[,3]; names(scores_age) <- de$summary$gene;

rich <- do_fgsea(scores_sex, pathways=immune_path)
rich_Age <- do_fgsea(scores_age, pathways=immune_path)


extract_pathway_dat <- function(res, this_pathway, this_name) {
	genes <- unlist(res$rich[unlist(res$rich[,1]) == this_pathway,8])
	names(genes) <- rep(this_name, length(genes));

	bar_point <- unlist(abs(log10(res$rich[unlist(res$rich[,1]) == this_pathway,"padj"]))); names(bar_point)<- this_name;
	return(list(genes=genes, point=bar_point))
}

LPS_Stim <- extract_pathway_dat(rich_Age, "GSE9988_LPS_VS_VEHICLE_TREATED_MONOCYTE_UP", "LPS Stimulation")
best_genes <- de_Young$gene[de_Young$gene %in% LPS_Stim$genes]

g = best_genes[3]
indivar_DE_vis(obj_3pr@assays$RNA@data[rownames(obj_3pr) == g,], obj_3pr@meta.data$donor, metadata=obj_3pr@meta.data, type="age")


# MAST + random effect DE #
require(Seurat)
require(MAST)
require(ggplot2)
#obj <- obj_full[,obj_full@meta.data$assay_type == "5pr"]

obj_sce <- as.SingleCellExperiment(obj)
tmp <- as.matrix(obj@assays$RNA@data); rownames(tmp) <- rownames(obj); colnames(tmp) <- colnames(obj)
obj_sca <- FromMatrix(tmp, colData(obj_sce), rowData(obj_sce), class="SingleCellAssay")

zlm.out2 <- zlm(~donor_sex + donor_age + donor_bmi + assay_type + (1|sample) + (1|donor), obj_sca, method="glmer", ebayes=FALSE) # changed May 31 to include assay type, and donor too. Helps a bit but not a lot.


zlm.out2 <- zlm(~donor_sex + (1|sample:donor_sex), obj_sca, method="glmer", ebayes=FALSE) 
coefAndCI <- as.data.frame(summary(zlm.out2, logFC=FALSE)$datatable)

res_sex <- coefAndCI[coefAndCI[,3] == "donor_sexM" & !is.na(coefAndCI[,6]),]
res_age <- coefAndCI[coefAndCI[,3] == "donor_age" & !is.na(coefAndCI[,6]),]

max_coef.age <- aggregate(res_age$coef, by=list(res_age$primerid), function(x) {x[abs(x) == max(abs(x))]})
max_coef.sex <- aggregate(res_sex$coef, by=list(res_sex$primerid), function(x) {x[abs(x) == max(abs(x))]})


zlm.lr.sex <- lrTest(zlm.out2, "donor_sex")
zlm.lr.age <- lrTest(zlm.out2, "donor_age")

pval.sex <- zlm.lr.sex[,,3][,"hurdle"]
pval.age <- zlm.lr.age[,,3][,"hurdle"]

sex.output <- cbind(max_coef.sex, pval.sex[match( max_coef.sex[,1], names(pval.sex))])
sex.output$qval <- p.adjust(sex.output[,3], method="fdr")
age.output <- cbind(max_coef.age, pval.age[match( max_coef.age[,1], names(pval.age))])
age.output$qval <- p.adjust(age.output[,3], method="fdr")

print(sum(age.output$qval < 0.05))
print(sum(sex.output$qval < 0.05))

# >7000 age, > 4000 sex out of 10,000 total genes - same if use normalized, 3' only, 5' only, or scaled data.

saveRDS(list(sex=sex.output, age=age.output), "Macrophage_MM_DE_5pr.rds")

## Sex ##
set.seed(910)
expr_mat <- obj@assays$RNA@scale.data
label <- obj@meta.data$donor_sex


de_sex <- t(apply(expr_mat, 1, function(x) {
		res = t.test(x[label=="F"], x[label=="M"])
		return(c(res$estimate, res$p.value))
		}))
de_sex <- cbind(de_sex, p.adjust(de_sex[,3], "fdr"))
de_sex <- cbind(de_sex, -1*log(de_sex[,3])*sign(de_sex[,1]-de_sex[,2]))
colnames(de_sex) <- c("F", "M", "pval", "qval", "score")

de_clean <- de_sex[!grepl("^RPS", rownames(de_sex)),]
de_clean <- de_clean[!grepl("^RPL", rownames(de_clean)),]

richments_de <- do_fgsea(de_clean[,5], reactome, fdr=0.05)

require(igraph)
names(V(richments_de$graph))

fix_names <- c("ROS Detox", "Oxidation", "Neutrophil degranulation", "Fatty acid", "Xenobiotics", 
	"Purine catabolism", "", "Protein local", "AA Met", "Carbohydrate Met", "Innate Immune", 
	"Stress response", "", "Lipid Met", "RNA Met", "SREBF", "mt-tRNA", "mt-rRNA", "rRNA", "tRNA")
G <- set.vertex.attribute(richments_de$graph, "name", value=fix_names)

png("Macrophage_DE_Reactome_Sex_pathways.png", width=8, height=8, units="in", res=300)
set.seed(2910)
plot(G, vertex.color=richments_de$vertex_col, vertex.size=richments_de$vertex_size)
dev.off()

## Peudobulk ##

pseudo_sex <- get_pseudobulk_means(expr_mat, obj@meta.data$donor_sex, obj@meta.data$sample)
is.F <- grepl("^F", colnames(pseudo_sex))
de_pseudo <- t(apply(pseudo_sex, 1, function(x) {
		res = t.test(x[is.F], x[!is.F])
		return(c(res$estimate, res$p.value))
		}))
de_pseudo <- cbind(de_pseudo, p.adjust(de_pseudo[,3], "fdr"))
de_pseudo <- cbind(de_pseudo, -1*log(de_pseudo[,3])*sign(de_pseudo[,1]-de_pseudo[,2]))
colnames(de_pseudo) <- c("F", "M", "pval", "qval", "score")


saveRDS(list(de=de_sex, richments=richments_de, fix_names = fix_names, 
		pseudo_de=de_pseudo, pseudo_table=pseudo_sex), 
		file="Macrophage_sex_de.rds")




## Rejection ##
set.seed(1672)
expr_mat <- obj@assays$RNA@scale.data
label <- obj@meta.data$trans.rejected


de_rej <- t(apply(expr_mat, 1, function(x) {
		res = t.test(x[label=="Y"], x[label=="N"]) 
		return(c(res$estimate, res$p.value))
		}))
de_rej <- cbind(de_rej, p.adjust(de_rej[,3], "fdr"))
de_rej <- cbind(de_rej, -1*log(de_rej[,3])*sign(de_rej[,1]-de_rej[,2]))
colnames(de_rej) <- c("Y", "N", "pval", "qval", "score")

de_clean <- de_rej[!grepl("^RPS", rownames(de_rej)),]
de_clean <- de_clean[!grepl("^RPL", rownames(de_clean)),]

richments_de <- do_fgsea(de_clean[,5], reactome, fdr=0.05)

require(igraph)
names(V(richments_de$graph))

fix_names <- c("Protein phosphorylation", "IGF", "Xenobiotics", "Oxidation", "Cytochrome P450", 
	"Complement cascade", "PPARA", "", "Vitamin Met", "Retinoid", "Clotting cascade", "", "",
	"Lipid Met", "Fatty acid", "", "mt-rRNA", "mt-tRNA", "Scavenger Receptors", "BRAF", 
	"Cell Cycle", "BCR Signalling", "Antigen presentation", "Interferon a/b",
	"RCERI", "Immunological synapse", "ER-Phagosome", "Nef", "TCR Phosphorylation",
	"HIV", "DAP12", "", "PD1", "HIV", "2nd messenger", "CD28", "Adaptive Immune",
	"TCR signaling", "Lymphoid/Non-Lymphoid", "TCR"
	)
G <- set.vertex.attribute(richments_de$graph, "name", value=fix_names)

png("Macrophage_DE_Reactome_Rejection_pathways.png", width=8, height=8, units="in", res=300)
set.seed(2910)
plot(G, vertex.color=richments_de$vertex_col, vertex.size=richments_de$vertex_size)
dev.off()


## Peudobulk ##

pseudo_reject <- get_pseudobulk_means(expr_mat, obj@meta.data$trans.rejected, obj@meta.data$sample)
reject <- grepl("^Y", colnames(pseudo_reject))
notreject <- grepl("^N", colnames(pseudo_reject))
de_pseudo <- t(apply(pseudo_reject, 1, function(x) {
		res = t.test(x[reject], x[notreject])
		return(c(res$estimate, res$p.value))
		}))
de_pseudo <- cbind(de_pseudo, p.adjust(de_pseudo[,3], "fdr"))
de_pseudo <- cbind(de_pseudo, -1*log(de_pseudo[,3])*sign(de_pseudo[,1]-de_pseudo[,2]))
colnames(de_pseudo) <- c("Reject", "NotReject", "pval", "qval", "score")


saveRDS(list(de=de_rej, richments=richments_de, fix_names = fix_names, 
		pseudo_de=de_pseudo, pseudo_table=pseudo_reject), 
		file="Macrophage_rejection_de.rds")




## Age ##
set.seed(8391)
expr_mat <- obj@assays$RNA@scale.data
label <- obj@meta.data$donor_age_group


de_age <- t(apply(expr_mat, 1, function(x) {
		res = t.test(x[label=="elderly"], x[label=="young"])
		return(c(res$estimate, res$p.value))
		}))
de_age <- cbind(de_age, p.adjust(de_age[,3], "fdr"))
de_age <- cbind(de_age, -1*log(de_age[,3])*sign(de_age[,1]-de_age[,2]))
colnames(de_age) <- c("elderly", "young", "pval", "qval", "score")

de_clean <- de_age[!grepl("^RPS", rownames(de_age)),]
de_clean <- de_clean[!grepl("^RPL", rownames(de_clean)),]

richments_de <- do_fgsea(de_clean[,5], reactome, fdr=0.05)

require(igraph)
names(V(richments_de$graph))

fix_names <- c(
	"Complement cascade", "", "Heme scavenging", "Steroid met", "Nucleobase catablism", "lipoprotein",
	"Cholesterol", "", "Protein phophorylation", "Gluconeogenesis", "Scavenger Receptor", "SREBF",
	"IGF", "", "Fatty acids", "Bile acid", "", "", "", "Cell-Vascular interaction", "TGF-beta down", 
	"Senescence", "mt-rRNA", "TGF-beta SMADs", "Cell death", "Autophagy", "NOTCH2", "", 
	"Neurodegeneration", "", "Oxidative Stress", "HIV", "Mitophagy", "Heat shock", "Steroid Chaperone",
	"Aggrephagy", "Endosomal Transport", "HSF1", "Interferon a/b", "Selective autophagy")
G <- set.vertex.attribute(richments_de$graph, "name", value=fix_names)

png("Macrophage_DE_Reactome_Age_pathways.png", width=8, height=8, units="in", res=300)
set.seed(2910)
plot(G, vertex.color=richments_de$vertex_col, vertex.size=richments_de$vertex_size)
dev.off()


## Peudobulk ##

pseudo_age <- get_pseudobulk_means(expr_mat, obj@meta.data$donor_age_group, obj@meta.data$sample)
young <- grepl("^young", colnames(pseudo_age))
elderly <- grepl("^elderly", colnames(pseudo_age))
de_pseudo <- t(apply(pseudo_age, 1, function(x) {
		res = t.test(x[elderly], x[young])
		return(c(res$estimate, res$p.value))
		}))
de_pseudo <- cbind(de_pseudo, p.adjust(de_pseudo[,3], "fdr"))
de_pseudo <- cbind(de_pseudo, -1*log(de_pseudo[,3])*sign(de_pseudo[,1]-de_pseudo[,2]))
colnames(de_pseudo) <- c("Elderly", "Young", "pval", "qval", "score")

saveRDS(list(de=de_age, richments=richments_de, fix_names = fix_names, 
		pseudo_de=de_pseudo, pseudo_table=pseudo_age), 
		file="Macrophage_age_de.rds")


### --------- Pathways Subtypes  -------------- ###
obj <- obj_full
cluster_col = "Subcluster_Manual"
obj <- obj[, (obj@meta.data[,cluster_col] !=  "Debris")]

table(obj@meta.data[,cluster_col])

ident1 = "InfSynap" # "MHCII" # "ResNonInf" "Activated", "InfSynap", "PhagoNonInf"
de <- run_wilcox(obj, obj@meta.data[,cluster_col], ident.1=ident1)
de <- de[order(de$log2fc, decreasing=T),]

write.table(de, paste("Macrophage", ident1, "wilcox_de_for_pathways.csv", sep="_"), sep=",")

source("C:/Users/tandrews/Documents/UHNSonya/FGSEA_enrichment_script.R")

tmp <- de[,"log2fc"]
names(tmp) <- rownames(de)
tmp <- tmp[!grepl("^RP[LS]", names(tmp))]
tmp[tmp == Inf] <- max(tmp[is.finite(tmp)])+1
tmp[tmp == -Inf] <- min(tmp[is.finite(tmp)])-1
tmp <- tmp[!is.na(tmp)]
res_immune2 <- do_fgsea(tmp, pathway=immune_path)
res_kegg <- do_fgsea(tmp, pathway=BaderKegg)
res_reactome <- do_fgsea(tmp, pathway=reactome)
res_msig <- do_fgsea(tmp, pathway=MSigAll)
res_hallmark <- do_fgsea(tmp, pathway=Hallmark_path)
res_bader <- do_fgsea(tmp, pathway=BaderMSig)


all_rich <- list(immune = res_immune2, kegg = res_kegg, msig=res_msig, react=res_reactome , hallmark=res_hallmark, bader=res_bader)
saveRDS(all_rich, paste("Macrophage", ident1, "patheways.rds", sep="_"))
#this_res <- res_sig
#plot(this_res$graph, vertex.color=this_res$vertex_col, vertex.size=this_res$vertex_size)
#tail(unlist(this_res$rich[,1]), 30)

detection_rate <- group_rowmeans(obj_full@assays$RNA@counts, obj_full@meta.data$Subcluster_Manual)


extract_pathway_dat <- function(res, this_pathway, this_name) {
	genes <- unlist(res$rich[unlist(res$rich[,1]) == this_pathway,8])
	names(genes) <- rep(this_name, length(genes));

	bar_point <- unlist(abs(log10(res$rich[unlist(res$rich[,1]) == this_pathway,"padj"]))); names(bar_point)<- this_name;
	return(list(genes=genes, point=bar_point))
}


### MHCII ###

all_rich <- readRDS("Macrophage_MHCII_patheways.rds")
res_immune2<- all_rich$immune
res_kegg <- all_rich$kegg
res_reactome <- all_rich$reactome
res_hallmark <- all_rich$hallmark
res_msig <- all_rich$msig

immune_cd16Mono <- extract_pathway_dat(res_immune2, "GSE34515_CD16_NEG_VS_POS_MONOCYTE_DN", "CD16- Monocyte")
immune_DC <- extract_pathway_dat(res_immune2, "GSE22886_DC_VS_MONOCYTE_UP", "DC")
kegg_antigen <- extract_pathway_dat(res_kegg, "Antigen processing and presentation%KEGG%hsa04612", "MHCII")
hallmark_myc <- extract_pathway_dat(res_hallmark, "HALLMARK_MYC_TARGETS_V1", "Myc targets")
hallmark_ifna <- extract_pathway_dat(res_hallmark, "HALLMARK_INTERFERON_ALPHA_RESPONSE", "IFNa response")
reactome_ifg <-  extract_pathway_dat(res_reactome, "Interferon gamma signaling", "IFNg signaling")


all_genes <- c(immune_cd16Mono$genes, kegg_antigen$genes, hallmark_myc$genes, hallmark_ifna$genes, reactome_ifg$genes)
all_genes <- all_genes[all_genes %in% rownames(detection_rate)[detection_rate[,iden1] > 0.1]]
all_genes <- all_genes[!duplicated(all_genes)]
DotPlot(obj_full, features=all_genes, group.by="Subcluster_Manual")

cluster_colours <- get_seurat_colours(obj_full, "Subcluster_Manual")
png("Macrophage_MHCII_pathways.png", width=4*0.8, height=3*0.8, units="in", res=150)
par(mar=c(4,7.5,1,1))
barplot(c(hallmark_ifna$point, reactome_ifg$point, immune_cd16Mono$point, kegg_antigen$point, hallmark_myc$point), col=cluster_colours[7], horiz=T, las=1, xlab="log10(p-value)")
dev.off()



### Retinol/Retinoid ###
all_rich <- readRDS("Macrophage_ResNonInf_patheways.rds")
res_immune2<- all_rich$immune
res_kegg <- all_rich$kegg
res_reactome <- all_rich$reactome
res_hallmark <- all_rich$hallmark
res_msig <- all_rich$msig

immune_earlyMono <- extract_pathway_dat(res_immune2, "GSE22886_DAY0_VS_DAY7_MONOCYTE_IN_CULTURE_DN", "Monocyte\nNaive") # highest in phago
immune_ppargKO <- extract_pathway_dat(res_immune2, "GSE25123_WT_VS_PPARG_KO_MACROPHAGE_ROSIGLITAZONE_STIM_UP", "PPARG KO")
react_lipoprotein <- extract_pathway_dat(res_reactome, "Plasma lipoprotein remodeling", "Lipoprotein\nRemodeling") # highest in debris
react_iron <- extract_pathway_dat(res_reactome, "Iron uptake and transport", "Iron Uptake")
react_water <- extract_pathway_dat(res_reactome, "Vasopressin regulates renal water homeostasis via Aquaporins", "Water")
msig_oxphos <- extract_pathway_dat(res_msig, "MOOTHA_VOXPHOS", "Oxidative\nPhosphorylation")
msig_comp <- extract_pathway_dat(res_msig, "BIOCARTA_COMP_PATHWAY", "Complement\npathway")
hallmark_secretion <- extract_pathway_dat(res_hallmark, "HALLMARK_PROTEIN_SECRETION", "Protein\nSecretion")
hallmark_androgen <- extract_pathway_dat(res_hallmark, "HALLMARK_ANDROGEN_RESPONSE", "androgen") 


all_genes <- c(immune_DC$genes, hallmark_secretion$genes, msig_comp$genes, msig_oxphos$genes, immune_ppargKO$genes, react_iron$genes, react_water$genes)
all_genes <- all_genes[all_genes %in% rownames(detection_rate)[detection_rate[,"Activated"] > 0.1]]
all_genes <- all_genes[!duplicated(all_genes)]
DotPlot(obj_full, features=all_genes, group.by="Subcluster_Manual")

cluster_colours <- get_seurat_colours(obj_full, "Subcluster_Manual")
png("Macrophage_Resident_pathways_big.png", width=4*2, height=3*2, units="in", res=150)
par(mar=c(4,7.5,1,1))
barplot(rev(c(immune_ppargKO$point, immune_earlyMono$point, msig_oxphos$point, msig_comp$point, hallmark_secretion$point, react_iron$point)), col=cluster_colours[10], horiz=T, las=1, xlab="log10(p-value)")
dev.off()

immune_DC <- extract_pathway_dat(res_immune2, "GSE22886_DC_VS_MONOCYTE_UP", "Dendritic")
immune_earlyMono <- extract_pathway_dat(res_immune2, "GSE22886_DAY0_VS_DAY7_MONOCYTE_IN_CULTURE_DN", "Monocyte") # highest in phago
msig_oxphos <- extract_pathway_dat(res_msig, "MOOTHA_VOXPHOS", "Ox-Phosph")
msig_comp <- extract_pathway_dat(res_msig, "BIOCARTA_COMP_PATHWAY", "Complement")
hallmark_secretion <- extract_pathway_dat(res_hallmark, "HALLMARK_PROTEIN_SECRETION", "Secretion")
react_iron <- extract_pathway_dat(res_reactome, "Iron uptake and transport", "Iron Uptake")
png("Macrophage_Resident_pathways_small.png", width=4*0.8, height=3*0.8, units="in", res=150)
par(mar=c(4,7.5,1,1))
barplot(rev(c( immune_earlyMono$point, immune_DC$point, msig_oxphos$point, msig_comp$point, hallmark_secretion$point, react_iron$point)), col=cluster_colours[10], horiz=T, las=1, xlab="log10(p-value)")
dev.off()




### Clearance/Phagocytosis ###
all_rich <- readRDS("Macrophage_PhagoNonInf_patheways.rds")
res_immune2<- all_rich$immune
res_kegg <- all_rich$kegg
res_reactome <- all_rich$reactome
res_hallmark <- all_rich$hallmark
res_msig <- all_rich$msig

immune_earlyMono <- extract_pathway_dat(res_immune2, "GSE22886_DAY0_VS_DAY7_MONOCYTE_IN_CULTURE_DN", "Monocyte\nNaive")
immune_IL2 <- extract_pathway_dat(res_immune2, "GSE36888_UNTREATED_VS_IL2_TREATED_TCELL_17H_UP", "IL2 Stim")
kegg_lyso <- extract_pathway_dat(res_kegg, "Lysosome%KEGG%hsa04142", "Lysosome")
kegg_phago <- extract_pathway_dat(res_kegg, "Phagosome%KEGG%hsa04145", "Phagosome")
kegg_glycan <- extract_pathway_dat(res_kegg, "Other glycan degradation%KEGG%hsa00511", "Glygan Degradation")
react_pd1 <- extract_pathway_dat(res_reactome, "PD-1 signaling", "PD1")
react_antigen <- extract_pathway_dat(res_reactome, "MHC class II antigen presentation", "Antigen Present")
react_collagen <- extract_pathway_dat(res_reactome, "Collagen degradation", "Collagen\nDegradation")
msig_dendritic <- extract_pathway_dat(res_msig, "LENAOUR_DENDRITIC_CELL_MATURATION_UP", "DC Maturation")


all_genes <- c(immune_IL2$genes, kegg_lyso$genes, kegg_phago$genes, msig_dendritic$genes, react_antigen$genes)
all_genes <- all_genes[all_genes %in% rownames(detection_rate)[detection_rate[,"Activated"] > 0.1]]
all_genes <- all_genes[!duplicated(all_genes)]
DotPlot(obj_full, features=all_genes, group.by="Subcluster_Manual")

cluster_colours <- get_seurat_colours(obj_full, "Subcluster_Manual")
png("Macrophage_Phagocytosis_pathways.png", width=4*0.8, height=3*0.8, units="in", res=150)
par(mar=c(4,7.5,1,1))
barplot(c(kegg_phago$point, react_antigen$point, msig_dendritic$point, kegg_lyso$point, immune_IL2$point), col=cluster_colours[9], horiz=T, las=1, xlab="log10(p-value)")
dev.off()




### InfSynap ###
all_rich <- readRDS("Macrophage_InfSynap_patheways.rds")
res_immune2<- all_rich$immune
res_kegg <- all_rich$kegg
res_reactome <- all_rich$reactome
res_hallmark <- all_rich$hallmark
res_msig <- all_rich$msig

reactome_ISyn <- extract_pathway_dat(res_reactome, "Translocation of ZAP-70 to Immunological synapse", "Immunological\nSynapse")
reactome_CD3Phos <- extract_pathway_dat(res_reactome, "Phosphorylation of CD3 and TCR zeta chains", "TCR\nPhosphorylation")
msig_IFNg <- extract_pathway_dat(res_msig, "DER_IFN_GAMMA_RESPONSE_UP", "IFNg\nResponse")
msig_foxp3 <- extract_pathway_dat(res_msig, "MARSON_FOXP3_CORE_DIRECT_TARGETS", "FOXP3 Targets")
immune_monocyte <- extract_pathway_dat(res_immune2, "GSE22886_DAY0_VS_DAY7_MONOCYTE_IN_CULTURE_UP", "Monocyte\nMaturation")
immune_naiveMono <- extract_pathway_dat(res_immune2, "GSE9988_LPS_VS_CTRL_TREATED_MONOCYTE_DN", "Mono naive")
hallmark_INFA <- extract_pathway_dat(res_hallmark, "HALLMARK_INTERFERON_ALPHA_RESPONSE", "IFNa")
hallmark_INFg <- extract_pathway_dat(res_hallmark, "HALLMARK_INTERFERON_GAMMA_RESPONSE", "IFNg\nResponse")
bader_myc <- extract_pathway_dat(res_bader, "HALLMARK_MYC_TARGETS_V1%MSIGDB_C2%HALLMARK_MYC_TARGETS_V1", "MYC Targets")
bader_tcr <- extract_pathway_dat(res_bader, "PID_TCR_PATHWAY%MSIGDB_C2%PID_TCR_PATHWAY", "TCR Pathway")
bader_IL8 <- extract_pathway_dat(res_bader, "PID_IL8_CXCR2_PATHWAY%MSIGDB_C2%PID_IL8_CXCR2_PATHWAY", "IL8-CXCR2")
bader_IL12 <- extract_pathway_dat(res_bader, "PID_IL12_2PATHWAY%MSIGDB_C2%PID_IL12_2PATHWAY", "IL12")


all_genes <- c( bader_tcr$genes, bader_IL8$genes, reactome_ISyn$genes, reactome_CD3Phos$genes, hallmark_INFA$genes,  bader_myc$genes, immune_naiveMono$genes)
all_genes <- all_genes[all_genes %in% rownames(detection_rate)[detection_rate[,"Activated"] > 0.1]]
all_genes <- all_genes[!duplicated(all_genes)]
DotPlot(obj_full, features=all_genes, group.by="Subcluster_Manual")

png("Macrophage_InfSynap_pathways.png", width=4*0.8, height=3*0.8, units="in", res=150)
par(mar=c(4,7.5,1,1))
barplot(c( bader_IL8$point,  reactome_ISyn$point, bader_tcr$point, hallmark_INFA$point, bader_myc$point), col=cluster_colours[6], horiz=T, las=1, xlab="log10(p-value)")
dev.off()



### Activated ###

all_rich <- readRDS("Macrophage_Activated_patheways.rds")
res_immune2<- all_rich$immune
res_kegg <- all_rich$kegg
res_reactome <- all_rich$reactome
res_hallmark <- all_rich$hallmark
res_msig <- all_rich$msig

immune_lps <- extract_pathway_dat(res_immune2, "GSE9988_LPS_VS_CTRL_TREATED_MONOCYTE_UP", "LPS Stimulation")
msig_tnf <- extract_pathway_dat(res_msig, "PHONG_TNF_TARGETS_UP", "TNF Targets")
hallmark_IL2 <- extract_pathway_dat(res_hallmark, "HALLMARK_IL2_STAT5_SIGNALING", "IL2-STAT5")
reactome_IL10 <- extract_pathway_dat(res_reactome, "Interleukin-10 signaling", "IL10 Signaling")
reactome_Chemo <- extract_pathway_dat(res_reactome, "Chemokine receptors bind chemokines", "Chemokine Signaling")
hallmark_IL12 <- extract_pathway_dat(res_bader, "PID_IL12_2PATHWAY%MSIGDB_C2%PID_IL12_2PATHWAY", "IL12 pathway")
hallmark_inf <- extract_pathway_dat(res_hallmark, "HALLMARK_INFLAMMATORY_RESPONSE", "Inflammation")



all_genes <- c( hallmark_IL12$genes,immune_lps$genes, msig_tnf$genes, hallmark_inf$genes, hallmark_IL2$genes)
all_genes <- all_genes[all_genes %in% rownames(detection_rate)[detection_rate[,"Activated"] > 0.1]]
all_genes <- all_genes[!duplicated(all_genes)]
DotPlot(obj_full, features=all_genes, group.by="Subcluster_Manual")

cluster_colours <- get_seurat_colours(obj, "Coarse_clusters")
png("Macrophage_Activated_pathways.png",  width=4*0.8, height=3*0.8, units="in", res=150)
par(mar=c(4,7.5,1,1))
barplot(c(hallmark_IL12$point,hallmark_IL2$point, hallmark_inf$point, msig_tnf$point, immune_lps$point), col=cluster_colours[1], horiz=T, las=1, xlab="log10(p-value)")
dev.off()

# ------------------- DE Vis ------------- #
source("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/IndividualVariability.R")
de_sex <- readRDS("Macrophage_sex_de.rds")
de_age <- readRDS("Macrophage_age_de.rds")
de_rej <- readRDS("Macrophage_rejection_de.rds")

obj <- readRDS("Macrophage_varimax_Subcluster.rds")
obj <- obj[,obj@meta.data$assay_type == "3pr"]
obj <- obj_all_genes[,obj_all_genes@meta.data$assay_type == "3pr"]

### Manual check individual genes ###
sonya_gene_list = c("CCL2", "MCP1", "CX3CL1", "CXCL12", "CCL5", "MYD88")
par(mfrow=c(2,3))
for(g in sonya_gene_list) {
	if (! (g %in% rownames(obj))) {next}
	sex_boxplot(obj@assays$RNA@scale.data[g,], obj@meta.data$donor, obj@meta.data, g)
}


type ="sex"

#top DE
if (type =="age") {
	de <- de_age
	up = "elderly"
	dn = "young"
}
if (type =="sex") {
	de <- de_sex
	up = "F"
	dn = "M"
}
if (type =="reject") {
	de <- de_rej
	up = "Y"
	dn = "N"
}

#de <- de_rej$de
de2 <- de$de
de2 <- de2[de2[,5] >0 & de2[,1] >0 | de2[,5] < 0 & de2[,2] >0,]

top_genes <- rownames(de2[order(de2[,5],decreasing=T),])[1:8] # F = N, F = M, F = young

png(paste("Macrophage", up, "Hi_genes.png", sep="_"), width=8, height=6, units="in", res=150)

par(mfrow=c(2,4))
for(g in top_genes) {
	indivar_DE_vis(obj@assays$RNA@scale.data[g,], obj@meta.data$donor, obj@meta.data, type=type)
#	indivar_DE_vis(obj@assays$RNA@data[g,], obj@meta.data$donor, obj@meta.data)

	title(main=g)
}
dev.off()


top_genes <- rownames(de2[order(de2[,5],decreasing=F),])[1:8]

png(paste("Macrophage", dn, "Hi_genes.png", sep="_"), width=8, height=6, units="in", res=150)

par(mfrow=c(2,4))
for(g in top_genes) {
	indivar_DE_vis(obj@assays$RNA@scale.data[g,], obj@meta.data$donor, obj@meta.data, type=type)
#	indivar_DE_vis(obj@assays$RNA@data[g,], obj@meta.data$donor, obj@meta.data)

	title(main=g)
}
dev.off()

#top pathways

#de = de_age

up_pathways <- de$richments$rich[de$richments$rich$NES > 0,] # F, Y, elderly
up_pathways <- up_pathways[rev(1:nrow(up_pathways)),]
dn_pathways <- de$richments$rich[de$richments$rich$NES < 0,] # M, N, young

png(paste("Macrophage", up, "Hi_pathways.png", sep="_"), width=8, height=6, units="in", res=150)

par(mfrow=c(2,4))
for(i in 1:min(8, nrow(up_pathways))) {
	genes <- unlist(up_pathways[i, "leadingEdge"])
	indivar_DE_vis(obj@assays$RNA@scale.data[genes,], obj@meta.data$donor, obj@meta.data, type=type)
	title(main=up_pathways$pathway[i])
}

dev.off()

png(paste("Macrophage", dn, "Hi_pathways.png", sep="_"), width=8, height=6, units="in", res=150)

par(mfrow=c(2,4))
for(i in 1:min(8, nrow(dn_pathways))) {
	genes <- unlist(dn_pathways[i, "leadingEdge"])
	indivar_DE_vis(obj@assays$RNA@scale.data[genes,], obj@meta.data$donor, obj@meta.data, type=type)
	title(main=dn_pathways$pathway[i])
}

dev.off()

## Random Forest ##
library(randomForest)
library(datasets)


Macrophage_genes_dot <- c(
		label_genes(c("MARCO", "CD5L", "SLC40A1", "FTL", "CD163", "SEPP1", "C1QC", "C1QB", "C1QA", "CTSB", "HMOX1", "VCAM1"), "NonInfammatory"), 
		label_genes(c("HLA-DRA", "HLA-DPA1", "HLA-DQB1", "HLA-DPB1", "CD74"), "MHCII"), "VSIG4", "NINJ1", "IL18", 
		label_genes(c("LYZ", "S100A4", "S100A6", "S100A8", "S100A9", "S100A12", "MNDA", "VCAN", "FCN1"), "Infammatory"),
		label_genes(c("ACP5", "PLD3", "PSAP", "CSTB", "LGMN"), "Lysosome"),"FTH1", "CD68", "APOE", "GLUL", 
		label_genes(c("FABP5", "GPNMB", "CD9", "SPP1", "APOC1"), "LAM"),
		"RBP7", "PLTP", label_genes(c("FOLR2", "TIMD4", "LYVE1", "FCER1G", "MS4A7", "TIMP1"),"Resident"), "CD14",
		label_genes(c("CXCL3", "THBS1", "NAMPT", "CXCL2", "CD83", "IL1B", "AREG", "CCL3","PLAUR", "SRGN"), "Activation"),
		"PLAC8", "CD54", label_genes(c("LST1", "IFITM3", "AIF1", "COTL1"), "Synapse"),
		label_genes(c("DNASE1L3", "FCN2", "CCL14", "FCN3", "SPARC", "CLEC1B", "ENG"), "LSECs")
		) #label_genes(c("ALB", "SERPINA1", "HP", "FGA"),"Hepto")

obj <- readRDS("Macrophage_varimax_Subcluster.rds")
cluster_col = "Coarse_clusters"

data <- obj@assays$RNA@data[rownames(obj) %in% Macrophage_genes_dot,]
data <- t(as.matrix(data))
require(Matrix)
data <- as.data.frame(data)

rf <- randomForest(x=data, y=factor(cell_types))

