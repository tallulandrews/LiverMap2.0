source("../../scripts/LiverMap2.0/My_R_Scripts.R")
source("SubColour_Scheme.R")
source("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/IndividualVariability.R")
source("C:/Users/tandrews/Documents/UHNSonya/scripts/LiverMap2.0/Colour_Scheme.R")
require(fgsea)
require(Seurat)
require(ggplot2)

immune_path <- gmtPathways("../../ExternalData/MSigDb_immune_signatures_c7.all.v7.1.symbols.gmt")
Hallmark_path <- gmtPathways("../../ExternalData/MSigDbHalmarkPathways.gmt")
MSigAll <- gmtPathways("../../ExternalData/MSigDb_curated_c2.all.v7.1.symbols.gmt")
reactome <- gmtPathways("../../ExternalData/ReactomePathways.gmt")


BaderMSig <- gmtPathways("../../ExternalData/BaderLab25Aug2020/Human_MSigdb_August_01_2020_symbol.gmt.txt")
BaderReact <- gmtPathways("../../ExternalData/BaderLab25Aug2020/Human_Reactome_August_01_2020_symbol.gmt.txt")



subsubcluster <- function(obj, cluster=0, cluster_col="Core_clusters") {
	set.seed(101)
	subset <- obj[,obj@meta.data[,cluster_col]==cluster]
	subset <- RunPCA(subset)
	subset <- FindVariableFeatures(subset, method="vst", nfeatures=1000)
	subset <- RunPCA(subset, features=VariableFeatures(subset), ndims=1:20)
	subset <- FindNeighbors(subset, dims=1:12)
	subset <- FindClusters(subset, resolution=0.5)
	subset <- RunUMAP(subset, dims=1:20)
	png(paste("SubSubcluster_NKT_", cluster, "_umap.png", sep=""), width=8, height=8, units="in", res=300)
	print(DimPlot(subset))
	dev.off();
	return(subset)
}


## Colour Scheme ##
col_CD8_CD3_Tcells=Cell_type_colours[Cell_type_colours[,1] == "gdTcells2",2];
col_IL7R_CD3_Tcells=Cell_type_colours[Cell_type_colours[,1] == "Tcell",2];
col_Myeloid=Cell_type_colours[Cell_type_colours[,1] == "NonInfMac",2];
col_Bcells=Cell_type_colours[Cell_type_colours[,1] == "Bcell",2];
col_Prolif=Cell_type_colours[Cell_type_colours[,1] == "gdTcells1",2];
col_Plasmablasts=Cell_type_colours[Cell_type_colours[,1] == "Eryth",2];
col_lrNK="#DB8E00";
col_cNK=Cell_type_colours[Cell_type_colours[,1] == "NKcells",2];
col_Hepatocyte1="#F8766D"
col_Hepatocyte2=Cell_type_colours[Cell_type_colours[,1] == "Hepatocyte",2]
col_Flush="black"



## Read in data ##
obj <- readRDS("NKT_varimax_Subcluster.rds")

obj_all_genes <- readRDS("AllGenes/NKT_harmony_Subcluster_Allgenes.rds")

obj_all_genes <- obj_all_genes[,match(obj@meta.data$cell_ID, obj_all_genes@meta.data$cell_ID)]


cluster_col = "Core_clusters"

# change CD3+ to CD4+ T cells 15/12/2021, changed Hepato & HepatoDoublets to Debris

manual_cluster_anno <- c("CD8+Tcell", "lrNKcell", "CD4+Tcell", "cNKcell", "Prolif", "Debris", 
		"MatBcell", "Debris", "Myeloid", "Flush", "Erythroblasts")
type_colour_scheme <- c(col_CD8_CD3_Tcells, col_lrNK, col_IL7R_CD3_Tcells, col_cNK, col_Prolif, col_Hepatocyte1,
		col_Bcells, col_Hepatocyte1, col_Myeloid, col_Flush, col_Hepatocyte2) # col_Plasmablasts
names(type_colour_scheme) <- manual_cluster_anno

obj@meta.data$Subcluster_Manual <- manual_cluster_anno[obj@meta.data[,cluster_col]]


cluster_order_levels <- c("CD4+Tcell", "CD8+Tcell", "lrNKcell", "cNKcell", "MatBcell", "Myeloid", "Prolif", "Erythroblasts", "Flush", "Debris")

#cluster_colours <- get_seurat_colours(obj, cluster_col)[c(3, 1, 2, 4, 7, 9, 5, 11, 10, 6, 8)]
#cluster_colours <- type_colour_scheme[c(3, 1, 2, 4, 7, 9, 5, 11, 10, 6, 8)]
cluster_colours <- type_colour_scheme[c(3, 1, 2, 4, 7, 9, 5, 11, 10, 6)]


obj@meta.data$Subcluster_Manual <- factor(obj@meta.data$Subcluster_Manual, 
	levels=cluster_order_levels)

names(cluster_colours) <- levels(obj@meta.data$Subcluster_Manual)

obj_all_genes@meta.data[,cluster_col] <- obj@meta.data[,cluster_col]
obj_all_genes@meta.data$Subcluster_Manual <- obj@meta.data[,"Subcluster_Manual"]
obj_full <- obj
#saveRDS(shiny_obj, "ShinyApps/NKT_LatticeObj.rds")
#saveRDS(obj@meta.data, "NKT_fullmetadata.rds")


suppl_tab1 <- get_cluster_summary(obj, samples="sample")
suppl_tab2 <- get_cluster_summary(obj, samples="donor")
write.table(suppl_tab2, file="NKT_cluster_summary.csv", sep=",")



### Object for Shiny ###
shiny_obj <- obj_all_genes;
detect_rate <- group_rowmeans(shiny_obj@assays$RNA@counts > 0, shiny_obj@meta.data$Subcluster_Manual)
exclude_genes <- apply(detect_rate, 1, max) < 0.005
shiny_obj@meta.data <- shiny_obj@meta.data[,c("nCount_RNA","nFeature_RNA","percent.mt","Phase", "donor", "sample", "donor_sex", "donor_age", "Subcluster_Manual")]
shiny_obj@reductions <- obj@reductions
shiny_obj@assays$RNA@scale.data <- matrix(0)
shiny_obj@assays$RNA@var.features <- c(0)
shiny_obj <- shiny_obj[!exclude_genes,]
expr_mat <- shiny_obj@assays$RNA@data
metadata <- shiny_obj@meta.data
metadata$UMAP1 <- shiny_obj@reductions$umap@cell.embeddings[,1]
metadata$UMAP2 <- shiny_obj@reductions$umap@cell.embeddings[,2]
metadata$cell_colour <- cluster_colours[match(shiny_obj@meta.data$Subcluster_Manual, names(cluster_colours))]

saveRDS(shiny_obj, "ShinyApps/NKT_LatticeObj.rds")
saveRDS(list(expr_mat=expr_mat, metadata=metadata), "ShinyApps/NKT_shiny.rds")


png("NKT_recoloured_UMAP.png", width=8, height=6, units="in", res=300)
DimPlot(obj_full, group.by="Subcluster_Manual")+scale_color_manual(values=type_colour_scheme)
dev.off()

## Marker Plots - Nice Versions
label_genes <- function(x, this_name) { names(x) <- rep(this_name, length(x)); return(x)}


NKT_genes <- c("CD8A", "CD3D", "TRAC", "TRBC2", "TRDC", "GNLY", "GZMB", "GZMA", "CCL5", 
			"NKG7", "FCGR3A", "FGFBP2", "CD8B", "IL7R", "CD74", "HLA-DRB1")
NKT_exhaustion = c("PDCD1", "CTLA4", "LAG3", "CD160", "CD274") 
NKT_Treg = c("ITGAE", "FOXP3")
NKT_genes_dot_old <- c("CD52","NCR1", "CXCR6", "CXCR3", "CD69", "PTPRC", "CD8A", "CD3D", "CD3E", "CD8B", "CD4", "TRAC", "TRBC2", "TRDC", "TRGC1", "TRGC2", 
				"IL32", "IL7R", "LTB", "IL17A", "IL18R1", "CD44", "KLRB1", "KLRC1", "KLRF1", "KLRK1", "CCL4", "CCL5", "NKG7",
				"FCGR3A", "FGFBP2", "IL2RB", "GZMA", "GZMB", "CSF2", "GNLY", "KLRD1",
				"CD74", "CD79A", "CD79B", "HLA-DRB1", "HLA-DRA", "AIF1", "PDCD1",
				Prolif_genes, Contamination_genes, RBC_genes)

NKT_genes_dot <- c("CD3D", "CD3E", "TRAC", "TRBC2", "TRDC", "TRGC1", "TRGC2", 
				"CD8A", "CD8B", "CCL3", "CCL4", "CCL5", "IL7R", "LTB", "KLRB1", "TPT1",
				"KLRC1", "KLRF1", "GZMK", "CMC1", "XCL1", "XCL2", "GZMB", "FCGR3A", "GNLY", "CXCR6", "CD69", "EOMES", "TBX21",
				"CD79A", "CD79B", "CD74", "HLA-DRB1", "HLA-DRA", "HLA-DPA1", "HLA-DQB1",
				"AIF1", "VSIG4", "LYZ", "VCAN", "CD163", "C1QC", "ITGAM", "ITGAE", "CST3",
				"MKI67", "BRIC5", "TOP2A", "CDK1", "HMGB2", "HBB", "HBA1", "HBA2", "HBD",
				"ALB", "SERPINA1", "APOA1", "FGA" 
				)

# Dendritic Cell marker from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4403526/
cDC1 <- CD45+("PTPRC") HLA-DR+("HLA-DRA") CD141+("THBD") CD123-("IL3RA")CD11c+("ITGAX") CD14-("CD14")
cDC2 <- CD45+("PTPRC") HLA-DR+("HLA-DRA") CD1c+("CD1C") CD123-("IL3RA") CD11c+("ITGAX") CD141-("THBD") CD14+ ("CD14")
pDC <- HLA-DR+("HLA-DRA") CD123+("IL3RA") CD11c-("ITGAX") CD303+("CLEC4C") CD304+ ("NRP1")
Mac <- CD68+("CD68") CD163+("CD163")
Monocytes <- CD14hi("CD14") CD16- ("FCGR3A")
# https://www.frontiersin.org/articles/10.3389/fimmu.2021.641240/full
cDC1 <- CD8A+, CD103+ (ITGAE), B220- (PTPRC), CD11b-, CD1C+
myeloid <- CD8A-, B220-, CD11B+ (ITGAM)

NKT_dendritic <- c("PTPRC", "HLA-DRA", "IL3RA", "ITGAX", "CD14", "THBD", "CD1C", "CLEC4C", "NRP1", "CD68", "CD163", "FCGR3A", "ITGAM", "CD8A", "ITGAE")


DotPlot(obj_full, group.by="Subcluster_Manual", features=c(NKT_dendritic, "IL10", "TGFB"))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

NKT_global_genes <- c(
			label_genes(c("CD3D", "CD3E", "TRAC", "TRBC2"), "abT"), 
			label_genes(c("TRDC", "TRGC1", "TRGC2", "GATA3", "CD7", "CTSW"), "gdT"), 
			label_genes(c("CD8A", "CD8B"), "CD8+"),
			label_genes(c("IL7R", "CD4", "CD40LG", "RCAN3", "LTB","S100A4", "SPOCK2" ), "CD4+ T"),
			label_genes(c("KLRC1", "KLRF1", "GZMK", "CMC1", "XCL1", "CD69", "EOMES"), "lrNK"),
			label_genes(c("XCL2", "GZMB", "FCGR3A", "GNLY","TBX21"), "cNK"),
			label_genes(c("CD74", "HLA-DRB1", "HLA-DRA", "HLA-DPA1", "HLA-DQB1"), "MHCII"),
			label_genes(c("CD79A", "CD79B"), "B"),
			label_genes(c("AIF1", "LYZ", "CD163", "C1QC", "CST3"), "Myeloid"),
			label_genes(c("HBB", "HBA1", "HBA2", "HBD"), "RBC"),
			label_genes(c("MKI67", "TOP2A", "CDK1", "HMGB2"), "CC"),
			label_genes(c("CYP3A4", "CYP2E1", "CYP1A2", "TTR"), "Hep")
			)
require(ggplot2)
png("Figure_NKT_ManualMarkers.png", width=12.5, height=4, units="in", res=300)
DotPlot(obj_full, group.by="Subcluster_Manual", features=NKT_global_genes)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


neutrophils <- DotPlot(obj_all_genes, features=c("FCGR3B", "CXCL8", "CXCR2", "CXCR1", "IFITM2", "CSF3R", "FPR1", "S100A11", "BASP1", "NAMPT", "G0S2"), group.by="Subcluster_Manual")


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

obj@meta.data$Bcell <- cells_AUC@assays@data$AUC["HumanBcell",]
obj@meta.data$Tmem <- cells_AUC@assays@data$AUC["cTmem",]
obj@meta.data$Naive_CD4 <- cells_AUC@assays@data$AUC["NaiveCD4",]
obj@meta.data$Treg <- cells_AUC@assays@data$AUC["Treg",]
obj@meta.data$CD8T <- cells_AUC@assays@data$AUC["CD8T",]
obj@meta.data$Treg <- cells_AUC@assays@data$AUC["Treg",]
obj@meta.data$ResNK <- cells_AUC@assays@data$AUC["HumanResNK",]
obj@meta.data$cNK <- cells_AUC@assays@data$AUC["HumancNK",]

png("NKT_Guilliams_AUCell.png", width=6*3/2, height=8, units="in", res=300)
FeaturePlot(obj, feature=c("cNK", "ResNK", "Treg", "Tmem", "Naive_CD4", "CD8T", "Bcell"))
dev.off()

## ------------ ##


#BioRad
cDC <- c("CD8A", "CD1C", "ITGAX", "ITGAM", "ITGAE", "LY75", "HLA-DRA")
pDC <- c("TLR7", "TLR9", "IL3RA", "LILRA4", "NRP1", "CLEC4C")

DotPlot(obj_all_genes, group.by="Subcluster_Manual", features=c("PTPRC", "HLA-DRA", "IL3RA", "ITGAX", "CD14", "THBD", "CLEC4C", "NRP1", "CD68", "CD163", "FCGR3A")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

obj_all_genes@reductions <- obj@reductions
FeaturePlot(obj_all_genes, features=pDC)
DotPlot(obj_all_genes, features=c(cDC, pDC), group.by="Subcluster_Manual")




## lr vs cNK genes
markers <- FindMarkers(obj_full, group.by="Subcluster_Manual", ident.1="cNKcell", ident.2="lrNKcell", logfc.threshold=-Inf)
DotPlot(obj_full, group.by="Subcluster_Manual", features=rownames(markers)[1:50])+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

write.table(markers, "Raw_cNKcell_vs_lrNKcell_markers.csv", sep=",", row.names=T, col.names=T)

pseudo_detect <- group_rowmeans(obj_full@assays$RNA@counts > 0, obj_full@meta.data$Subcluster_Manual)
pseudo_mean <- group_rowmeans(obj_full@assays$RNA@data, obj_full@meta.data$Subcluster_Manual)

is.clean <- ! colnames(pseudo_detect) %in% c("Debris", "Flush", "Prolif")

lr_l2fc <- log2(pseudo_mean[,"lrNKcell"]/ (apply(pseudo_mean[,is.clean & colnames(pseudo_detect) != "lrNKcell"], 1, max)) )
c_l2fc <- log2(pseudo_mean[,"cNKcell"]/ (apply(pseudo_mean[,is.clean & colnames(pseudo_detect) != "cNKcell"], 1, max)) )
lr_detect <- pseudo_detect[,"lrNKcell"]- apply(pseudo_detect[,is.clean & colnames(pseudo_detect) != "lrNKcell"], 1, max)
c_detect <- pseudo_detect[,"cNKcell"] - apply(pseudo_detect[,is.clean & colnames(pseudo_detect) != "cNKcell"], 1, max)

markers[rownames(markers) %in% names(lr_detect)[lr_detect > 0.1],]

write.table(markers[rownames(markers) %in% names(lr_detect)[lr_detect > 0.1],], "ForLewis_lrNK_Markers.csv", sep=",", row.names=T, col.names=T)
write.table(markers[rownames(markers) %in% names(c_detect)[c_detect > 0.1],], "ForLewis_cNK_Markers.csv", sep=",", row.names=T, col.names=T)

## vs https://www.cell.com/cell-reports/pdfExtended/S2211-1247(19)31372-5

rapid_cycling <- c("PTTG1", "PTPN7", "BCL7C", "TMSB10", "GZMA", "XCL2", "CEBPD", "ENTPD1", "MT-CO3", "HSPD1", "FBXO6", "LAIR2", "CCND2", "GSTP1",
		"ASB2", "CLDND1", "NME1", "GYG1", "GAPDH", "TESC", "LST1", "VIM", "LTB", "CD52", "ITGB7", "NDFIP2", "GZMK", "HLA-DQB1", "MRPL52", "KLRC1",
		"LGALS1", "RCBTB2", "XCL1", "LTA", "C15orf48", "CAPG", "COTL1",
		"KLRC1", "GZMK", "TNFSF10", "ENTPD", "ASB2", "CAPG", "COTL1", "NME1", "CCND2", "FBXO6", "BCL7C", "CEBPD", "HSPD1", "FBXO6", "TESC")
slow_cycling <- c("FGFBP2", "CCL5", "MALAT1", "S100A4", "KLRF1", "FCGR3A", "LITAF", "KLRD1", "CST7", "CMC1", "BTG1", "PTPRC", "TIMP1", "KLRB1", "ZFP36L2",
			"MYOM2", "PLAC8", "HLA-E", "ANXA1", "PTGER2", "S100A6", "EMP3", "NFKBIA", "APMAP", 
			"LITAF", "ZFP36L2", "S100A4", "TIMP1", "FGFBP2", "FCGR3A", "PTPRC", "CCL5", "KLRF1")
DotPlot(obj_full, group.by="Subcluster_Manual", features=c(label_genes(unique(rapid_cycling), "rapid"), label_genes(unique(slow_cycling), "slow")))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

sig_de <- markers[markers$p_val_adj < 0.05,]
scatter_de_plot_x <- sig_de$avg_log2FC # L2FC
scatter_de_plot_y <- sig_de$pct.1-sig_de$pct.2 # diff in detection

png("NKT_cNK_lrNK_vs_CellPaperCycling.png", width=4, height=4, units="in", res=300)
par(mar=c(4,4,1,1))
plot(scatter_de_plot_x, scatter_de_plot_y, pch=16, col="grey50", xlab="Log2 Fold Change (cNK/lrNK)", ylab="Change in Detection Rate (cNK-lrNK)")
points(scatter_de_plot_x[rownames(sig_de) %in% rapid_cycling], scatter_de_plot_y[rownames(sig_de) %in% rapid_cycling], pch=16, col="red", cex=1.5)
points(scatter_de_plot_x[rownames(sig_de) %in% slow_cycling], scatter_de_plot_y[rownames(sig_de) %in% slow_cycling], pch=16, col="blue", cex=1.5)
abline(v=0, h=0, lty=2)
top_genes <- c("CMC1", "S100A4", "FCGR3A", "GZMK", "XCL1", "GSTP1", "XCL2", "FGFBP2")
indeces <- match(top_genes, rownames(sig_de))
text(scatter_de_plot_x[indeces], scatter_de_plot_y[indeces], top_genes, pos=c(4, 4, 3, 4, 4, 2, 4, 3), cex=0.75, offset=0.3)
legend("bottomright", c("rapid", "slow"), pch=16, pt.cex=1.5, cex=1, col=c("red", "blue"), bty="n")
dev.off()

png("NKT_cNK_lrNK_vs_CellPaperCycling_resized.png", width=8, height=4, units="in", res=300)
par(mar=c(4,4,1,1))
plot(scatter_de_plot_x, scatter_de_plot_y, pch=16, col="grey50", xlab="Log2 Fold Change (cNK/lrNK)", ylab="Change in Detection Rate (cNK-lrNK)")
points(scatter_de_plot_x[rownames(sig_de) %in% rapid_cycling], scatter_de_plot_y[rownames(sig_de) %in% rapid_cycling], pch=16, col="red", cex=1.5)
points(scatter_de_plot_x[rownames(sig_de) %in% slow_cycling], scatter_de_plot_y[rownames(sig_de) %in% slow_cycling], pch=16, col="blue", cex=1.5)
abline(v=0, h=0, lty=2)
top_genes <- c("CMC1", "S100A4", "FCGR3A", "GZMK", "XCL1", "GSTP1", "XCL2", "FGFBP2")
indeces <- match(top_genes, rownames(sig_de))
text(scatter_de_plot_x[indeces], scatter_de_plot_y[indeces], top_genes, pos=c(4, 4, 3, 4, 4, 2, 4, 3), cex=0.75, offset=0.3)
legend("bottomright", c("rapid", "slow"), pch=16, pt.cex=1.5, cex=1, col=c("red", "blue"), bty="n")
dev.off()


# Functional enrichments? for sig DE genes?

## Indi Var Barplot##
require(Seurat)
require(cluster)
require(gplots)
set.seed(4302)

obj_clean <- obj[,! (obj@meta.data[,"Subcluster_Manual"] %in% c("Flush", "Debris"))]

donor_col="donor"
cluster_col="Subcluster_Manual"
min.cells=10
ns <- table(obj_clean@meta.data[,donor_col])
mat <- get_pseudobulk_means(obj_clean@assays$RNA@data, factor(obj_clean@meta.data[,cluster_col]), obj_clean@meta.data[,donor_col])
d <- as.dist(1-cor(mat, method="pearson"))
id <- factor(unlist(lapply(strsplit(colnames(mat), "_"), function(x){x[[1]]})))
scores <- c();
stderr <- c();
for (i in unique(id)) {
	scores <- c(scores, mean(as.matrix(1-d)[id==i, id==i]))
      stderr <- c(stderr, sd(as.matrix(1-d)[id==i, id==i])/sqrt(sum(id==i)^2)) # Fix 25 May 2021
}

png("NKT_recoloured_crosscorrelation_barplot.png", width=8*0.5, height=5*0.5, units="in", res=300)
par(mar=c(4,6,1,1))
bar_loc <- barplot(scores, name=unique(id), xlab="Cross Donor Correlation (average)", xlim=c(0,1), las=1, col=cluster_colours, horiz=T)
arrows(scores, bar_loc,  scores+2*stderr, bar_loc, angle=90, len=0)
#legend("topright", lty=1, c("95% CI"), bty="n")
dev.off()

## Frequency Plots

obj <- obj[,obj$Subcluster_Manual != "Flush"]
obj <- obj[,obj$Subcluster_Manual != "Hepato"]
obj <- obj[,obj$Subcluster_Manual != "HepDoublets"]

obj@meta.data$Subcluster_Manual <- manual_cluster_anno[obj@meta.data[,cluster_col]]
png("NKT_Sex_Proportions.png", width=6, height=8, units="in", res=300)
barplot_freq(obj, cluster_col="Subcluster_Manual", colours=colours_sex, metadata="donor_sex")
dev.off()
png("NKT_Age_Proportions.png", width=6, height=8, units="in", res=300)
barplot_freq(obj, cluster_col="Subcluster_Manual", colours=colours_age, metadata="donor_age_group")
dev.off()
png("NKT_Rejection_Proportions.png", width=6, height=8, units="in", res=300)
barplot_freq(obj, cluster_col="Subcluster_Manual", colours=colours_reject, metadata="trans.rejected")
dev.off()

obj@meta.data$Subcluster_Manual <- factor(obj@meta.data$Subcluster_Manual, levels=cluster_order_levels[cluster_order_levels %in% unique(obj@meta.data$Subcluster_Manual)])

png("NKT_donor_Frequency.png", width=8*0.75, height=5*0.75, units="in", res=300)
tab <- indivar_freq(obj@meta.data$Subcluster_Manual, obj@meta.data$donor, obj@meta.data)
dev.off()

donor2sex <- apply(table(obj@meta.data$donor, obj@meta.data$donor_sex), 1, function(x){levels(obj@meta.data$donor_sex)[which(x==max(x))]})



indi_expr(obj, cluster_col="Subcluster_Manual", tag="NKT", 
		bar_col=cluster_colours[names(cluster_colours) %in% unique(obj$Subcluster_Manual)], width=12)


## DE
for (i in unique(obj@meta.data$Subcluster_Manual)) {
	de <- FindMarkers(obj_all_genes, group.by="Subcluster_Manual", ident.1=i,  logfc.threshold=-Inf)
	write.table(de, file=paste("SubCluster_DE_Tables/NKT_", i, "_DE.csv", sep=""), sep=",", row.names=T, col.names=T, quote=T)
}



#### DO DE - MM - NEBULA ####

for(assay_type in c("5pr", "3pr")) {

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

write.table(output, paste("NKT_Nebula", assay_type, "DE.csv", sep="_"), sep=",", row.names=FALSE, col.names=TRUE, quote=TRUE)
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

	write.table(output, paste("NKT_Nebula", assay_type, "Subcluster", subtype, "DE.csv", sep="_"), sep=",", row.names=FALSE, col.names=TRUE, quote=TRUE)

}





#### DO DE - MM ####
# MAST + random effect DE #
require(Seurat)
require(MAST)
require(ggplot2)
obj <- obj_full[,obj_full@meta.data$assay_type == "5pr"]

obj_sce <- as.SingleCellExperiment(obj)
tmp <- as.matrix(obj@assays$RNA@scale.data); rownames(tmp) <- rownames(obj); colnames(tmp) <- colnames(obj)
obj_sca <- FromMatrix(tmp, colData(obj_sce), rowData(obj_sce), class="SingleCellAssay")

zlm.out2 <- zlm(~donor_sex + donor_age + (1/sample), obj_sca)
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

saveRDS(list(sex=sex.output, age=age.output), "NKT_MM_DE_5pr.rds")


## Sex ##
set.seed(38201)
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

fix_names <- c("Oxidation", "Xenobiotics", "Lipid Met", "",
	"AA Met", "PPARalpha", "Cytochrome P450", "", "Plasma Lipoprotein", "",
	"Peroxisomal protein", "Scavenger Receptors", "Bile Met", "", "", "Steroid Met",
	"", "heme scavenging", "Peroxisomal lipid", "Sulphur AA", "Nervous system", "TB infection", "Antigen Presentation",
	"DAP12", "Leishmania infection", "MHC I", "Phagocytosis", "Nef", "Axon guidance", "E3 ubiquitin ligase",
	"", "Interferon gamma", "Adaptive Immune", "Infection", "HIV", "Interferon a/b", "", "", "ER-Phagosome",
	"Lymphoid/NonLymphoid interactions")
G <- set.vertex.attribute(richments_de$graph, "name", value=fix_names)

png("NKT_DE_Reactome_Sex_pathways.png", width=8, height=8, units="in", res=300)
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
		file="NKT_sex_de.rds")



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

fix_names <- c(
	"Cytochrome P450", "tRNA", "Xenobiotics", "MT-rRNA", "IGF & IGFBPs", "protein phosphorylation",
	"", "MT-tRNA", "Oxidation", "Complement cascade", "Clotting cascade", "Bile", "Plasma lipoprotein",
	"", "ABC transporters", "NR1H2 & NR1H3", "PDGF", "SLC mmb. transport", "Steroid Met", "Extracellular Matrix",
	"FCGR3A-phagocytosis", "", "Parasite infection", "PCGR-phagocytosis", "TCR signaling", "ROS-detoxification",
	"SRP targetting", "DNA Damage Recognition", "DAP12", "Lymphoid-NonLymphoid cell",
	"EPH-Ephrin", "", "Interferon", "HIV", "", "Interferon gamma",
	"ER-Phagosome", "Neutrophil degraulation", "Translation", "Antigen presentation")
G <- set.vertex.attribute(richments_de$graph, "name", value=fix_names)

png("NKT_DE_Reactome_Rejection_pathways.png", width=8, height=8, units="in", res=300)
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
		file="NKT_rejection_de.rds")

## Age ##
set.seed(3771)
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
	"Xenobiotics", "Oxidation", "Plasma lipoprotein", "Cytochrome P450", "Vitamin Met", "", "SLC-transport",
	"Steroid Met", "Bile", "Transport", "SREBF", "", "", "Clotting Cascade", "SREBP", "Visual phototransduction",
	"NR1H2/NR1H3", "Xenobiotics", "PPARalpha", "", "", "ABC transporters", "protein phosphorylation", 
	"Lipid Met", "AA Met", "Nuclear Receptors", "Translation")
G <- set.vertex.attribute(richments_de$graph, "name", value=fix_names)

png("NKT_DE_Reactome_Age_pathways.png", width=8, height=8, units="in", res=300)
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
		file="NKT_age_de.rds")


## ----- Special Annotation ----- ##
# Monocytes
dendritic_markers = c("CD1C", "CLEC7A", "CLEC6A", "XCR1", "CLEC9A", "THBD", "CLEC4C", "NRP1", "IL3RA", "CD14", "CD209", 
				"F13A1", "CD16", "CXCR1", "SLAN", "SECISBP2L", "FCGR3A", "HLA-DRA", "HLA-DPA1", "HLA-DQB1", "HLA-DPB1")
DotPlot(obj_all_genes, group.by="Core_clusters", features=dendritic_markers)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Dendritic Cells / Monocytes
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7186229/
DC_genes <- c(label_genes(c("CLEC9A", "PTTG1", "HIST1H4C"), "cDC1"),
		  label_genes(c("HLA-DQA1", "HLA-DQB1", "CD1C", "PPA1", "FCER1A"), "cDC2"),
		  label_genes(c("LILRA4", "TCF4", "CCDC50", "IRF8", "BCL11A", "ITM2C", "C12orf75", "IRF7", "SPIB"), "pDC"),
		  label_genes(c("LGALS2", "VIM", "CLEC4E", "CRIP1", "S100A8", "S100A9", "S100A12", "RNASE2", "PLBD1", 
					"CLU", "RETN", "VCAN", "MNDA", "LYZ"), "CD14Mono"),
		  label_genes(c("HLA-DPB1", "HLA-DPA1", "CD74", "HLA-DRB1", "MS4A7", "C1QA", "C1QB", "C1QC", 
					"FCGR3A", "RHOC", "SMIM25", "CSF1R"), "CD16Mono")
		  )
png("NKT_Supple_DCMono_dotplot.png", width=12.5, height=4, units="in", res=300)
DotPlot(obj_all_genes, features=DC_genes, group.by="Subcluster_Manual") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

## ----- Special subclustering  ------ ##
NKT_genes_dot_good <- c(label_genes(c("CD3D", "CD3E", "TRAC", "TRBC2", "TRDC", "TRGC1", "TRGC2", "CD8A", "CD8B", "CCL3", "CCL4", "CCL5", "IL7R"), "Tcells"), 
				label_genes(c("KLRB1", "TPT1", "KLRC1", "KLRF1", "GZMK", "CMC1", "XCL1", "XCL2"), "lrNK"), 
				label_genes(c("GZMB", "FCGR3A", "GNLY"), "cNK"),
				label_genes(c("CD79A", "CD79B", "CD74", "HLA-DRB2", "HLA-DRA"), "Bcells"),
				label_genes(c("AIF1", "VSIG4", "LYZ", "VCAN", "CD163", "C1QC"), "Myeloid"),
				label_genes(c("MKI67", "BRIC5", "TOP2A", "CDK1", "HMGB2"), "Prolif"), 
				label_genes(c("HBB", "HBA1", "HBA2", "HBD"), "Eryth"),
				label_genes(c("ALB", "SERPINA1", "APOA1", "FGA"),"Hepato") 
				)
# "LTB","ITGAM", "ITGAE", "CST3",

# Prolif

require(Seurat)
require(ggplot2)
set.seed(109)
prolif <- obj[,obj@meta.data$Core_clusters == "4"]
prolif <- FindVariableFeatures(prolif, method="vst", nfeatures=1000)
prolif <- RunPCA(prolif, features=VariableFeatures(prolif), ndims=1:20)
prolif <- FindNeighbors(prolif, dims=1:12)
prolif <- FindClusters(prolif, resolution=0.5)
png("SubSubcluster_NKT_4_umap.png", width=5.5, height=5.5, units="in", res=300)
DimPlot(prolif)
dev.off();

png("NKT_ProlifAnno_Dotplot.png", width=11, height=4, units="in", res=300)
DotPlot(prolif, features=NKT_genes_dot_good, group.by="seurat_clusters")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


set.seed(819)
de_out <- list();
for (c in levels(prolif$seurat_clusters)) {
	de <- FindMarkers(prolif, ident.1=c)
	de_out[[c]] <- de
}
saveRDS(list(obj=prolif, de=de_out), file="Proliferating_NKT_SubSubclustering_4.rds")

# Prolif Percent Plot
prolif <- readRDS("Proliferating_NKT_SubSubclustering_4.rds")
prolif <- prolif$obj
prolif_ncells <- table(prolif@meta.data$seurat_clusters)
#obj <- obj_full
full_ncells <- table(obj@meta.data$Subcluster_Manual)

cell_types <- c("T cells", "lrNK", "cNK", "B cells", "Monocytes")

Percent_T <- prolif_ncells[2]/sum(full_ncells[1], full_ncells[2])
Percent_lrNK <- (prolif_ncells[1]+prolif_ncells[3])/full_ncells[3]
Percent_cNK <- prolif_ncells[5]/full_ncells[4]
Percent_B <- prolif_ncells[6]/full_ncells[5]
Percent_Mono <- prolif_ncells[7]/full_ncells[6]

nT <- prolif_ncells[2] + full_ncells[1] + full_ncells[2]
nlrNK <- prolif_ncells[1]+prolif_ncells[3] + full_ncells[3]
ncNK <- prolif_ncells[5]+full_ncells[4]
nB <- prolif_ncells[6]+full_ncells[5]
nMono <- prolif_ncells[7]+full_ncells[6]

rownames(Cell_type_colours) <- Cell_type_colours[,1]


Percents <- c(Percent_T, Percent_lrNK, Percent_cNK, Percent_B, Percent_Mono)
Ns <- c(nT, nlrNK, ncNK, nB, nMono)
png("NKT_Percent_Prolifering.png", width=5, height=3.5*0.7, units="in", res=300)
par(mar=c(4,5,1,1))
bar_loc <- barplot(Percents,
		name=c("T cells", "lrNK", "cNK", "B cells", "Myeloid"), 
		xlab="Percent Proliferating", xlim=c(0,0.25), las=1, 
		col=c(col_IL7R_CD3_Tcells, 
			col_lrNK, 
			col_cNK, 
			col_Bcells, 
			col_Myeloid), horiz=T)
arrows( Percents, bar_loc, Percents+2*( Percents*(1-Percents)/sqrt(Ns)), bar_loc, angle=90, length=0)
legend("bottomright", lty=1, c("95% CI"), bty="n")
dev.off()

# Flush

require(Seurat)
require(ggplot2)
set.seed(109)
flush <- obj[,obj@meta.data$Core_clusters == "9"]
flush <- FindVariableFeatures(flush, method="vst", nfeatures=1000)
flush <- RunPCA(flush, features=VariableFeatures(flush), ndims=1:20)
flush <- FindNeighbors(flush, dims=1:12)
flush <- FindClusters(flush, resolution=0.5)
png("SubSubcluster_NKT_9_umap.png", width=8, height=8, units="in", res=300)
DimPlot(flush)
dev.off();

DotPlot(flush, features=NKT_genes_dot_good, group.by="seurat_clusters")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

set.seed(819)
de_out <- list();
for (c in levels(flush$seurat_clusters)) {
	de <- FindMarkers(flush, ident.1=c)
	de_out[[c]] <- de
}
saveRDS(list(obj=flush, de=de_out), file="Flush_NKT_SubSubclustering_9.rds")

cluster_anno<- c("gdT", "cNK", "Prolif"); names(cluster_anno) <- c("0","1","2")
flush <- RenameIdents(flush, cluster_anno)

tmp <- table(obj@meta.data$Subcluster_Manual, obj@meta.data$sample)
grep("Flush", colnames(tmp))
png("Flush_proportion.png", width=4, height=3, units="in", res=200)
par(mar=c(6,4,1,1))
barplot(rbind(tmp[,14], rowSums(tmp[,-14])), las=2, col=c( type_colour_scheme["Flush"], "grey65") )
legend("topright", fill=c("grey65", type_colour_scheme["Flush"]), c("Caudate", "Flush"), bty="n")
dev.off()

## Flush from other computer:
## Flush ##
png("Flush_summary_stats.png", width=8, height=4.5, units="in", res=300)
par(mfrow=c(1,2))
par(mar=c(6,4,1,1))
is.Flush <- grepl("Flush", obj@meta.data$sample)
tab <- table(is.Flush, obj@meta.data$Subcluster_Manual)
tab <- tab[,match(names(cluster_colours), colnames(tab))]
Prop_Flush <- tab[2,]/colSums(tab)
sdev_Flush <- sqrt(Prop_Flush*(1-Prop_Flush)/colSums(tab))
ats <- barplot(Prop_Flush*100, las=2, ylab="Flush (%)", col=cluster_colours, ylim=c(0,80))
arrows(ats, Prop_Flush*100, ats, (Prop_Flush+sdev_Flush*2)*100, len=0)
legend("topleft", lty=1, c("95% CI"), bty="n")

par(mar=c(0,2,3,2))
Prop_type <- tab[2,]
Prop_type <- Prop_type[match(names(cluster_colours), names(Prop_type))]
pie(Prop_type, col=cluster_colours, main="Flush cells")
#ats <- pie(Prop_type), las=2, xlab="Distribution of Flush cells (%)")
#arrows(ats, Prop_Flush*100, ats, (Prop_Flush+sdev_Flush*2)*100, len=0)
dev.off()

# DE within subclusters
all_de <- list()
n_de <- vector()
for (type in names(cluster_colours)) {
	subset <- obj[,obj@meta.data$Subcluster_Manual == type]
	if (sum(subset@meta.data$sample == "C51_Flush") < 5 || sum(subset@meta.data$sample == "C51") < 5){next;}
	de <- run_wilcox(subset, subset@meta.data$sample, ident.1="C51_Flush", ident.2="C51")
	n_de <- c(n_de, sum(de$q.value < 0.05, na.rm=T))
	all_de[[type]] <- de
}

#Flush DE
tab <- matrix(0, nrow=length(all_de), ncol=2)
for (i in 1:length(all_de)) {
	out <- all_de[[i]]
	up <- out[out$log2fc > 2 & is.finite(out$log2fc) & out$q.value < 0.05,]
	dn <- out[out$log2fc < -2 & is.finite(out$log2fc) & out$q.value < 0.05,]
	tab[i,] <- c(nrow(up), nrow(dn))
}

rownames(tab) <- names(all_de)

png("Flush_vs_C51_DE_log2FC.png", width=6*0.8, height=4*0.8, units="in", res=300)
par(mar=c(6,4,1,1))
barplot(t(tab[c(-8,-9),]), las=2, col=c("red", "blue"), ylab="Genes L2FC > 2")
legend("topright", fill=c("red", "blue"), c("Up", "Down"), bty="n")
dev.off()





## ILC subclustering - IL7R Cluster ##
ILC_genes <- unique(c("IKZF2", "ID2", "TNFSF4", "NCR1", "GATA3","NCAM1", "ITGA4", "CXCR6", "CD69", "KLRC2", "KLRB1", "TBX21", "EOMES", "PTPRC", "FCGR3A", "NCR2", "PTGDR2","KIT", "IL7R", "CD14", "CD2", "CD3D", "ITGAM", "FUT4", "CD14", "CD19", "CD46", "IL4RA"))
ILC_genes <- c( label_genes(c("NCAM1", "NCR1", "ITGAE", "CD69", "EOMES", "IFNG", "PRF1"), "ILC1"), "IL7R",
			label_genes(c("PTGDR2", "KLGR1", "ST2", "IL2RA", "KLRB1", "GATA3", "BCL11B", "IL5", "IL9", "IL13", "AREG", "IL10"), "ILC2"),
			label_genes(c("NCR2", "KIT", "CCR6", "IL22", "IL17A", "IL26", "CSF2", "TNF"), "ILC3"),
			label_genes(c("SELL", "CCR7", "S1PR1", "GZMB", "GZMK", "KLRF1"), "NKcell"),
			label_genes(c("CD3D", "TRAC", "TRBC1", "TRBC2"), "Tcell"))
			
# https://www.cell.com/cell/pdf/S0092-8674(16)30993-X.pdf
ILC_genes2 <- c( label_genes(c("TBX21", "IFNG", "IL21R", "CCL5", "CCL4", "CCL3", "XCL1"), "ILC1"),
			label_genes(c("GATA3", "HES1", "LMO4", "KLF4", "AREG", "CCL1", "CSF2", "IL4", "IL5", "IL13"), "ILC2"),
			label_genes(c("RORC", "TCF7", "FOXS1", "BATF3", "IL22", "CX3CL1", "IL17F", "GPX1"), "ILC3"),
			label_genes(c("GZMA", "HPOX", "EPAS1"), "ILC1-ILC2"),
			label_genes(c("CXCL2", "CXCL3", "ARG1"), "ILC2-ILC3"),
			label_genes(c("VEGFA", "GZMB", "NCR1"), "ILC1-ILC3"))

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7186229/
ILC_genes3 <- c("LST1", "FXYD5", "OTUD5", "SPINK2", "IL1R1", "ALDOC", "IL4I1", "TNFRSF4", "IL7R", "TNFSF13B")

DotPlot(obj_all_genes, features=ILC_genes, group.by="Subcluster_Manual")





## Correlate to Profiles
profiles <- read.table("C:/Users/tandrews/Documents/UHNSonya/ExternalData/SpectrumRegulatoryLandscapeIntestinalILCs_TableS4.csv", sep=",", header=T)
source("C:/Users/tandrews/Documents/UHNSonya/ExternalData/Ensembl_Stuff.R")
gene_names <- profiles[,1]
gene1 <- sapply(strsplit(gene_names, ";"), function(x){x[[1]]})
gene2 <- sapply(strsplit(gene_names, ";"), function(x){if(length(x)> 1) {return(x[[2]])}else{return(c(""))}})
gene1 <- General_Map(gene1, in.org="Mmus", out.org="Hsap", in.name="symbol", out.name="symbol")
gene2 <- General_Map(gene2, in.org="Mmus", out.org="Hsap", in.name="symbol", out.name="symbol")
gene_names_human <- gene1
gene_names_human[gene1 == ""] <- gene2[gene1==""]
non_unique <- unique(gene_names_human[duplicated(gene_names_human)])
gene_names_human[gene_names_human %in% non_unique] <- ""

profiles_h <- profiles[gene_names_human != "",2:15];
rownames(profiles_h) <- gene_names_human[gene_names_human != ""]

top_genes <- sort(names(tail(sort(apply(profiles_h, 1, var)), 500)))
#heatmap(as.matrix(profiles_h[rownames(profiles_h) %in% top_genes,]))

require(qlcMatrix)
use.genes <- top_genes[top_genes %in% rownames(obj_all_genes)]
ref <- as.matrix(profiles_h[match(use.genes, rownames(profiles_h)),])
query <- obj_all_genes@assays$RNA@data[match(use.genes, rownames(obj_all_genes)),]
cor_mat <- corSparse(ref, query)
matched_type <- apply(cor_mat, 2, function(x){ colnames(ref)[x==max(x)]})
obj_all_genes@meta.data$ILC_anno <- matched_type
obj_all_genes@meta.data$ILC_score <- apply(cor_mat, 2, function(x){ max(x)})


obj_all_genes@reductions <- obj@reductions
DimPlot(obj_all_genes, group.by="ILC_anno")
FeaturePlot(obj_all_genes, features="ILC_score")

round(table(obj_all_genes@meta.data$ILC_anno[obj_all_genes@meta.data$Subcluster_Manual == "lrNKcell"])/sum(obj_all_genes@meta.data$Subcluster_Manual == "lrNKcell"), digits=2)
round(table(obj_all_genes@meta.data$ILC_anno[obj_all_genes@meta.data$Subcluster_Manual == "cNKcell"])/sum(obj_all_genes@meta.data$Subcluster_Manual == "lrNKcell"), digits=2)

#subsubcluster2 <- subsubcluster(obj, cluster=2)
#saveRDS(subsubcluster2,  "NKT_Subsubcluster_2.rds")

subsubcluster2 <- readRDS("NKT_Subsubcluster_2.rds")
all_gene_subsubcluster2 <- obj_all_genes[,obj_all_genes@meta.data[,"Core_clusters"]==2]
all_gene_subsubcluster2@meta.data$Subsubcluster <- subsubcluster2@meta.data$seurat_clusters

DotPlot(all_gene_subsubcluster2, group.by="Subsubcluster", features=ILC_genes)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
DimPlot(subsubcluster2, group.by="seurat_clusters")

markers <- FindMarkers(all_gene_subsubcluster2, group.by="Subsubcluster", ident.1="0", test.use="wilcox", logfc.threshold=-Inf)
head(markers[markers[,2] > 0,], 40)

# Top markers: 
# 0 = ?
# 1 = FOS, JUN, DUSP1 - dissociation / immune response
# 2 = APOA1, APOC3, RBP4 - Hepatocyte contamination
# 3 = GNB2L1, SELL, BIRC3, SOCS3 - Interact with endothelial?
# 4 = CD40LG, CXCR6, GZMK, GZMA - inflammatory
# 5 = SLC4A10, KLRB1, CD160, NKG7 - MAIT?
# 6 = CCR6, KLRB1, CD40LG <--- ILC?
# 7 = TRAT1, ALB, APOB, APOA1
DotPlot(all_gene_subsubcluster2, group.by="Subsubcluster", 
		features=c("CD3D", "CD4", "CD8A", "CD8B", "GNB2L1", "SELL", "BIRC3", "SOCS2", "CXCR6", "KLRB1", "CD69", "EOMES", "ITGA4", "GATA3", 
		"ID2", "CD40LG", "CD160", "CCR6", "NKG7", "GZMK", "GZMA", "SLC4A10")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

FeaturePlot(all_gene_subsubcluster2, c("CD40LG", "SLC4A10", "CXCR6", "CD8A", "CD8B"))


## ILC subclustering - lrNK Cluster ##

set.seed(101)
subsubcluster_lrNK <- subsubcluster(obj_all_genes, cluster="lrNKcell", cluster_col="Subcluster_Manual")
subsubcluster_lrNK@meta.data$Subsubcluster <- subsubcluster_lrNK@meta.data$seurat_clusters
saveRDS(subsubcluster_lrNK, "NKT_Subsubcluster_lrNKcells.rds")

subsubcluster_cNK <- subsubcluster(obj_all_genes, cluster="cNKcell", cluster_col="Subcluster_Manual")
subsubcluster_cNK@meta.data$Subsubcluster <- subsubcluster_cNK@meta.data$seurat_clusters
saveRDS(subsubcluster_cNK, "NKT_Subsubcluster_cNKcells.rds")


#subsubcluster_NK <- readRDS("NKT_Subsubcluster_lrNKcells.rds")

png("NKT_lrNK_subsubcluster_ILC_Dotplot.png", width=8, height=3, units="in", res=300)
DotPlot(subsubcluster_lrNK, group.by="Subsubcluster", features=ILC_genes2)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

png("NKT_cNK_subsubcluster_ILC_Dotplot.png", width=8, height=3, units="in", res=300)
DotPlot(subsubcluster_cNK, group.by="Subsubcluster", features=ILC_genes)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


# Profile 2 profile

png("NKT_lrNK_subsubcluster_ILC_UMAP.png", width=4, height=4, units="in", res=300)
DimPlot(subsubcluster_lrNK)
dev.off()

png("NKT_cNK_subsubcluster_ILC_UMAP.png", width=4, height=4, units="in", res=300)
DimPlot(subsubcluster_cNK)
dev.off()

for(i in 0:7) {
	markers <- FindMarkers(subsubcluster_lrNK, group.by="Subsubcluster", ident.1=i, test.use="wilcox", logfc.threshold=-Inf)
	write.table(markers, file=paste("NKT_lrNK_subsubcluster", i, "DE.csv", sep="_"), sep=",")
	print(head(markers[markers[,2] > 0,], 20))
}
# Top markers: 
# 0 = FOS, DUSP1, ZFP36, JUNB, JUN, CD69, TSC22D3, IER2, AREG, FOS - Early response genes
# 1 = hepato
# 2 = GNB2L1, NKG7, KLRC1, ATP5G2, ATP5L, KLRB1, FYB
# 3 = XCL1, ALB, hepato
# 4 = TYMS, COTL1, ACTB, ACTG1, CDT1, TRDC, KIAA0101, GAPDH
# 5 = IL2RB, FGL1, hepato
# 6 = KLRC1, AMBP, FGL1, hepato
# 7 = TRDC, CD3D, CD3G IL32



for(i in 0:7) {
	markers <- FindMarkers(subsubcluster_cNK, group.by="Subsubcluster", ident.1=i, test.use="wilcox", logfc.threshold=-Inf)
	write.table(markers, file=paste("NKT_cNK_subsubcluster", i, "DE.csv", sep="_"), sep=",")
	print(head(markers[markers[,2] > 0,], 20))
}

# Top markers: 
# 0 = JUN, FOS, DUSP, HBB, IGKC
# 1 = APOA1, ALB, APCO3, HP, VTN
# 2 = PTGDS**, S100B**, FGFBP2, MALAT, JUN
# 3 = MYOM2, FGG, APOB, FN1, APOH, FGB
# 4 = LAIR2, RBM3, CD8A, NKG7, NEL, PTPN7
# 5 = PTGDS, ATP5F1E, C3, VTN, APOA1
# 6 = KIR2DL3, WDR44, FAM50A, ACTB, PTAR1
# 7 = TRDC, APOA1, VTN, ALB, TF, FTCD



DotPlot(subsubcluster_lrNK, group.by="Subsubcluster", 
		features=c("FOS", "JUN", "AREG", "IER2", "GNB2L1", "NKG7", "KLRC1", "KLRB1", "XCL1", "ALB", "TYMS", "ACTB", "ACTG1", "TRDC", "CD3D", "CD3G", "IL32", "IL2RB", "FGL1")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Cluster Cor
subsubcluster_cNK_profile <- group_rowmeans(subsubcluster_cNK@assays$RNA@data, subsubcluster_cNK@meta.data$Subsubcluster)
subsubcluster_lrNK_profile <- group_rowmeans(subsubcluster_lrNK@assays$RNA@data, subsubcluster_lrNK@meta.data$Subsubcluster)

use.genes <- top_genes[top_genes %in% rownames(obj_all_genes)]
ref <- as.matrix(profiles_h[match(use.genes, rownames(profiles_h)),])
query <- subsubcluster_cNK_profile[match(use.genes, rownames(subsubcluster_cNK_profile)),]
cor_mat <- Hmisc::rcorr(ref, query)

png("cNK_vs_ILC_profiles.png", width=4, height=6, units="in", res=300)
gplots::heatmap.2(cor_mat$r[1:14, 15:ncol(cor_mat$r)], scale="none", dendrogram="none",
				col=colorRampPalette(c("white", "black"))(10), trace="none", breaks=seq(from=0.30, to=0.5, length=11),
				density.info="none", key.title="", key.xlab="Correlation (r)", key.ylab="", Rowv=FALSE)
dev.off()

use.genes <- top_genes[top_genes %in% rownames(obj_all_genes)]
ref <- as.matrix(profiles_h[match(use.genes, rownames(profiles_h)),])
query <- subsubcluster_lrNK_profile[match(use.genes, rownames(subsubcluster_lrNK_profile)),]
cor_mat <- Hmisc::rcorr(ref, query)

png("lrNK_vs_ILC_profiles.png", width=4, height=6, units="in", res=300)
gplots::heatmap.2(cor_mat$r[1:14, 15:ncol(cor_mat$r)], scale="none", dendrogram="none",
				col=colorRampPalette(c("white", "black"))(10), trace="none", breaks=seq(from=0.30, to=0.5, length=11),
				density.info="none", key.title="", key.xlab="Correlation (r)", key.ylab="", Rowv=FALSE)
dev.off()


## Add TCR data ##
TCR_files <- Sys.glob("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/TCR_BCR/*.rds")
is.MAIT <- rep(FALSE, ncol(subsubcluster2))

for(f in TCR_files) {
	tmp <- readRDS(f)
	potential_MAITs <- tmp$tcr[,"v_gene"] %in% c("TRAV1", "TRAV1-2") & 
					tmp$tcr[,"j_gene"] %in% c("TRAJ33", "TRAJ12", "TRACJ20") & 
					tmp$tcr[,"c_gene"] == "TRAC"
	id <- unlist(strsplit(f, "/")); id <- unlist(strsplit(id[length(id)], "_TCR"))[1]
	barcode <- sub("-1", "", tmp$tcr[potential_MAITs,1])
	is.MAIT[subsubcluster2@meta.data$sample == id & subsubcluster2@meta.data$cell_barcode %in% barcode] <- TRUE
}
subsubcluster2@meta.data$seurat_clusters[is.MAIT]


is.MAIT.allNKT <- rep(FALSE, ncol(obj))

for(f in TCR_files) {
	tmp <- readRDS(f)
	potential_MAITs <- tmp$tcr[,"v_gene"] %in% c("TRAV1", "TRAV1-2") & 
					tmp$tcr[,"j_gene"] %in% c("TRAJ33", "TRAJ12", "TRACJ20") & 
					tmp$tcr[,"c_gene"] == "TRAC"
	id <- unlist(strsplit(f, "/")); id <- unlist(strsplit(id[length(id)], "_TCR"))[1]
	barcode <- sub("-1", "", tmp$tcr[potential_MAITs,1])
	is.MAIT.allNKT[obj@meta.data$sample == id & obj@meta.data$cell_barcode %in% barcode] <- TRUE

	print(table(tmp$tcr$raw_clonotype_id))
	print(dim(tmp$tcr))

}
n_TCRs_total = 153+6457+5511+46

png("NKT_MAIT_from_TCR_by_Subcluster.png", width=3, height=4, units="in", res=150)
barplot(table(obj@meta.data$Subcluster_Manual[is.MAIT]), las=2, col=cluster_colours)
dev.off()



## NK Markers ##

extract_pathway_dat <- function(res, this_pathway, this_name) {
	genes <- unlist(res$rich[unlist(res$rich[,1]) == this_pathway,8])
	names(genes) <- rep(this_name, length(genes));

	bar_point <- unlist(abs(log10(res$rich[unlist(res$rich[,1]) == this_pathway,"padj"]))); names(bar_point)<- this_name;
	return(list(genes=genes, point=bar_point))
}


# liver resident vs typical NK cells.
cluster_col = "Subcluster_Manual"
ident1 = "lrNKcell"
ident2 = "cNKcell"

mk_lrNK <- run_wilcox(obj, obj@meta.data[,cluster_col], ident.1=ident1)
mk_tyNK <- run_wilcox(obj, obj@meta.data[,cluster_col], ident.1=ident2)
mk_lrNK_vs_tyNK <- run_wilcox(obj, obj@meta.data[,cluster_col], ident.1=ident1, ident.2=ident2)

head(mk_lrNK_vs_tyNK[order(mk_lrNK_vs_tyNK$log2fc, decreasing=T),], 20)

mk_lrNK_vs_tyNK <- mk_lrNK_vs_tyNK[order(mk_lrNK_vs_tyNK$log2fc, decreasing=T),]

source("C:/Users/tandrews/Documents/UHNSonya/FGSEA_enrichment_script.R")

saveRDS(mk_lrNK_vs_tyNK, file="mk_lrNK_vs_tyNK.rds")
mk_lrNK_vs_tyNK <- readRDS("mk_lrNK_vs_tyNK.rds")

tmp <- mk_lrNK_vs_tyNK[,"log2fc"]
names(tmp) <- rownames(mk_lrNK_vs_tyNK)
tmp <- tmp[!grepl("^RP[LS]", names(tmp))]
tmp[tmp == Inf] <- max(tmp[is.finite(tmp)])+1
tmp[tmp == -Inf] <- min(tmp[is.finite(tmp)])-1
tmp <- tmp[!is.na(tmp)]
res_immune2 <- do_fgsea(tmp, pathway=immune_path)
res_kegg <- do_fgsea(tmp, pathway=BaderKegg)
res_react <- do_fgsea(tmp, pathway=reactome)
res_hallmark <- do_fgsea(tmp, pathway=Hallmark_path)

genes1 <- unlist(res_kegg$rich[which(res_kegg$rich[,1] == "Natural killer cell mediated cytotoxicity%KEGG%hsa04650"),8])
DotPlot(obj, features=genes1, group.by="Subcluster_Manual")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


require(igraph)
names(V(res_immune2$graph))

fix_names <- c(
	"", "DN PD1 lo", "DN PD1 low CD8+Tcell", "UP IL4 Stim Mac", "", "UP Infection", "", "IL7R+CD8+Tcell", "", "DN PD1 hi", 
	"", "DN ACD Stim CD4+Tcell", "", "Chronic Infection CD8+Tcell", "CIITA downstream DC", "UP CD11C low DecidualMac", "", "DN ACD Stim CD4+Tcell", "", "UP LPS Stim splenocytes",
	"", "NKcell", "", "", "UP TFH CD4+Tcell", "CD4+Tcell", "UP IL4+PPARG Stim Mac", "Tcell", "DN LPS Stim MEF", "DN CD11C low DecidualMac", 
	"NKcell", "NKcell", "DN IL4 Stim Mac", "UP PD1 low CD8+Tcell", "")

G <- set.vertex.attribute(res_immune2$graph, "name", value=fix_names)

png("lrNK_vs_cNK_fgsea_immune_pathways.png", width=8, height=8, units="in", res=300)
set.seed(2910)
plot(G, vertex.color=res_immune2$vertex_col, vertex.size=res_immune2$vertex_size)
dev.off()

# Key pathways: "GSE26495_PD1HIGH_VS_PD1LOW_CD8_TCELL_DN" vs "GSE26495_PD1HIGH_VS_PD1LOW_CD8_TCELL_UP"  
#			"GSE22342_CD11C_HIGH_VS_LOW_DECIDUAL_MACROPHAGES_DN" vs "GSE22342_CD11C_HIGH_VS_LOW_DECIDUAL_MACROPHAGES_UP" 
#			"GSE25123_CTRL_VS_IL4_STIM_MACROPHAGE_DN" vs "GSE25123_CTRL_VS_IL4_STIM_MACROPHAGE_UP"  

genes1 <- unlist(res_immune2$rich[which(res_immune2$rich[,1] == "GSE26495_PD1HIGH_VS_PD1LOW_CD8_TCELL_DN"),8])
genes2 <- unlist(res_immune2$rich[which(res_immune2$rich[,1] == "GSE26495_PD1HIGH_VS_PD1LOW_CD8_TCELL_UP"),8])
genes3 <- unlist(res_immune2$rich[which(res_immune2$rich[,1] == "GSE22342_CD11C_HIGH_VS_LOW_DECIDUAL_MACROPHAGES_DN"),8])
genes4 <- unlist(res_immune2$rich[which(res_immune2$rich[,1] == "GSE22342_CD11C_HIGH_VS_LOW_DECIDUAL_MACROPHAGES_UP"),8])

DotPlot(obj, features=unique(c(genes1, genes2, genes3, genes4)), group.by="Subcluster_Manual")
DotPlot(obj, features=unique(genes1), group.by="Subcluster_Manual")

saveRDS(res_immune2, "lrNK_vs_cNK_immunepathway_rich.rds")

res_immune2 <- readRDS("lrNK_vs_cNK_immunepathway_rich.rds")

# Key pathways
my_key_pathways <- c();

## NK & PD-1 story

interferons <- read.delim("C:/Users/tandrews/Documents/UHNSonya/ExternalData/GeneCards_interferons.csv", sep=",", header=T)
PD1_direct <- c("PDCD1", "CD274", "PDCD1LG2")
PD1_pathway <- res_immune2$rich[grep("PD1", unlist(res_immune2$rich[,1])),]
PD1_hi_vs_lo_UP <- unlist(PD1_pathway[4,8]); names(PD1_hi_vs_lo_UP) <- rep("PD-1 UP", length(PD1_hi_vs_lo_UP))
PD1_hi_vs_lo_DN <- unlist(PD1_pathway[2,8]); names(PD1_hi_vs_lo_DN) <- rep("PD-1 DN", length(PD1_hi_vs_lo_DN))

my_key_pathways <- c(my_key_pathways, PD1_pathway[2,1], PD1_pathway[4,1])

PD1_related <- c(PD1_direct, PD1_hi_vs_lo_UP, PD1_hi_vs_lo_DN);
detection_rate <- group_rowmeans(obj_all_genes@assays$RNA@counts, obj_all_genes@meta.data$Subcluster_Manual)
detection_rate_PD1 <- detection_rate[match(PD1_related, rownames(detection_rate)),]
PD1_related_good <- unlist(PD1_related[detection_rate_PD1[,3] > 0.05 | detection_rate_PD1[,4] > 0.05])

require(ggplot2)
png("NKT_PD1_Dotplot.png", width=8.5, height=5, units="in", res=300)
DotPlot(obj, features=PD1_related_good, group.by="Subcluster_Manual")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

## NK & CD11C
CD11C_direct <- c("ITGAX")
CD11C_pathway <- res_immune2$rich[grep("CD11C", unlist(res_immune2$rich[,1])),]
CD11C_hi_vs_lo_UP <- unlist(CD11C_pathway[1,8]); names(CD11C_hi_vs_lo_UP) <- rep("CD11C UP", length(CD11C_hi_vs_lo_UP))
CD11C_hi_vs_lo_DN <- unlist(CD11C_pathway[2,8]); names(CD11C_hi_vs_lo_DN) <- rep("CD11C DN", length(CD11C_hi_vs_lo_DN))

my_key_pathways <- c(my_key_pathways, CD11C_pathway[1,1], CD11C_pathway[2,1])

CD11C_related <- c(CD11C_direct, CD11C_hi_vs_lo_UP, CD11C_hi_vs_lo_DN);
detection_rate <- group_rowmeans(obj_all_genes@assays$RNA@counts, obj_all_genes@meta.data$Subcluster_Manual)
detection_rate_CD11C<- detection_rate[match(CD11C_related, rownames(detection_rate)),]
CD11C_related_good <- unlist(CD11C_related[detection_rate_CD11C[,3] > 0.05 | detection_rate_CD11C[,4] > 0.05])

require(ggplot2)
png("NKT_CD11C_Dotplot.png", width=8.5, height=5, units="in", res=300)
DotPlot(obj, features=CD11C_related_good, group.by="Subcluster_Manual")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


## NK & IL4
IL4_direct <- c("IL4")
IL4_pathway <- res_immune2$rich[grep("IL4", unlist(res_immune2$rich[,1])),]
IL4_hi_vs_lo_UP <- unlist(IL4_pathway[1,8]); names(IL4_hi_vs_lo_UP) <- rep("IL4 UP", length(IL4_hi_vs_lo_UP))
IL4_hi_vs_lo_DN <- unlist(IL4_pathway[10,8]); names(IL4_hi_vs_lo_DN) <- rep("IL4 DN", length(IL4_hi_vs_lo_DN))

my_key_pathways <- c(my_key_pathways, IL4_pathway[1,1], IL4_pathway[10,1])

IL4_related <- c(IL4_direct, IL4_hi_vs_lo_UP, IL4_hi_vs_lo_DN);
detection_rate <- group_rowmeans(obj_all_genes@assays$RNA@counts, obj_all_genes@meta.data$Subcluster_Manual)
detection_rate_IL4 <- detection_rate[match(IL4_related, rownames(detection_rate)),]
IL4_related_good <- unlist(IL4_related[detection_rate_IL4[,3] > 0.1 | detection_rate_IL4[,4] > 0.1])

require(ggplot2)
png("NKT_IL4_Dotplot.png", width=10.5, height=5, units="in", res=300)
DotPlot(obj, features=IL4_related_good, group.by="Subcluster_Manual")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

my_key_pathways <- unique(unlist(my_key_pathways))
pvals <- abs(log10(res_immune2$rich[match(my_key_pathways, unlist(res_immune2$rich[,1])),3]))
direction <- sign(res_immune2$rich[match(my_key_pathways, unlist(res_immune2$rich[,1])),6])
names <- c("PD1- CD8+Tcell", "PD1+ CD8+Tcell", 
		"CD11C+ Macrophage", "CD11C- Macropahge", 
		"IL4+ Macrophage", "IL4- Macrophage")

dat_table <- cbind(names, pvals*direction)
dat_table <- dat_table[c(2,1,3,4,5,6),]

png("NKT_key_pathways_barplot.png", width=8*0.7, height=4.5*0.7, units="in", res=300)
par(mar=c(4,9,1,1))
barplot(unlist(abs(dat_table[,2])), col=c(col_lrNK, col_cNK, col_cNK, col_lrNK, col_cNK, col_lrNK), names=unlist(dat_table[,1]), 
		xlim=c(0, 7), xlab="log10(p value)", horiz=T, las=2)
legend("topright", fill=c(col_lrNK, col_cNK), c("lrNK", "cNK"), bty="n", ncol=2)
dev.off()

# ------------------- DE Vis ------------- #
source("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/IndividualVariability.R")
de_sex <- readRDS("NKT_sex_de.rds")
de_age <- readRDS("NKT_age_de.rds")
de_rej <- readRDS("NKT_rejection_de.rds")

obj <- readRDS("NKT_varimax_Subcluster.rds")

type ="reject"

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

png(paste("NKT", up, "Hi_genes.png", sep="_"), width=8, height=6, units="in", res=150)

par(mfrow=c(2,4))
for(g in top_genes) {
	indivar_DE_vis(obj@assays$RNA@scale.data[g,], obj@meta.data$donor, obj@meta.data, type=type)
#	indivar_DE_vis(obj@assays$RNA@data[g,], obj@meta.data$donor, obj@meta.data)

	title(main=g)
}
dev.off()


top_genes <- rownames(de2[order(de2[,5],decreasing=F),])[1:8]

png(paste("NKT", dn, "Hi_genes.png", sep="_"), width=8, height=6, units="in", res=150)

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

png(paste("NKT", up, "Hi_pathways.png", sep="_"), width=8, height=6, units="in", res=150)

par(mfrow=c(2,4))
for(i in 1:min(8, nrow(up_pathways))) {
	genes <- unlist(up_pathways[i, "leadingEdge"])
	indivar_DE_vis(obj@assays$RNA@scale.data[genes,], obj@meta.data$donor, obj@meta.data, type=type)
	title(main=up_pathways$pathway[i])
}

dev.off()

png(paste("NKT", dn, "Hi_pathways.png", sep="_"), width=8, height=6, units="in", res=150)

par(mfrow=c(2,4))
for(i in 1:min(8, nrow(dn_pathways))) {
	genes <- unlist(dn_pathways[i, "leadingEdge"])
	indivar_DE_vis(obj@assays$RNA@scale.data[genes,], obj@meta.data$donor, obj@meta.data, type=type)
	title(main=dn_pathways$pathway[i])
}

dev.off()


