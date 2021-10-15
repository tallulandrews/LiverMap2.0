
dir = "/cluster/projects/macparland/TA/LiverMap2.0/RawData/Analysis/EmptyDrops"
#immune_genesets <- read.table("Immune_pathways_for_sex_differences.csv", sep=",", stringsAsFactors=F, header=TRUE)



Spatial_Portal <- c("ALDOB", "APOA1", "HAMP1", "AGT", "FGB", "CYP2A7", "HAL", "FGA", "SDS")
Spatial_Central <- c("CYP3A4", "CYP2E1", "CYP1A2", "ADH1B", "GLUL", "ADH4", "DCXR", "CES1")
Spatial_RBC <- c("HBA2", "HBA1", "HBB", "HBD", "AHSP", "HBM", "CA1", "SLC44A1", "SLC25A37", "GYPC", "SAA1", "C1QB", "SAA2", "C1QA")
Spatial_cluster12 <- c("CCDC152", "LRRC75A", "MAP3K12", "PHKG1", "ACTB", "ABI2", "AFMID", "GGCX", "MT1X", "TAT")
Spatial_cluster9 <- c("PTGDS", "CCL21", "C7", "FBLN1", "MGP", "TAGLN", "VIM", "IGFBP7", "RNASE1", "S100A6", "TMSB10", "GPX3", "TIMP3", "MYL9")

Map_Portal <- c("CYP3A5", "MLXIPL", "ACSM2B", "SERPINF2", "ITIH1", "MST1")
Map_Central <- c("ANG", "C9", "CHCHD10", "ECHS1", "GATM", "GSTA1", "HSD17B6", "LEAP2", "LINC00844", "MT1F", "MT1M", "SAA4", "SULT2A1", "TAT")
Map_Tcell <- c("CCL5", "NKG7", "HCST", "GZMA", "KLRB1", "CORO1A", "CD7", "CST7", "CCL4", "CD3D", "KLRD1", "PTPRC", "CD48", "CD3E", "TRAC")
Map_cvLSEC <- c("DNASE1L3", "FCN3", "FCN2", "CRHBP", "CLEC4G", "CCL14", "GNG11", "ENG", "STAB1", "CELC1B", "CLEC4M", "LYVE1")
Map_InfMac <- c("LYZ", "AIF1","LST1", "CTSS", "S100A9", "TYROBP", "FCER1G", "S100A4", "S100A6", "S100A8", "VCAN", "CD163")
Map_NonInfMac <- c("C1QC", "C1QB", "C1QA", "CD163", "MS4A7", "AIF1", "HLA-DRA", "FCER1G", "MS4A4", "FCGR3A", "CD68", "MARCO")
Map_NK <- c("STMN1", "TYMS", "HMGB2", "KIAA0101", "H2AFX", "SMC4", "RAC2", "H2AFV", "NKG7", "XCL2", "IL2RB", "GZMA")
Map_AntiB <- c("MZB1", "JCHAIN", "DERL3", "FKBP11", "ITM2C", "ISG20", "SEC11C", "PIM2", "TNFRSF17", "IGKC", "IGHG1", "IGHG3")
Map_Chol <- c("FZYD2", "ANXA4", "CD24", "SPP1", "KRT7", "DEFB1", "AQP1", "TM4SF4", "ELF3", "KRT8", "KRT19", "KRT18", "EPCAM")
Map_Stellate <- c("DCN", "IGFBP7", "BGN", "QSOX1", "COL6A2", "COL6A1", "TAGLN", "COL3A1", "COLEC11", "COL1A1", "MYL9", "COL1A2")
Map_Bcell <- c("HLA-DRA", "CD52", "CD19A", "CD37", "HLA-DPB1", "HLA-DPA1", "HLA-DQB1", "CD19B", "CD74")
Map_Eryth <- c("SNCA", "ALAS2", "AHSP", "CA1", "GYPA","HBD", "SLC4A1", "HBM", "EPB42", "HBA1", "HBA2")


all_gene_lists <- list(Spatial_Portal=Spatial_Portal, Spatial_Central=Spatial_Central, Spatial_RBC=Spatial_RBC, Spatial_cluster12=Spatial_cluster12, Spatial_cluster9=Spatial_cluster9, Map_Portal=Map_Portal, Map_Central=Map_Central, Map_Tcell=Map_Tcell, Map_cvLSEC=Map_cvLSEC, Map_InfMac=Map_InfMac, Map_NonInfMac=Map_NonInfMac, Map_NK=Map_NK, Map_AntiB=Map_AntiB, Map_Chol=Map_Chol, Map_Stellate=Map_Stellate, Map_Bcell=Map_Bcell, Map_Eryth=Map_Eryth)

obj <- readRDS("Merged_EmptyOnly_obj_Map2.2_ImportedClusters_ManualAnno.rds")
obj_3pr <- obj[,obj@meta.data$assay_type == "3pr"]
obj_5pr <- obj[,obj@meta.data$assay_type == "5pr"]

Paul_collaborator <- c(rownames(obj)[grep("WNT", rownames(obj))], rownames(obj)[grep("HIF", rownames(obj))], rownames(obj)[grep("FZD", rownames(obj))])
# Temporary!
all_gene_lists=list(Paul_Custom_list=Paul_collaborator)

for(i in 1:length(all_gene_lists)) {
	require(Seurat)
	require(ggplot2)
	this_name <- names(all_gene_lists)[i]
	gene_set <- unlist(all_gene_lists[[i]])
	png(paste(this_name, "3pr_geneset_dotplot.png", sep="_"), width=8, height=8, units="in", res=300)
	print(DotPlot(obj_3pr, features=gene_set, group.by="Coarse_Manual_Anno") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
	dev.off()
	png(paste(this_name, "5pr_geneset_dotplot.png", sep="_"), width=8, height=8, units="in", res=300)
	print(DotPlot(obj_5pr, features=gene_set, group.by="Coarse_Manual_Anno") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
	dev.off()
}

