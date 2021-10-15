AntiB_genes <- c("IGHG2", "IGLL5", "IGHA2", "IGHGP", 
			"IGHM", "IGHG1", "IGHA1", "IGHG3", 
			"IGKC", "IGHG4", "IGLC2", "IGLC3")

Cholangiocyte_genes <- c("EPCAM", "MUC1", "MUC20", "MUC3A", "KRT19", "KRT9",
				"KRT18", "KRT8", "TROP2", "TROP1",
				 "FGFR2", "TM4SF4", "CLDN1",
				"ANXA4", "MKI67", "TFF3")


Endo_genes <- c("EPCAM", "PECAM1", "VWF", "CLEC4G", "CLEC4M", 
			"CD34", "CD14", "LYVE1", "RSPO3", "WNT2", "COL1A2", "TFF3",
			"ENG", "STAB2", "CLDN5", "SPARCL1", "RBP7")

Stellate_genes <- c("EPCAM", "COL1A1", "DCN", "ACTA2", "VWF", "PODXL", "SOCS3", "CCL2", "GJA1")

Macrophage_genes <- c("MARCO", "CD5L", "C1QC", "C1QB", "C1QA", 
				"CD163", "HLA-DRA", "HLA-DPA1", "CD74", "FABP5",
				"PLAUR", "LYZ", "S100A4", "S100A8", "VCAN", "FCN1")
				
"CXCL8", "LYZ", "MARCO", "S100A8", "S100A9", "TLR2", "TLR4");

#From Map 1.0
Endo_cvLSEC_vs_ppLSEC <- c("CTSD", "CTSL", "CLEC1B", "MS4A6A", 
				"STAB1", "CLEC4G", "CRHBP", "DNASEiL3", "FCN2", "FCN3")
Endo_ppLSEC_vs_cvLSEC_portEndo <- c("MGP", "VIM", "ADIRF", "SPARCL1", "CLU", "S100A6", 
							"CD9", "CLEC14A", "AQP1", "TM4DF1")


require(Seurat)
obj <- readRDS("Macrophage_harmony_Subcluster.rds")
FeaturePlot(obj, features=Macrophage_genes)

DimPlot(obj, group.by="Coarse_clusters", label=T)

DimPlot(obj, group.by="assay_type")

DimPlot(obj, group.by="Phase")
DimPlot(obj, group.by="donor_sex")
DimPlot(obj, group.by="donor_age_group")

de <- readRDS("Macrophage_subcluster_pseudobulkDE_results.rds")

score_tab <- c()

for (cluster in names(de)) {
	#cluster <- "1"
	score <- de[[cluster]]
	score <- sqrt((score[,1]-score[,3])*(score[,2]-score[,4]))*sign(score[,1]-score[,3])
	score[is.na(score)] <- 0;
	if (is.null(dim(score_tab))) {
		score_tab <- cbind(score_tab, score)
	} else {
		score <- score[match(rownames(score_tab), names(score))]
		score_tab <- cbind(score_tab, score)
	}
}
colnames(score_tab) <- names(de)

#heatmap(cor(score_tab), distfun=function(x){as.dist(1-x)}, scale="none")


cluster="3"
score <- sort(score_tab[,cluster])
score <- score[!grepl("^RPL", names(score))]
score <- score[!grepl("^RPS", names(score))]
head(score, 10)
tail(score, 10)

tail(score_tab[order(score_tab[,cluster]),], 20)

score2 <- score_tab[,cluster]-apply(score_tab[,colnames(score_tab) != cluster], 1, max)
FeaturePlot(obj, names(tail(sort(score), 12)))

require(fgsea)
immune_path <- gmtPathways("../../ExternalData/MSigDb_immune_signatures_c7.all.v7.1.symbols.gmt")
Hallmark_path <- gmtPathways("../../ExternalData/MSigDbHalmarkPathways.gmt")
MSigAll <- gmtPathways("../../ExternalData/MSigDb_curated_c2.all.v7.1.symbols.gmt")
reactome <- gmtPathways("../../ExternalData/ReactomePathways.gmt")

res <- fgsea(MSigAll, score, minSize=15, maxSize=1000)
res <- res[!is.na(res$pval) & res$padj < 0.05,]
res <- res[order(res$NES),]