label_genes <- function(x, this_name) { names(x) <- rep(this_name, length(x)); return(x)}

Contamination_genes <- label_genes(c("ALB", "SERPINA1", "APOA1", "FGA", "CYP3A5", "CYP2D6", "ASGR1"), "Hepato");
Prolif_genes <- label_genes(c("MKI67", "BIRC5", "TOP2A", "CDK1", "HMGB2", "CDC20", "CDK4"), "Cellcycle"; 
RBC_genes <- label_genes(c("HBB", "HBA1", "HBA2", "HBD"), "Eryth");

Dendritic <- c( label_genes(c("HLA-DRA", "HLA-DRB5", "HLA-DRB1"), "MHCII"), 
		    label_genes(c("CLEC7A", "CLEC6A", "ESAM", "CD4", "ITGAM", "CD1C"), "cDC"),
		    label_genes(c("CLEC9A", "XCR1", "THBD"), "crossDC"),
		    label_genes(c("CLEC4C", "NRP1", "IL3RA"), "pDC"),
		    label_genes(c("CD14", "CD209", "F13A1", "CD16", "CXCR1", "SLAN", "SECISBP2L"), "Monocyte"),
		    label_genes(c("ITGAX", "ANPEP", "CD33"), "mDC")
		  )


dendritic_gene_set <- c("HLA-DRA", "HLA-DRB5", "HLA-DRB1"); names(dendritic_gene_set) <- "MHCII"
classic_dendritic<- c("CLEC7A", "CLEC6A", "ESAM", "CD4", "ITGAM", "CD1C"); names(classic_dendritic) <- "cDC"
cross_presenting <- c("CLEC9A", "XCR1", "THBD"); names(cross_presenting) <- "crossDC"
pDC <- c("CLEC4C", "NRP1", "IL3RA"); names(pDC) <- "pDC"
monocyte <- c("CD14", "CD209", "F13A1", "CD16", "CXCR1", "SLAN", "SECISBP2L")
myeloid_MDCs <- c("ITGAX", "ANPEP", "ITGAM", "CD33")

dendritic_gene_neg_set <- c("TRAC", "CD3D", "CD3E", "CD19", "MS4A1", "CD79A", "CD79B", "NCAM1")


## AntiB ##
AntiB_genes <- label_genes(c("IGHG2", "IGLL5", "IGHA2", "IGHGP", 
			"IGHM", "IGHG1", "IGHA1", "IGHG3", 
			"IGKC", "IGHG4", "IGLC2", "IGLC3"), "BCR")
AntiB_genes_dot <- c(AntiB_genes, Prolif_genes, Contamination_genes)
AntiB_cluster_Annotations <- c("Plasmablasts", "IgG+IgK+", "IgA+IgK+", "IgG+IgL+", "IgA+IgL+", "Naive", "Hep/Doublet")

Cholangiocyte_genes <- label_genes(c("EPCAM", "MUC1", "MUC20", "MUC3A", "KRT19", "KRT9",
				"KRT18", "KRT8", "TROP2", "TROP1",
				 "FGFR2", "TM4SF4", "CLDN1",
				"ANXA4", "MKI67", "TFF3"), "Chol")

Cholangiocyte_genes_dot <- c(label_genes(c("MUC1", "MUC5B", "MUC3A", "TFF3", "SCGB3A1", "SPINK1", "LGALS2", "PIGR", "SLPI", "LYZ"), "Mucus"),
		"ANXA4", "FXYD2", "RPL3", "EEF1A1", label_genes(c("KRT7", "KRT8", "KRT18", "KRT19"), "Keratins"), "DEFB1", "GNB2L1", "AMBP", "NEAT1",  
		label_genes(c("APOC3", "APOA2", "APOA1", "APOC1"), "ApoLipo"), "HP", "MT2A", "AGXT", "EPCAM", "CLDN1")


Endo_genes_SN_SC <- c("CTGF", "FCGR2B", "S100A13", "FCN2", "FCN3", "LYVE1", "STAB1", "STAB2", "CLEC1B", "CLEC4G", "CLEC4M", "CRHBP",
				"F8", "CALCRL", "SPARCL1", "TM4SF1", "CLEC14A", "ID1", "IGFBP7", "VWF",
"ENG", "PECAM1", "RAMP3", "INMT", "DNASE1L3", "LIFR", "TIMP3", "C7"
				)
Endo_genes <- c("EPCAM", "PECAM1", "VWF", "CLEC4G", "CLEC4M", 
			"CD34", "CD14", "LYVE1", "RSPO3", "WNT2", "COL1A2", "TFF3",
			"ENG", "STAB2", "CLDN5", "SPARCL1", "RBP7")
Endo_cvLSEC_vs_ppLSEC <- c("CTSD", "CTSL", "CLEC1B", "MS4A6A", 
				"STAB1", "CLEC4G", "CRHBP", "DNASE1L3", "FCN2", "FCN3")
Endo_ppLSEC_vs_cvLSEC_portEndo <- c("MGP", "VIM", "ADIRF", "SPARCL1", "CLU", "S100A6", 
							"CD9", "CLEC14A", "AQP1", "TM4DF1")

Endo_genes_dot <- unique(c(label_genes(c("FCN2", "FCN3", "STAB1", "STAB2", "CLEC1B", "CLEC4G", "CLEC4M", "CRHBP"), "cvLSEC"),
			label_genes(c("SPARCL1", "TM4SF1", "CLEC14A", "ID1", "IGFBP7", "MGP", "ADRIF", "S100A6", "AQP1"), "ppLSEC"), "VWF", 
			label_genes(c("ENG", "PECAM1",  "DNASE1L3", "TIMP3", "LIFR", "C7"), "Endo") "RAMP3",
			label_genes(c("RSPO3", "ACKR1", "WNT2"), "cvEndo"), "INMT", "PLAC8",
			label_genes(c("PODXL", "PLVAP", "CD34"), "Arterial"), "GSN", "RBP7",
			"CCL21", "S100A6", "CST3", "ALB", "APOA1"))

#Endo_genes_SN_SC, "ACKR1", "WNT2", "RSPO3", 
#		Endo_cvLSEC_vs_ppLSEC, Endo_ppLSEC_vs_cvLSEC_portEndo, 
#		"CCL21", "MYL9", "FABP4", "IGFBP5", "COL1A2", "LGALS1", "CRYAB", 
#		"PODXL", "JAG2", "RBP7", "EBF`", "PLVAP", "GSN", "CD34", "ADAM15", "TIMP3"))

Stellate_SC_SN_genes1 <- c("CYP2C9", "CYP2B6", "CYP3A5", "ALDH1A2", "AOX1", "PLIN1", "ADIPOR1", "SOT1", 
			"HGF", "RBP1", "LRAT", "ADAMTSL2", "PDE3B", "ADAMTS2", "SREBF1", "PID1", "LAMB1", "DCN",
			"TGFB1", "PDGFR2", "ACTA2", "PDGFRB", "TGFBR2", "TGFBR1", "SMAD2", "SMAD3", "MYL9", "ACTG1", "VCL", 
			"COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "COL6A3", "IGFBP3", "TAGLN", "SPON2",
			"PPARG", "ADAMTS1", "IGFBP5", "CSRP2", "ACSL4", "IGF1", "CPT1A", "COL5A3", "SLC1A2", "CXCL2", "SMAD4")
Stellate_SC_SN_genes2 <- c("LRAT", "RBP1", "RARA", "RXRA", "RARB", "PLIN2", "ADIPOR1", "DCN", "HGF", "CYP2C9", "CYP2B6", "CYP3A5",
			"ALDH1A2", "AOX1", "ADH4", "ADH1B", "PDE3B", "NR1H4", "PPARG", "EPO", "PLG", "BAMBI", "LAMB1", "PECK1", "FBP1",
		"G6PC", "DGAT1", "PNPL,A2", "PPARGC1A", "HLA-DRA", "HLA-DRB1", "HLA-DMB", "HLA-DPA1", "HLA-DPB1", "CD74", "CD1D",
		"CTSS", "CLIP1", "TAPBP")
Stellate_SC_SN_genes3 <- c("ACTA2", "FN1", "SPARC", "SPRACL1", "TAGLN", "VIM", "MYL9", "ACTG1", "SPON2",
		"VCL", "ITGAV", "COL1A1", "COL1A1", "COL5A2", "COL4A1", "COL3A1", "COL1A2", "COL5A3", "COL4A2", "ACTB", "DDR2", "BGN",
		"GAS6", "SERPINH1", "TNFRSF10A", "CCDC88A", "PFKP", "PKM", "SLC16A3", "LDHA", "GLS", "LIPA", "ALOX5", "HIF1A", "MMP23B",
		"MMP14", "MMP15", "TIMP1", "TIMP2", "TIMP3", "ATG10", "ATG7", "ATG14", "ADAMTS1", "ADMATSL2", "ADAMTS2", "PCOLCE2", 
		"CTSK", "PDE4D", "SRBF1", "SCAP", "SMAD3", "SLF2", "PID1", "LGALS3", "SMO", "PTCH1", "PPARD", "POXO1", "PTEN")
Stellate_SC_SN_genes4 <- c("IFNLR1", "IFNGR1", "IFNAR1", "IFNGR2", "C5AR1", "LBP", "MYD88", "CD14", "TLR1", "TLR2", "IL18R1",
		"IL17RE", "IL17RC", "IL17RB", "IL1RAP", "IL4R", "IL6R", "IL27RA", "IL22RA1", "IL15RA", "IRF7", "IRF1", "CSF1", 
		"TNFSF10", "IL17RA", "IL32", "IL18", "CXCL8", "CCDC88A", "SMAD7", "IL1R1", "NLRP1", "P2RX7", "NLRC5", "CCL2", 
		"CXCL2", "CCL5", "PSTPIP1","CXCL1", "CCL28", "CX3CL1", "CCL4", "CXCL12", "CXCL16", "IL15", "IL7", "IL16", 
		"IL6ST", "CD274", "ICAM1", "VCAM1")
Stellate_SC_SN_genes5 <- c("AEBP1", "SPP1", "IL33", "ICAM1", "VCAM1", "CTGF", "TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", 
		"AGT", "AGTR1", "VEGFA", "VEGFB", "FLT1", "KDR", "IGF1", "PDGFA", "PDGFC", "PDGFD", "PDGFRA", "PDGFRB", "JAG1",
		"FAS", "CEBPA", "MICA", "IGFBP5", "PPARG", "PPARA", "GATA4", "KLF2", "SOCS3", "TP53", "ICAM1", "IL10RB", "IL22RA1")

Stellate_scatter_genes <- c("MYL9", "ACTA2", "TAGLN", "COL1A1", "COL1A2",  "COL6A2")

Stellate_genes <- c("EPCAM", "COL1A1", "CCL2","ACTA2", "MYL9", "DCN", "ACTA2", "VWF", "PODXL", "SOCS3",  "GJA1")
Stellate_fiber_genes <- c("COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A3", "TAGLN", "MYL9", "ACTA2", "DCN")
Stellate_genes_dot <- c("COL1A1", "COL1A2", "ACTA2","TAGLN", "MYL9",
				"CCL2", "SOD2", "SOCS3", 
				"PDGFRB", "PDGFRA", "IFIT3", "EDNRB",  "RSPO3", "DCN", 
				"SNAP25", "KCNIP4", "SYT1", "RBFOX1", "STMN1", "CADM2", "NRXN3", "STXBP5L", "OPCML", "FAM155A", "UCHL1", "NRXN1", "PCDH7", "NTM",
				"COLEC11", "BRPF3", "PTH1R", "PKD2", "RASAL2", "FAM65C", "DNAJB12", "CSEI1L", "ZNF266",
				 "C7", "IGFBP3", "VWF", "PODXL",  "CCL2", "GJA1", Prolif_genes, Contamination_genes)

Stellate_genes_dot <- c("DCN", "MYL9", "ACTG1", "ACTB", "ACTA2", "COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2",
	"TAGLN", "HGF", "RBP1", "LRAT", "ADAMTSL2", "ANGPTL6", "COLEC11", "JUND", "FOSB", "EGR1",
	"CD59", "GSN", "NEAT1", "ANXA2", "CXCL14", "CXCL12", "FBG", "SAA1", "SERPINA1", "MT-CO2", "MT-CO1", 
			"C7", "FBLN5", "CCL2",  "ROSB", "EFEMP1",
			"CRISPLD2",
			"TPM2", "MGP", "MYH11", "TAGLN", "LGALS1", "MFGE8","MYL9",
			"VWF", "MMRN1", "MMRN2", "CD9", "PECAM1", "SOX18", "TM4SF1",
			"SEMA3B", "PMEPA1", "AP1S2", "CRYAB", "S100A4", "IL11RA", "STMN1")

"CD9", "DCN", "VWF", "MMRN1", "MMRN2", 
		"TAGLN", "CDH5", "S100A4", "NEAT1", "GSN", "ANXA2", 
		"CD59", "C7", "FOSB", "IGFBP3", "FBLN5", "CRP", "FGB", 
		"FGA", "SAA1", "SERPINA1", "FGL1", "IGFBP3", "IGFBP7", 
		"RBP1", "CXCL14", "COLEC11", "BGN", "ACTB", "ACTA2", 
		"CXCL12", "CYR61", "COL3A1", "QSOX1", "RELN", "COL1A1", 
		"PTH1R", "DBH", "ANGPTL6")

Macrophage_genes <- c("MARCO", "CD5L", "C1QC", "C1QB", "C1QA", 
				"CD163", "HLA-DRA", "HLA-DPA1", "CD74", "FABP5",
				"PLAUR", "LYZ", "S100A4", "S100A8", "VCAN", "FCN1")

Macrophage_genes_dot <- c("MARCO", "CD5L", "C1QC", "C1QB", "C1QA", "FCGR3A", "FCER1G",
				"CD163", "HLA-DRA", "HLA-DPA1", "HLA-DQB1", "HLA-DPB1", "VSIG4", "CD74", 
				"GPNMB", "ACTP5", "FABP5", "SPP1", "FABP4", "TREM2", "LGALS3", "CTSB", "PSAP", "APOE", "APOC1", "NPC2",
				 "LYZ",  "S100A4", "S100A8", "S100A9", "S100A12", "MNDA", "VCAN", "FCN1", "FTH1", "CD68", 
				"PLAUR", "SRGN", "AREG", "THBS1", "CXCL3", "IL1B", "CCL3", "LILRA5", "ILRB2", "PLAC8", "CD52", 
				 "LST1", Contamination_genes)

Macrophage_genes_dot <- c(
		"MARCO", "CD5L", "LYVE1", "SLC40A1", "FTL", "CD163", "SEPP1", "C1QC", "C1QB", "C1QA", "CTSB", "HMOX1", "VCAM1", 
		"HLA-DRA", "HLA-DPA1", "HLA-DQB1", "HLA-DPB1", "CD74", "VSIG4", 
		"LYZ", "S100A4", "S100A6", "S100A8", "S100A9", "S100A12", "MNDA", "VCAN", "FCN1",
		"FABP5", "ACP5", "PLD3", "FTH1", "CD68", "APOE", "PSAP", "CSTB", "LGMN",
		"RBP7", "FOLR2", "FCER1G", "MS4A7", "TIMP1",
		"JUND", "FOS", "NFKBIA", "ACTG1", "CD14", "CXCL3", "THBS1", "NAMPT", "CXCL2", "CD83", "IL1B",
		"PLAUR", "SRGN", "AREG", "THBS1", "CXCL3", "IL1B", "CCL3",
		"PLAC8", "CD54", "LST1", "IFITM3", "AIF1", "COTL1",
		"DNASE1L3", "FCN2", "CCL14", "FCN3", "SPARC", "CLEC1B", "ENG",
		"ALB", "SERPINA1", "APOA1", "HP", "FGA")

NKT_genes <- c("CD8A", "CD3D", "TRAC", "TRBC2", "TRDC", "GNLY", "GZMB", "GZMA", "CCL5", 
			"NKG7", "FCGR3A", "FGFBP2", "CD8B", "IL7R", "CD74", "HLA-DRB1")
NKT_exhaustion = c("PDCD1") 
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

Hepatocyte1_genes <- c("GLUL", "CYP3A4", "CYP2E1", "CYP1A2", "ADH4", "DCXR", "CYP2A6", "CYP2A7", "HAL", "ALDOB", "APOA1", "CYP1A2", "FABP1", "FGB", "HBB", "HMGB2", "TUBA1B", "SAA1", "MT1H")
Hepatocyte1_genes_dot <- unique(c(Hepatocyte1_genes, Prolif_genes, RBC_genes,  Contamination_genes))

Hepatocyte2_genes_dot <- unique(c("CYP2A7", "CYP1A2", "GLUL", "CYP3A4",  "ADH4", "DCXR", "CYP2A6", "ALDOB", "APOA1", 
					  "FABP1", "FGB", "EPO", "EPOR", Prolif_genes, RBC_genes,  Contamination_genes))

simpson_index <- function(obj, cluster_col="Coarse_clusters") {
	tmp <- table(obj@meta.data[,cluster_col], obj@meta.data$sample)
	prob <- tmp/rowSums(tmp)
	return(rowSums((prob)^2))
}

files <- Sys.glob("*harmony_Subcluster.rds")
for (f in files) {
   print(f)
   obj2 <- readRDS(f);
   print(round(simpson_index(obj2, "Coarse_clusters"), digits=4))
   tmp <- table(obj2@meta.data[,"Coarse_clusters"], obj2@meta.data$sample)
   print( apply(tmp, 1, function(x){colnames(tmp)[which(x==max(x))]}) )
}

#From Map 1.0


require(Seurat)
require(ggplot2)

obj <- readRDS("Stellate_harmony_Subcluster.rds")
obj_all_genes <- readRDS("AllGenes/Stellate_harmony_Subcluster_Allgenes.rds")

FeaturePlot(obj, Cholangiocyte_genes)

FeaturePlot(obj, Stellate_fiber_genes)

DotPlot(obj, features=unique(Stellate_genes_dot), group.by="Fine_clusters")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


DotPlot(obj, features=unique(NKT_genes_dot), group.by="Coarse_clusters", split.by="assay_type", cols=c("red", "blue"))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
DotPlot(obj, features=unique(c(Hepatocyte1_genes_dot)), group.by="Coarse_clusters")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
FeaturePlot(obj, features=c("nCount_RNA", "nFeature_RNA"))
DimPlot(obj, group.by="Fine_clusters", label=T)

DimPlot(obj, group.by="assay_type")

DimPlot(obj, group.by="Phase")
DimPlot(obj, group.by="donor_sex")
DimPlot(obj, group.by="donor_age_group")

tmp1<-FindMarkers(obj, ident.1="0", min.pc=0.001, logfc.threshold=0, group.by="Fine_clusters")
head(tmp1)
tmp2<-FindMarkers(obj, ident.1="5", ident.2="6", min.pc=0.001, logfc.threshold=0, group.by="Coarse_clusters")
tmp3<-FindMarkers(obj, ident.1="6", ident.2="12", min.pc=0.001, logfc.threshold=0, group.by="Coarse_clusters")


de <- readRDS("NKT_subcluster_pseudobulkDE_results.rds")

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


cluster="1"
score <- sort(score_tab[,cluster])
score <- score[!grepl("^RPL", names(score))]
score <- score[!grepl("^RPS", names(score))]
head(score, 10)
tail(score, 10)

tail(score_tab[order(score_tab[,cluster]),], 20)

score2 <- score_tab[,cluster]-apply(score_tab[,colnames(score_tab) != cluster], 1, min)
FeaturePlot(obj, names(sort(score2, decreasing=T))[c(1:12)+12*0])


### Check for specific markers
exhaust_genes = c("LAG3", "CD160", "TBX21", "EOMES", "IL2RA", "CD127", "FOXP3", "CTALA4", "PDCD1", "TIGIT")

source("../../scripts/LiverMap2.0/My_R_Scripts.R")
# % expressed 
prop_expressed = group_rowmeans(obj@assays$RNA@counts > 0, obj@meta.data$Core_clusters)
require(gplots)
prop_expressed[rownames(prop_expressed) %in% exhaust_genes,]
heatmap.2(prop_expressed[rownames(prop_expressed) %in% exhaust_genes,], scale="none", trace="none")






require(fgsea)
immune_path <- gmtPathways("../../ExternalData/MSigDb_immune_signatures_c7.all.v7.1.symbols.gmt")
Hallmark_path <- gmtPathways("../../ExternalData/MSigDbHalmarkPathways.gmt")
MSigAll <- gmtPathways("../../ExternalData/MSigDb_curated_c2.all.v7.1.symbols.gmt")
reactome <- gmtPathways("../../ExternalData/ReactomePathways.gmt")
baderWP <- gmtPathways("../../ExternalData/BaderLab25Aug2020/Human_WikiPathways_August_01_2020_symbol.gmt.txt")
baderIOB <- gmtPathways("../../ExternalData/BaderLab25Aug2020/Human_IOB_August_01_2020_symbol.gmt.txt")
baderMSigdb <- gmtPathways("../../ExternalData/BaderLab25Aug2020/Human_MSigdb_August_01_2020_symbol.gmt.txt")
baderNetPath <- gmtPathways("../../ExternalData/BaderLab25Aug2020/Human_NetPath_August_01_2020_symbol.gmt.txt")

immune_path2 <- immune_path[!grepl("_KO_", names(immune_path))]

res <- fgsea(Hallmark_path, score, minSize=15, maxSize=1000)
res <- res[!is.na(res$pval) & res$padj < 0.05,]
res <- res[order(res$NES),]

#score <- rowSums(score_tab[,c(8,15,13)])

out <- do_fgsea(score, pathways=reactome)

do_fgsea <- function(scored_genes, pathways=MSigAll, fdr=0.05, nmax=20, remove_ribo=TRUE){	
	if (remove_ribo) {
	scored_genes <- scored_genes[!grepl("^RP"
	res <- fgsea(pathways, scored_genes, minSize=15, maxSize=1000, eps=0.00001)
	res <- res[!is.na(res$pval) & res$padj < fdr,]
	res <- res[order(res$NES),]
	if (nrow(res) > nmax) {
		res_pos <- data.frame(res[unlist(res$NES) >0,])
		res_pos <- res_pos[!is.na(unlist(res_pos[,1])),]
		res_neg <- data.frame(res[unlist(res$NES) <0,])
		res_neg <- res_neg[!is.na(unlist(res_neg[,1])),]
		res_pos <- res_pos[order(abs(unlist(res_pos$NES))),]
		res_neg <- res_neg[order(abs(unlist(res_neg$NES))),]
		res <- rbind(res_pos[1:min(nrow(res_pos), nmax),], res_neg[1:min(nrow(res_neg), nmax),])
		res <- res[order(res$NES),]
	}
		
	size <- abs(res$NES)
	colour <- sign(res$NES)
	col_palette <- c("dodgerblue", "grey50", "firebrick")
	gene_lists <- res[,"leadingEdge"]
	sim_mat <- matrix(0, nrow=nrow(res), ncol=nrow(res))
	for (i in 1:nrow(res)) {
		for (j in i:nrow(res)) {
			int <- length(intersect(unlist(gene_lists[i]), unlist(gene_lists[j])))
			uni <- length(union(unlist(gene_lists[i]), unlist(gene_lists[j])))
			sim_mat[i,j] <- int/uni
			sim_mat[j,i] <- int/uni
			colnames(sim_mat) <- unlist(res[,1])
			rownames(sim_mat) <- unlist(res[,1])
		}
	}
	require(igraph)
	G <- simplify(graph_from_adjacency_matrix(sim_mat > 0.1, mode="undirected"))
	plot(G, vertex.color=col_palette[colour+2], vertex.size=size*5, edge.width=2)
	res$cluster <- components(G)$membership
	return(list(rich=res, graph=G, vertex_col = col_palette[colour+2], vertex_size = size*5))
}

























# B-cells
set.seed(0329)
cluster = "0"
score2 <- score_tab[,cluster]-apply(score_tab[,colnames(score_tab) != cluster], 1, max)
out <- do_fgsea(score2, baderMSigdb)
v_names <- c("Coagulation", "Extracellular_Matrix", "Intrinsic_Clotting", 
		"Matrisome", "FOXA2_FOXA3", "Xenobiotic Metabolism",
		"Complement_Pathway", "Integrin2", "Matrisome", "Extracellular_Glycoproteins",
		"Core_Matrisome", "Bile Acid", "TNF_NFkB", "Mitotic_Spindle",
		"Unfolded_Protein_Response", "HDAC_ClassII", "Translation Initiation",
		"Spermatogenesis", "BARD1", "E2F", "Cellcycle", "Myc_targets", "CDK_DNA_replication",
		"MTORC1", "DNA_Repair", "ATR", "Aurora_A", "PLK1", "Oxidative_Phosphorylation",
		"FOXM1", "Proteasome", "MYC_active", "Aurora_B")
		
out$graph_renamed <-  set.vertex.attribute(out$graph, "name", value=v_names)
plot(out$graph_renamed, vertex.color=out$vertex_col, vertex.size=out$vertex_size, edge.width=2)

set.seed(9201)
cluster = "1"
score <- sort(score_tab[,cluster]); score <- score[!grepl("^RPL", names(score))]; score <- score[!grepl("^RPS", names(score))]
score2 <- score_tab[,cluster]-apply(score_tab[,colnames(score_tab) != cluster], 1, max)
out <- do_fgsea(score, reactome)

set.seed(9201)
cluster = "2"
score <- sort(score_tab[,cluster]); score <- score[!grepl("^RPL", names(score))]; score <- score[!grepl("^RPS", names(score))]
score2 <- score_tab[,cluster]-apply(score_tab[,colnames(score_tab) != cluster], 1, max)
out <- do_fgsea(score, reactome)

set.seed(9201)
cluster = "5"
score <- sort(score_tab[,cluster]); score <- score[!grepl("^RPL", names(score))]; score <- score[!grepl("^RPS", names(score))]
score2 <- score_tab[,cluster]-apply(score_tab[,colnames(score_tab) != cluster], 1, max)
out <- do_fgsea(score, reactome)



# This doesn't work
require(GSVA)
#path_score_tab <- gsva(score_tab, reactome)
