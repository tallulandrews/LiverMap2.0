#### ----- Set Up ----- ####

require("Seurat")
source("~/scripts/LiverMap2.0/Colour_Scheme.R")

dir="/cluster/projects/macparland/TA/PostReview_Disease_vs_Healthy_Map"

args <- commandArgs(trailingOnly=TRUE)

##### Input #####
if (args[1] == "SC") {
	outname <- "SC_Integrated_Map";
} else if (args[1] == "SN") {
	outname <- "SN_Integrated_Map";
} else {
	print("Error:Must provide either SC or SN as argument!")
	exit()
}

#### ----- Annotate with Map2.0 Markers ----- ####

label_genes <- function(x, this_name) { names(x) <- rep(this_name, length(x)); return(x)}

Global_Marker_Table <- read.table("Map2.0_Global_Marker_Table.csv", sep=",", header=TRUE)

# Specific gene plots

key_genes <- c("CES1", "BRD1", "BRD2", "BRD3", "BRD4", "BRDT");

png(paste(outname, "_keyfeatures.png", sep=""), width=12, height =8, units="in", res=300)
FeaturePlot(merged_obj, reduction="umap", features=c(key_genes))
dev.off()

#### ----- Annotate with Map2.0 Markers ----- ####

label_genes <- function(x, this_name) { names(x) <- rep(this_name, length(x)); return(x)}

Cholangiocyte_genes_dot <- c(
		label_genes(c("MUC1", "MUC5B", "MUC3A", "TFF3", "SCGB3A1", "SPINK1", "LGALS2", "PIGR", "SLPI", "LYZ"), "Mucus"),
		label_genes(c("ANXA4", "FXYD2", "RPL3", "EEF1A1", "DEFB1", "GNB2L1", "AMBP", "NEAT1", "HP", "MT2A", "AGXT", "EPCAM", "CLDN1"), "Cholangio"),
		label_genes(c("KRT7", "KRT8", "KRT18", "KRT19"), "KeratinsChol"), 
		label_genes(c("APOC3", "APOA2", "APOA1", "APOC1"), "ApoLipoChol") )

Endo_genes_dot <- c(label_genes(c("FCN2", "FCN3", "STAB1", "STAB2", "CLEC1B", "CLEC4G", "CLEC4M", "CRHBP"), "cvLSEC"),
			label_genes(c("SPARCL1", "TM4SF1", "CLEC14A", "ID1", "IGFBP7", "MGP", "ADRIF", "S100A6", "AQP1"), "ppLSEC"), 
			label_genes(c("RSPO3", "ACKR1", "WNT2"), "cvEndo"), 
			label_genes(c("PODXL", "PLVAP", "CD34"), "Arterial"), 
			label_genes(c("VWF","INMT", "PLAC8","GSN", "RBP7", "RAMP3"), "MiscEndo"),
			label_genes(c("ENG", "PECAM1",  "DNASE1L3", "TIMP3", "LIFR", "C7"), "Endo") )

Macrophage_genes_dot <- c(
		label_genes(c("MARCO", "CD5L", "SLC40A1", "FTL", "CD163", "SEPP1", "C1QC", "C1QB", "C1QA", "CTSB", "HMOX1", "VCAM1"), "NonInfamMac"), 
		label_genes(c("HLA-DRA", "HLA-DPA1", "HLA-DQB1", "HLA-DPB1", "CD74"), "MHCII"), 
		label_genes(c("LYZ", "S100A4", "S100A6", "S100A8", "S100A9", "S100A12", "MNDA", "VCAN", "FCN1"), "InfamMac"),
		label_genes(c("FABP5", "ACP5", "PLD3", "PSAP", "CSTB", "LGMN"), "Lysosome"),
		label_genes(c("FTH1", "CD68", "APOE", "RBP7", "PLTP", "VSIG4", "NINJ1", "IL18",  "CD14","PLAC8", "CD54"), "MiscMac"),
		label_genes(c("FOLR2", "TIMD4", "LYVE1", "FCER1G", "MS4A7", "TIMP1"),"ResidentMac"), 
		label_genes(c("CXCL3", "THBS1", "NAMPT", "CXCL2", "CD83", "IL1B", "AREG", "CCL3","PLAUR", "SRGN"), "ActivatedMac"),
		label_genes(c("LST1", "IFITM3", "AIF1", "COTL1"), "SynapseMac"),
		label_genes(c("FCGR3B", "CXCL8", "CXCR2", "CXCR1", "IFITM2", "CSF3R", "FPR1", "S100A11", "BASP1", "G0S2"), "Neutrophil")
		)
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
			label_genes(c("HBB", "HBA1", "HBA2", "HBD"), "RBC")
			)

Stellate_genes_dot <- c(label_genes(c("COL1A1", "COL1A2", "ACTA2","TAGLN", "MYL9"), "Fiber"),
				label_genes(c("COLEC11", "PTH1R", "RELN", "VIPR1", "HGF", "RBP1"), "Stellate"),
				label_genes(c("CCL2", "SOD2", "SOCS3", "PDGFRA", "FOSB", "FBLN5","JUNB"), "AP1+"),
				label_genes(c("CD9", "CRYAB", "PMEPA1", "RGCC"), "VSMC")
				)

AntiB_genes <- label_genes(c("IGHG2", "IGLL5", "IGHA2", "IGHGP", 
			"IGHM", "IGHG1", "IGHA1", "IGHG3", 
			"IGKC", "IGHG4", "IGLC2", "IGLC3"), "BCR")

Contamination_genes <- label_genes(c("ALB", "SERPINA1", "APOA1", "FGA", "CYP3A5", "CYP2D6", "ASGR1"), "Hepato");
Prolif_genes <- label_genes(c("MKI67", "BIRC5", "TOP2A", "CDK1", "HMGB2", "CDC20", "CDK4"), "Cellcycle"); 

gene_lists <- list(AntiB=AntiB_genes, Stellate=Stellate_genes_dot, Chol=Cholangiocyte_genes_dot, Endo=Endo_genes_dot, Mac=Macrophage_genes_dot, NKT=NKT_global_genes, Hepato=Contamination_genes, Prolif=Prolif_genes)


Global_Marker_Table <- read.table("Map2.0_Global_Marker_Table.csv", sep=",", header=TRUE)

auto_anno_map2.0 <- function(de, title_id="Cluster"){
	markers <- de[de$p_val_adj < 0.05 & de$avg_logFC > 0,]
	scored_marker_table <- markers[match(Global_Marker_Table[,1], rownames(markers)), "avg_logFC"];
	scored_marker_table[is.na(scored_marker_table)] <- 0;
	type_scores <- aggregate(scored_marker_table, by=list(Global_Marker_Table[,2]), mean)
	barplot(type_scores[,2], names=type_scores[,1], main=title_id, las=2)
	return(type_scores)
}

all_cluster_de <- list()
files <- Sys.glob(paste(outname,"*ClusterDE.csv", sep=""))
for (f in files) {
	tab <- read.table(f, sep=",")
	tmp <- unlist(strsplit(f, "_"))
	c <- tmp[length(tmp)-1]
	all_cluster_de[[c]] <- tab
}

all_groups <- unique(as.character(Global_Marker_Table[,2]))
all_scores <- matrix(0, nrow=length(all_groups), ncol=length(all_cluster_de))
rownames(all_scores) <- all_groups;
colnames(all_scores) <- names(all_cluster_de);

pdf(paste(outname, "Map2.0_anno_scores.pdf", sep="_"), width=8, height=4)
for (c in names(all_cluster_de)) {
	anno <- auto_anno_map2.0(all_cluster_de[[c]], paste("Cluster", c))
	anno <- anno[match(rownames(all_scores), anno[,1]),2]
	all_scores[,which(colnames(all_scores)==c)] <- anno	
}
dev.off()

require(gplots)
png(paste(outname, "type_scores_heatmap.png", sep="_"), width=8, height=8, units="in", res=300)
heatmap.2(all_scores, trace="none", scale="column")
dev.off()

pdf(paste(outname, "ref_marker_gene_Dotplot.pdf", sep="_"), width=6, height=10)
for (geneset_i in 1:length(gene_lists)) {
	require(ggplot2)
	#png(paste(outname, names(gene_lists)[geneset_i], "Dotplot.png", sep="_"), 
	#	width=12, height=7, units="in", res=300)
	print(DotPlot(obj, features=gene_lists[[geneset_i]]) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))) + ggtitle(main=names(gene_lists)[geneset_i])
	#dev.off()
}
dev.off()

#Check Neutrophils
require(Seurat)
png(paste(outname, "neutrophil_markers.png", sep="_"), width=8*3/2, height=8, units="in", res=300)
FeaturePlot(obj, features=c("FCGR3B","ITGAX", "IFITM2", "CXCR1", "CXCR2", "CSF3R"))  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# Check Naive T
NaiveTup <- label_genes(c("CCR7", "SELL", "IL7R", "IL2RG"), "NaiveCD4Up") # CD62L, CD127, CD132
NaiveTdn <- label_genes(c("IL2RA", "CD44", "CD69", "HLA-DRA"), "NaiveCD4Dn") #CD25, CD44, CD69, CD45RO, HLA-DRA
png(paste(outname, "NaiveT_markers.png", sep="_"), width=8*3/2, height=8, units="in", res=300)
DotPlot(obj, features=c(NaiveTup, NaiveTdn)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


# Check dendritic
general <- label_genes(c("HLA-DRA", "HLA-DPA1", "HLA-DQB1", "HLA-DPB1"),"Dendritic")
cDC <- label_genes(c("CD8A", "CD1C", "ITGAX", "ITGAM", "ITGAE", "LY75", "CLEC9A", "BATF3", "XCR1", "CLEC10A"), "cDC")
pDC <- label_genes(c("TLR7", "TLR9", "IL3RA", "LILRA4", "NRP1", "CLEC4C", "SCT"), "pDC")
png(paste(outname, "DC_markers.png", sep="_"), width=10*3/2, height=8, units="in", res=300)
DotPlot(obj, features=c(general, cDC, pDC))  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


## Post-Annotation Plots ##
cluster_anno = "cell_type"
require(ggplot2)
png(paste("Suppl", outname, "neutrophil_markers.png", sep="_"), width=6, height=6, units="in", res=300)
DotPlot(obj, features=c("FCGR3B","ITGAX", "IFITM2", "CXCR1", "CXCR2", "CSF3R"), group.by=cluster_anno)  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

png(paste("Suppl", outname, "NaiveT_markers.png", sep="_"), width=6, height=6, units="in", res=300)
DotPlot(obj, features=c(NaiveTup, NaiveTdn), group.by=cluster_anno) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

png(paste("Suppl", outname, "DC_markers.png", sep="_"), width=8*4/3, height=6, units="in", res=300)
DotPlot(obj, features=c(general, cDC, pDC), group.by=cluster_anno)  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


Edge_genes <- c("APOE", "CES1", "MT1G", "AKR1B10", "APOA2", "ANGPTL8", "SQSTM1", "IL32", "SOD2", "HLA-A", "FN1", "P4HB")
png(paste(outname, "edge_markers.png", sep="_"),  width=10*3/2, height=8, units="in", res=300)
print(FeaturePlot(obj, features=Edge_genes) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()


edge_genes <- c("AKR1B10", "TIMP1", "COL1A1", "HLA-A", "MT1G",
                "MT1H", "THY1", "FASN", "COL4A1", "GPNMB",
                "IL32", "LCN2", "COL4A2", "TXN", "SQSTM1", "PLA2G2A")
png(paste(outname, "Edge_signature.png", sep="_"), width=10*1.2*2, height=10*2, units="in", res=150)
FeaturePlot(obj, features=edge_genes)  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
png(paste(outname, "Edge_signature_dotplot.png", sep="_"), width=8, height=8, units="in", res=150)
DotPlot(obj, features=edge_genes, group.by=cluster_anno)  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

Flow_genes <- c("MRC1", "PTPRC", "HLA-DRA", "LYZ", "CD68")
png(paste(outname, "Flowgenes_dotplot.png", sep="_"), width=8, height=8, units="in", res=150)
DotPlot(obj, features=Flow_genes, group.by=cluster_anno)  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# Immune cell population markers
genes=c("EOMES", "CMC1", "KLRG1", "KLRB1", "GZMK", "GZMB", "XCL2", "XCL1", "CD8A", "CD8B", "CD3D", "CD3E", "TRAC", "TRDC", "TRGC2", "IL7R", "IGHM", "IGHD", "IGHA1", "IGKC", "IL2RA", "CTLA4", "FOXP3", "TIGIT", "GATA3", "IL32")
#reordering = c("AntiB", "MatB", "MatB/CD4+", "CD4+T", "Stellate", "ppLSEC", "cvLSEC","Fibroblast", "Kuppfer", "LAM-like", "Monocyte", "Neutrophil", "CD8+T", "gdT", "NKT", "lrNK", "cNK", "Cholangiocyte", )
reordering = c("AntiB", "MatB", "MatB/CD4+", "CD4+T", "Treg", "gdT", "CD8+T", "NKT", "lrNK", "cNK", "Stellate", "ppLSEC", "cvLSEC","Fibroblast", "Cholangiocyte", "Kupffer", "LAM-like", "Monocyte-derived", "Neutrophil", "cDC", "pDC", "Prolif", "RBC", "C-Hepato", "P-Hepato")
obj@meta.data$cell_type <- factor(obj@meta.data$cell_type, levels=reordering)
png(paste(outname, "Immune_markers.png", sep="_"), width=8, height=6, units="in", res=300)
DotPlot(obj[,obj@meta.data$cell_type %in% c("AntiB", "MatB", "Treg", "CD8+T", "gdT", "NKT", "lrNK", "cNK", "MatB/CD4+", "CD4+T")], features=c(genes), group.by=cluster_anno)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


#### ----- END ----- ####

### TLM Plots ###
obj <- readRDS(paste(outname, "_Annotated.rds", sep=""))
Global_Marker_Table <- read.table("Map2.0_Global_Marker_Table.csv", sep=",", header=TRUE)

require(ggplot2)
new_colour_scheme <- unlist(cell_type_cols)[order(names(cell_type_cols))]
names(new_colour_scheme) <- sort(names(cell_type_cols))
new_colour_scheme <- new_colour_scheme[names(new_colour_scheme) %in% obj@meta.data$cell_type]
png(paste(outname, "TLM_unlabelled_umap.png", sep="_"), width=8*1, height=6*0.8, units="in", res=300)
print(DimPlot(obj, reduction="umap", group.by=cluster_anno, pt.size=.1, label=FALSE) +
        scale_color_manual(values=new_colour_scheme))
dev.off()

png(paste(outname, "TLM_labelled_umap.png", sep="_"), width=8*1, height=6*0.8, units="in", res=300)
print(DimPlot(obj, reduction="umap", group.by=cluster_anno, pt.size=.1, label=TRUE) +
        scale_color_manual(values=new_colour_scheme))

dev.off()

png(paste(outname, "TLM_IFN.png", sep="_"), width=8*1, height=6*0.8, units="in", res=300)
print(VlnPlot(obj, "IFNG", group.by=cluster_anno) +
        scale_color_manual(values=new_colour_scheme))
dev.off()

png(paste(outname, "_TGFB.png", sep="_"), width=8*1, height=6*0.8, units="in", res=300)
print(VlnPlot(obj, c("TGFB1", "TGFB2", "TGFB3"), group.by=cluster_anno) +
        scale_color_manual(values=new_colour_scheme))
dev.off()


### ------------ ###


