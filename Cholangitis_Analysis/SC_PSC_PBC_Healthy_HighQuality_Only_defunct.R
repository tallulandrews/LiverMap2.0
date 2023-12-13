#### ----- Set Up ----- ####

require("Seurat")
source("~/scripts/LiverMap2.0/Colour_Scheme.R")

dir="/cluster/projects/macparland/TA/PostReview_Disease_vs_Healthy_Map"

##### Input #####
PSC = c(
	"/cluster/projects/macparland/TA/PSC/Processed/PSC018_5pr_caudate_EmptyOnly.rds",
	"/cluster/projects/macparland/TA/PSC/Processed/PSC014X_5pr_caudate_EmptyOnly.rds",
	"/cluster/projects/macparland/TA/PSC/Processed/PSC019_Caudate_5pr_EmptyOnly.rds",
	"/cluster/projects/macparland/TA/PSC/Processed/PSC024_Caudate_5pr_EmptyOnly.rds",
	"/cluster/projects/macparland/TA/PSC/Processed/PSC012_SC_5pr_EmptyOnly.rds",
	"/cluster/projects/macparland/TA/PSC/Processed/PSC016_SC_5pr_EmptyOnly.rds"
	)

PBC = c("/cluster/projects/macparland/TA/PBC/Processed/PBC005_Frozen_5prV2_EmptyOnly.rds",
	"/cluster/projects/macparland/TA/PBC/Processed/PBC001_Frozen_5pr_V2_EmptyOnly.rds")

Healthy = c("/cluster/projects/macparland/TA/LiverMap2.0/Cleaned/C58_5pr_EmptyOnly.rds",
	    "/cluster/projects/macparland/TA/LiverMap2.0/Cleaned/C59_5pr_EmptyOnly.rds",
	    "/cluster/projects/macparland/TA/LiverMap2.0/Cleaned/C61_5pr_EmptyOnly.rds",
	    "/cluster/projects/macparland/TA/LiverMap2.0/Cleaned/C63_5pr_reseq_EmptyOnly.rds",
	    "/cluster/projects/macparland/TA/LiverMap2.0/Cleaned/C64_5pr_EmptyOnly.rds",
	    "/cluster/projects/macparland/TA/LiverMap2.0/Cleaned/C70_5pr_reseq_EmptyOnly.rds")

all_samples <- c(PSC, PBC, Healthy)


outname <- "SC_Integrated_Map";

sample_names <- strsplit(all_samples, "/");
sample_names <- sapply(sample_names, function(x){x[length(x)]})
sample_names <- sub("_EmptyOnly.rds", "", sample_names)

## Filters: ##
n_features = 750
n_umi = 1000
scale_size = 2000

#### ----- Merge ----- ####

if (!file.exists(paste(outname, "_Merged.rds", sep=""))) {
set.seed(3921)

merged_obj <- NULL;
merged_scale_mat <- c();
hvgs <- c();

for (i in 1:length(all_samples)) {
	n <- sample_names[i];
	print(n);
	donor_id <- sapply(strsplit(n, "_"), function(x){x[1]})
	obj <- readRDS(all_samples[i])
	print(dim(obj))
	obj <- obj[,
		Matrix::colSums(obj@assays$RNA@counts> 0) > n_features &
		Matrix::colSums(obj@assays$RNA@counts) > n_umi]
	print(dim(obj));
	obj <- NormalizeData(obj, scale.factor=scale_size);
	obj@meta.data$sample <- obj@meta.data$orig.ident
	obj@meta.data$donor <- donor_id
	obj@meta.data$assay_type <- "3pr"
	obj@meta.data$sample_type <- "caudate"
	pheno = "?"
	if (all_samples[i] %in% Healthy) {
		pheno <- "Healthy"
	}
	if (all_samples[i] %in% PSC) {
		pheno <- "PSC"
	}
	if (all_samples[i] %in% PBC) {
		pheno <- "PBC"
	}
	obj@meta.data$Phenotype <- pheno
	if (grepl(all_samples[i], "_5pr_")) {
		obj@meta.data$assay_type <- "5pr"
	}
	if (grepl(all_samples[i], "_PBMCs_")) {
		obj@meta.data$sample_type <- "PBMCs"
	}
	#obj2 <- SCTransform(obj)

	# Merge
	hvgs <- c(hvgs, VariableFeatures(obj));
	if (is.null(merged_obj)) {
		merged_obj <- obj
		merged_scale_mat <- obj@assays$RNA@scale.data
		next;
	} else {
		# All genes
		merged_obj <- merge(merged_obj, obj, add.cell.ids=c("", n), project="SC_PSC_only");
		merged_genes <- rownames(merged_obj);

		merged_scale_mat <- merged_scale_mat[match(merged_genes, rownames(merged_scale_mat)),]
		merged_scale_mat[is.na(merged_scale_mat)] <- 0;

		scale_mat <- obj@assays$RNA@scale.data;
		scale_mat <- scale_mat[match(merged_genes, rownames(scale_mat)),]
		scale_mat[is.na(scale_mat)] <- 0;
		
		merged_scale_mat <- cbind(merged_scale_mat, scale_mat);
		rownames(merged_scale_mat) <- merged_genes;
		print(dim(merged_obj))
		print(dim(merged_scale_mat))
	}

}

hvgs <- unique(hvgs[duplicated(hvgs)])
VariableFeatures(merged_obj) <- hvgs;
colnames(merged_scale_mat) <- colnames(merged_obj)
merged_obj@assays$RNA@scale.data <- merged_scale_mat;
saveRDS(merged_obj, paste(outname, "_Merged.rds", sep=""))

}

#### ----- END ----- ####

#### ----- Integrate with Harmony ----- ####
if (!file.exists(paste(outname, "_Harmony.rds", sep=""))) {
require("harmony")
set.seed(728)

merged_obj <- readRDS(paste(outname, "_Merged.rds", sep=""))

merged_obj <- RunPCA(merged_obj, features=VariableFeatures(merged_obj))
merged_obj <- RunUMAP(merged_obj, seed.use=42, dims=1:30)
merged_obj <- FindNeighbors(merged_obj, reduction="pca", dims=1:30)
merged_obj <- FindClusters(merged_obj, reduction="pca", resolution=0.5, method="igraph")


png(paste(outname, "_indiscaled_sample_umap.png", sep=""), width=9, height =6, units="in", res=300)
DimPlot(merged_obj, reduction="umap", group.by="orig.ident", pt.size=0.1, label=TRUE)
dev.off();

png(paste(outname, "_indiscaled_clusters_umap.png", sep=""), width=9, height =6, units="in", res=300)
DimPlot(merged_obj, reduction="umap", group.by="seurat_clusters", pt.size=0.1, label=TRUE)
dev.off();

png(paste(outname, "_indiscaled_pheno_umap.png", sep=""), width=9, height =6, units="in", res=300)
DimPlot(merged_obj, reduction="umap", group.by="Phenotype", pt.size=0.1, label=TRUE)
dev.off();

merged_obj@meta.data$marker_labs[merged_obj@meta.data$marker_labs == "PortalHep"]  <- "Hepatocyte"
png(paste(outname, "_indiscaled_autoanno_umap.png", sep=""), width=9, height =6, units="in", res=300)
Type_DimPlot(merged_obj,type_col="marker_labs", cluster_col="marker_labs")
dev.off();



# Harmony
set.seed(2910)

merged_obj <- RunHarmony(merged_obj, c("sample"), plot_convergence=TRUE)
merged_obj <- RunUMAP(merged_obj, seed.use=42, dims=1:30, reduction="harmony")
merged_obj <- FindNeighbors(merged_obj, reduction="harmony", dims=1:30)
merged_obj <- FindClusters(merged_obj, reduction="harmony", resolution=2, method="igraph")


png(paste(outname, "_harmony_indiscaled_sample_umap.png", sep=""), width=9, height =6, units="in", res=300)
DimPlot(merged_obj, reduction="umap", group.by="orig.ident", pt.size=0.1, label=TRUE)
dev.off();

png(paste(outname, "_harmony_indiscaled_pheno_umap.png", sep=""), width=9, height =6, units="in", res=300)
DimPlot(merged_obj, reduction="umap", group.by="Phenotype", pt.size=0.1, label=TRUE)
dev.off();

png(paste(outname, "_harmony_indiscaled_clusters_umap.png", sep=""), width=9, height =6, units="in", res=300)
DimPlot(merged_obj, reduction="umap", group.by="seurat_clusters", pt.size=0.1, label=TRUE)
dev.off();

merged_obj@meta.data$marker_labs[merged_obj@meta.data$marker_labs == "PortalHep"]  <- "Hepatocyte"
png(paste(outname, "_harmony_indiscaled_autoanno_umap.png", sep=""), width=9, height =6, units="in", res=300)
Type_DimPlot(merged_obj,type_col="marker_labs", cluster_col="marker_labs")
dev.off();
saveRDS(merged_obj, file=paste(outname, "_Harmony.rds", sep=""))
}

### Integration with CCA ###

if (!file.exists(paste(outname, "_CCA.rds", sep=""))) {
set.seed(728)

merged_obj <- readRDS(paste(outname, "_Merged.rds", sep=""))
obj_list <- SplitObject(merged_obj, split.by="orig.ident")
hvgs <- VariableFeatures(merged_obj)

integration.anchors <- FindIntegrationAnchors(object.list=obj_list, ancho.features=hvgs)
merged_obj <- IntegrateData(anchorset=integration.anchors)

# Cluster
set.seed(2910)

merged_obj <- RunUMAP(merged_obj, seed.use=42, dims=1:30)
merged_obj <- FindNeighbors(merged_obj, dims=1:30)
merged_obj <- FindClusters(merged_obj, resolution=2, method="igraph")

png(paste(outname, "_cca_indiscaled_sample_umap.png", sep=""), width=9, height =6, units="in", res=300)
DimPlot(merged_obj, reduction="umap", group.by="orig.ident", pt.size=0.1, label=TRUE)
dev.off();

png(paste(outname, "_cca_indiscaled_pheno_umap.png", sep=""), width=9, height =6, units="in", res=300)
DimPlot(merged_obj, reduction="umap", group.by="Phenotype", pt.size=0.1, label=TRUE)
dev.off();

png(paste(outname, "_cca_indiscaled_clusters_umap.png", sep=""), width=9, height =6, units="in", res=300)
DimPlot(merged_obj, reduction="umap", group.by="seurat_clusters", pt.size=0.1, label=TRUE)
dev.off();

merged_obj@meta.data$marker_labs[merged_obj@meta.data$marker_labs == "PortalHep"]  <- "Hepatocyte"
png(paste(outname, "_cca_indiscaled_autoanno_umap.png", sep=""), width=9, height =6, units="in", res=300)
Type_DimPlot(merged_obj,type_col="marker_labs", cluster_col="marker_labs")
dev.off();
saveRDS(merged_obj, file=paste(outname, "_CCA.rds", sep=""))
}
exit()

# Specific gene plots

key_genes <- c("CES1", "BRD1", "BRD2", "BRD3", "BRD4", "BRDT");

png(paste(outname, "_keyfeatures.png", sep=""), width=12, height =8, units="in", res=300)
FeaturePlot(merged_obj, reduction="umap", features=c(key_genes))
dev.off()


# Cluster DE
require(Seurat)
set.seed(5128)
obj <- readRDS(paste(outname, "_Harmony.rds", sep=""))

tag = outname

all_cluster_de=list()

for(c in obj@meta.data$seurat_clusters) {
	outfile <- paste(tag, c, "ClusterDE.csv", sep="_")
#	if (file.exists(outfile)) {next;}
        DE <- FindMarkers(obj, ident.1=c, group.by="seurat_clusters", test.use="wilcox")
	all_cluster_de[[c]]<-DE
        write.table(DE, file=outfile, sep=",")
}

#### ----- END ----- ####

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

Global_Marker_Table <- cbind(c(Prolif_genes, Contamination_genes, AntiB_genes, Stellate_genes_dot, 
				NKT_global_genes, Macrophage_genes_dot, Endo_genes_dot, Cholangiocyte_genes_dot), 
			     c(as.character(names(Prolif_genes)), as.character(names( Contamination_genes)), 
				as.character(names(AntiB_genes)), as.character(names(Stellate_genes_dot)),
				as.character(names(NKT_global_genes)), as.character(names(Macrophage_genes_dot)), 
				as.character(names(Endo_genes_dot)), as.character(names(Cholangiocyte_genes_dot)))
			    )
colnames(Global_Marker_Table) <- c("Gene", "Type")
#write.table(Global_Marker_Table, "Map2.0_Global_Marker_Table.csv", sep=",", row.names=FALSE)

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

for (geneset_i in 1:length(gene_lists)) {
	require(ggplot2)
	png(paste(outname, names(gene_lists)[geneset_i], "Dotplot.png", sep="_"), 
		width=12, height=7, units="in", res=300)
	print(DotPlot(obj, features=gene_lists[[geneset_i]]) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
	dev.off()
}

#Check Neutrophils
require(Seurat)
png(paste(outname, "neutrophil_markers.png", sep="_"), width=8*3/2, height=8, units="in", res=300)
FeaturePlot(obj, features=c("FCGR3B","ITGAX", "IFITM2", "CXCR1", "CXCR2", "CSF3R"))  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# Check Naive T
NaiveTup <- c("CCR7", "SELL", "IL7R", "IL2RG") # CD62L, CD127, CD132
NaiveTdn <- c("IL2RA", "CD44", "CD69", "HLA-DRA") #CD25, CD44, CD69, CD45RO, HLA-DRA
png(paste(outname, "NaiveT_markers.png", sep="_"), width=8*3/2, height=8, units="in", res=300)
DotPlot(obj, features=c(NaiveTup, NaiveTdn)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


# Check dendritic
general <- c("HLA-DRA", "HLA-DPA1", "HLA-DQB1", "HLA-DPB1")
cDC <- c("CD8A", "CD1C", "ITGAX", "ITGAM", "ITGAE", "LY75", "CLEC9A", "BATF3", "XCR1", "CLEC10A")
pDC <- c("TLR7", "TLR9", "IL3RA", "LILRA4", "NRP1", "CLEC4C", "SCT")
png(paste(outname, "DC_markers.png", sep="_"), width=10*3/2, height=8, units="in", res=300)
DotPlot(obj, features=c(general, cDC, pDC))  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


# Manual Annotations
cluster2anno <- list( "0"="MatB",
			"1" = "cNK",
			"2" = "CD4+T", # Naive CD4+T
			"3" = "cNK",
			"4" = "Treg", # GATA3+CD4+T = T helper2
			"5" = "P-Hepato",
			"6" = "CD8+T",
			"7" = "lrNK",
			"8" = "Neutrophil",
			"9" = "Monocyte-derived",
			"10" = "MatB",
			"11" = "cvLSEC",
			"12" = "C-Hepato",
			"13" = "gdT",
			"14" = "CD4+T",
			"15" = "NKT",
			"16" = "P-Hepato",
			"17" = "ppLSEC",
			"18" = "LAM-like",
			"19" = "Prolif",
			"20" = "MatB/CD4+",
			"21" = "cDC",
			"22" = "Kupffer",
			"23" = "Stellate",
			"24" = "cNK",
			"25" = "AntiB",
			"26" = "Monocyte-derived",
			"27" = "Treg",
			"28" = "RBC",
			"29" = "Kupffer",
			"30" = "MatB",
			"31" = "Fibroblast",
			"32" = "Cholangiocyte",
			"33" = "pDC",
			"34" = "MatB",
			"35" = "MatB")

require(RColorBrewer)
blues <- brewer.pal(5, "Blues")
greens <- brewer.pal(5, "Greens")

source("~/scripts/LiverMap2.0/Colour_Scheme.R")
# Get same colours as used for subclustering
cell_type_cols <- list(
	"CD8+T"=Cell_type_colours[Cell_type_colours[,1] == "gdTcells2",2],
	"CD4+T"=Cell_type_colours[Cell_type_colours[,1] == "Tcell",2],
	"CD4+T"=Cell_type_colours[Cell_type_colours[,1] == "Tcell",2],
	"NKT"=greens[3],
	"Treg"=greens[5],
	"MatB"=Cell_type_colours[Cell_type_colours[,1] == "Bcell",2],
	"MatB/CD4+"=Cell_type_colours[Cell_type_colours[,1] == "Bcell",2],
	"AntiB"=Cell_type_colours[Cell_type_colours[,1] == "AntiBcell",2],
	"Cellcycle"=Cell_type_colours[Cell_type_colours[,1] == "gdTcells1",2],
	"gdT"=Cell_type_colours[Cell_type_colours[,1] == "gdTcells1",2],
	"lrNK"="#DB8E00",
	"cNK"=Cell_type_colours[Cell_type_colours[,1] == "NKcells",2],
	"miscNK"=Cell_type_colours[Cell_type_colours[,1] == "NKcells",2],
	"Kupffer"=Cell_type_colours[Cell_type_colours[,1] == "NonInfMac",2],
	"Monocyte-derived"=Cell_type_colours[Cell_type_colours[,1] == "InfMac",2],
	"cDC2"=blues[2],
	"cDC"=blues[2],
	"pDC"=blues[3],
	"LAM-like"=blues[5],
	"Neutrophil"="navy",
	"RBC"=Cell_type_colours[Cell_type_colours[,1] == "Eryth",2],
	"Prolif"="black",
	"Hepato"=Cell_type_colours[Cell_type_colours[,1] == "Hepatocyte",2],
	"C-Hepato"=Cell_type_colours[Cell_type_colours[,1] == "Hepatocyte",2],
	"P-Hepato"=Cell_type_colours[Cell_type_colours[,1] == "Eryth",2],
	"Stellate"=Cell_type_colours[Cell_type_colours[,1] == "Stellate",2],
	"Fibroblast"="grey50",
	"Cholangiocyte"=Cell_type_colours[Cell_type_colours[,1] == "Cholangiocyte",2],
	"cvLSEC"=Cell_type_colours[Cell_type_colours[,1] == "cvLSECs",2],
	"ppLSEC"=Cell_type_colours[Cell_type_colours[,1] == "PortalLSECs",2],
	"Arterial"=Cell_type_colours[Cell_type_colours[,1] == "Portalendo",2])
	
convert_using_list <- function(vector_to_convert, conversion_list) {
	coverted <- sapply(vector_to_convert, function(x){conversion_list[[x]]})
	return(coverted)
}
obj@meta.data$cell_type <- convert_using_list(as.character(obj@meta.data$seurat_clusters), cluster2anno)

require(ggplot2)
new_colour_scheme <- unlist(cell_type_cols)[order(names(cell_type_cols))]
names(new_colour_scheme) <- sort(names(cell_type_cols))
new_colour_scheme <- new_colour_scheme[names(new_colour_scheme) %in% obj@meta.data$cell_type]
png(paste(outname, "manual_anno_umap.png", sep="_"), width=8*1, height=6*0.8, units="in", res=300)
DimPlot(obj, reduction="umap", group.by="cell_type", pt.size=.1, label=FALSE) + 
	scale_color_manual(values=new_colour_scheme)
dev.off()

saveRDS(obj, paste(outname, "_Annotated.rds", sep=""))

## Post-Annotation Plots ##
require(ggplot2)
png(paste("Suppl", outname, "neutrophil_markers.png", sep="_"), width=6, height=6, units="in", res=300)
DotPlot(obj, features=c("FCGR3B","ITGAX", "IFITM2", "CXCR1", "CXCR2", "CSF3R"), group.by="cell_type")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

png(paste("Suppl", outname, "NaiveT_markers.png", sep="_"), width=6, height=6, units="in", res=300)
DotPlot(obj, features=c(NaiveTup, NaiveTdn), group.by="cell_type") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

png(paste("Suppl", outname, "DC_markers.png", sep="_"), width=8*4/3, height=6, units="in", res=300)
DotPlot(obj, features=c(general, cDC, pDC), group.by="cell_type")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
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
DotPlot(obj, features=edge_genes, group.by="cell_type")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

Flow_genes <- c("MRC1", "PTPRC", "HLA-DRA", "LYZ", "CD68")
png(paste(outname, "Flowgenes_dotplot.png", sep="_"), width=8, height=8, units="in", res=150)
DotPlot(obj, features=Flow_genes, group.by="cell_type")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# Immune cell population markers
genes=c("EOMES", "CMC1", "KLRG1", "KLRB1", "GZMK", "GZMB", "XCL2", "XCL1", "CD8A", "CD8B", "CD3D", "CD3E", "TRAC", "TRDC", "TRGC2", "IL7R", "IGHM", "IGHD", "IGHA1", "IGKC", "IL2RA", "CTLA4", "FOXP3", "TIGIT", "GATA3", "IL32")
#reordering = c("AntiB", "MatB", "MatB/CD4+", "CD4+T", "Stellate", "ppLSEC", "cvLSEC","Fibroblast", "Kuppfer", "LAM-like", "Monocyte", "Neutrophil", "CD8+T", "gdT", "NKT", "lrNK", "cNK", "Cholangiocyte", )
reordering = c("AntiB", "MatB", "MatB/CD4+", "CD4+T", "Treg", "gdT", "CD8+T", "NKT", "lrNK", "cNK", "Stellate", "ppLSEC", "cvLSEC","Fibroblast", "Cholangiocyte", "Kupffer", "LAM-like", "Monocyte-derived", "Neutrophil", "cDC", "pDC", "Prolif", "RBC", "C-Hepato", "P-Hepato")
obj@meta.data$cell_type <- factor(obj@meta.data$cell_type, levels=reordering)
png(paste(outname, "Immune_markers.png", sep="_"), width=8, height=6, units="in", res=300)
DotPlot(obj[,obj@meta.data$cell_type %in% c("AntiB", "MatB", "Treg", "CD8+T", "gdT", "NKT", "lrNK", "cNK", "MatB/CD4+", "CD4+T")], features=c(genes), group.by="cell_type")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
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
print(DimPlot(obj, reduction="umap", group.by="cell_type", pt.size=.1, label=FALSE) +
        scale_color_manual(values=new_colour_scheme))
dev.off()

png(paste(outname, "TLM_labelled_umap.png", sep="_"), width=8*1, height=6*0.8, units="in", res=300)
print(DimPlot(obj, reduction="umap", group.by="cell_type", pt.size=.1, label=TRUE) +
        scale_color_manual(values=new_colour_scheme))

dev.off()

png(paste(outname, "TLM_IFN.png", sep="_"), width=8*1, height=6*0.8, units="in", res=300)
print(VlnPlot(obj, "IFNG", group.by="cell_type") +
        scale_color_manual(values=new_colour_scheme))
dev.off()


### ------------ ###


#### ----- Subcluster ----- ####

subcluster <- function(obj, clusters=c(0), cluster_col="seurat_clusters") {
	set.seed(101)
	subsample_cells <- obj[,obj@meta.data[,cluster_col] %in% clusters]
	subsample_cells <- RunPCA(subsample_cells)
	subsample_cells <- FindVariableFeatures(subsample_cells, method="vst", nfeatures=1000)
	subsample_cells <- RunPCA(subsample_cells, features=VariableFeatures(subsample_cells), ndims=1:20)
	subsample_cells <- FindNeighbors(subsample_cells, dims=1:12)
	subsample_cells <- FindClusters(subsample_cells, resolution=0.8)
	subsample_cells <- RunUMAP(subsample_cells, dims=1:20)
	png(paste("SubSubcluster_", outname, "_", paste(clusters, collapse="_"), "_umap.png", sep=""), width=8, height=8, units="in", res=300)
	print(DimPlot(subsample_cells))
	dev.off();
	return(subsample_cells)
}
obj <- readRDS(paste(outname, "_Annotated.rds", sep=""))

### For TLM ###
subcluster_prolif <- subcluster(obj, clusters="Prolif", cluster_col="cell_type")

### ---------- ###


subcluster_prolif <- subcluster(obj, clusters=19) #CD8+T & gdT
subcluster_macrophage <- subcluster(obj, clusters=c(8,9,18,21,22,26,29,33))
subcluster_NKT <- subcluster(obj, clusters=c(1,2,3,4,6,7,13,14,15,20,24,27))
subcluster_LSEC <- subcluster(obj, clusters=c(11,17))
subcluster_AntiB <- subcluster(obj, clusters=c(25))

subclustering_list <- list("Prolif_19-Subcluster"=subcluster_prolif,
			"NKT-Subcluster"=subcluster_NKT,
			"LSEC-Subcluster"=subcluster_LSEC,
			"Mac-Subcluster"=subcluster_macrophage,
			"AntiB-Subcluster"=subcluster_AntiB)

# Subcluster Anno

Global_Marker_Table<- read.delim("Map2.0_Global_Marker_Table.csv", sep=",")
all_groups <- unique(c(Global_Marker_Table[,2]))

for( i in 1:length(subclustering_list) ) {
	tag = names(subclustering_list)[i]
	subcluster_obj <- subclustering_list[[i]]

	sub_cluster_de=list()
	all_groups <- unique(c(Global_Marker_Table[,2]))
	all_scores <- matrix(0, nrow=length(all_groups), ncol=length(unique(subcluster_obj@meta.data$seurat_clusters)))
	rownames(all_scores) <- all_groups;
	colnames(all_scores) <- sort(unique(subcluster_obj@meta.data$seurat_clusters));

	for(c in subcluster_obj@meta.data$seurat_clusters) {
	        DE <- FindMarkers(subcluster_obj, ident.1=c, group.by="seurat_clusters", test.use="wilcox")
		sub_cluster_de[[c]]<-DE
	        anno <- auto_anno_map2.0(sub_cluster_de[[c]], paste("Cluster", c))
	        anno <- anno[match(rownames(all_scores), anno[,1]),2]
	        all_scores[,which(colnames(all_scores)==c)] <- anno
	}

	saveRDS(list(obj=subcluster_obj, 
			de=sub_cluster_de, 
			anno_scores=all_scores), 
		      file=paste(tag, ".rds", sep=""))
}

#### ----- END ----- ####

#### ----- Compare to Healthy ----- ####

# Load 5' Healthy map.

HealthyMap_RDS <- "/cluster/projects/macparland/TA/LiverMap2.0/RawData/Analysis/EmptyDrops/Merged_EmptyOnly_obj_Map2.2_ImportedClusters_ManualAnno_dimreduce.rds"

require("Seurat")
source("~/scripts/LiverMap2.0/My_R_Scripts.R")

mergedobj <- readRDS(HealthyMap_RDS)
healthy_5pr <- mergedobj[,mergedobj@meta.data$assay_type == "5pr"]
healthy_5pr_save <- mergedobj[,mergedobj@meta.data$assay_type == "5pr"]
rm(mergedobj)

obj <- readRDS(paste(outname, "_Annotated.rds", sep=""))
healthy_5pr@meta.data$full_annotation <- as.character(healthy_5pr@meta.data$Coarse_Manual_Anno)
for (f in Sys.glob("/cluster/projects/macparland/TA/LiverMap2.0/RawData/Analysis/EmptyDrops/Subcluster/FinalMetadata/*_fullmetadata.rds")) {
	subcluster_meta <- readRDS(f)
	matched_anno <- as.character(subcluster_meta$Subcluster_Manual)[match(healthy_5pr@meta.data$cell_ID, subcluster_meta$cell_ID)]
	matched_anno[is.na(matched_anno)] <- healthy_5pr@meta.data$full_annotation[is.na(matched_anno)]
	healthy_5pr@meta.data$full_annotation <- matched_anno
}

# Manual Fixups
healthy_5pr@meta.data$full_annotation[healthy_5pr@meta.data$Coarse_Manual_Anno == "AntiBcell"] <- "AntiB"
healthy_5pr@meta.data$full_annotation[healthy_5pr@meta.data$Coarse_Manual_Anno == "Stellate"] <- "Stellate"
healthy_5pr@meta.data$full_annotation[healthy_5pr@meta.data$Coarse_Manual_Anno == "Cholangiocyte"] <- "Cholangiocyte"
healthy_5pr@meta.data$full_annotation[healthy_5pr@meta.data$Coarse_Manual_Anno == "MHCII"] <- "MHCII/DC"
healthy_5pr@meta.data$full_annotation[healthy_5pr@meta.data$full_annotation == "MatBcells"] <- "MatB"
healthy_5pr@meta.data$full_annotation[healthy_5pr@meta.data$full_annotation == "Resident-Kupffer"] <- "Kupffer"
obj@meta.data$cell_type[obj@meta.data$cell_type == "Monocyte-derived"] <- "Monocyte"
obj@meta.data$cell_type[obj@meta.data$cell_type == "cDC2"] <- "MHCII/DC"
obj@meta.data$cell_type[obj@meta.data$cell_type == "cDC"] <- "MHCII/DC"
obj@meta.data$cell_type[obj@meta.data$cell_type == "Hepato"] <- "C-Hepato"
healthy_5pr@meta.data$full_annotation[healthy_5pr@meta.data$full_annotation == "Central"] <- "C-Hepato"
#healthy_5pr@meta.data$full_annotation[healthy_5pr@meta.data$full_annotation == "PortalHep"] <- "P-Hepato"
healthy_5pr@meta.data$full_annotation[healthy_5pr@meta.data$full_annotation == "Portal"] <- "P-Hepato" # fiexed 14 Jun 2022
healthy_5pr@meta.data$full_annotation[healthy_5pr@meta.data$full_annotation == "Prolif"] <- "Cellcycle"
healthy_5pr@meta.data$full_annotation <- sub("cell", "", healthy_5pr@meta.data$full_annotation)

# Synced Maps
# Get same colours as used for subclustering

tmp <- strsplit(obj@meta.data$orig.ident, "_")
obj@meta.data$donor <- unlist(lapply(tmp, function(x){x[[1]]}))

# Save output for cellphonedb
saveRDS(list(healthy_5pr=healthy_5pr, psc=obj), file=paste(outname, "Objects_forCellCellInterations.rds", sep="_"))

x <- readRDS(paste(outname, "Objects_forCellCellInterations.rds", sep="_"))
source("~/scripts/LiverMap2.0/My_R_Scripts.R")
healthy_5pr <- x$healthy_5pr
obj <- x$psc

# Random Figures #

png(paste(outname, "CD28_psc.png", sep="_"), width=6, height=6, unit="in", res=100)
FeaturePlot(x$psc, "CD28")
dev.off()
png(paste(outname, "CD28_healthy.png", sep="_"), width=6, height=6, unit="in", res=100)
FeaturePlot(x$healthy, "CD28")
dev.off()

## Custom_thing - TLM ##
# Manual survey of cytokine interactions (from KEGG):
# CCL21, CCL19 -> CCR7, ACKR4
# CXCL10 -> CXCR3
# IFNG -> IFNGR1, IFNGR2
# IL2 -> IL2RA, IL2RB, IL2RG, IL2RB, IL2RG
# IL21R -> IL21
# CXCR6 -> CXCL16
# CCL3 -> CCR5, CCR1
# CCL4 -> CCR5
# CCL5 -> CCR3, CCR5, CCR4, CCR1
# CX3CR1 -> CX3CL1
# IL32 ->
# IL6 -> IL6R, IL6ST
# CXCL8 -> CXCR1, CXCR2

# Pathway IL2 -> NK/T -> IL32 -> Macs -> TNFalpha
genes <- c("IFNG", "IFNA1", "TNF", "CXCL8", "CXCL10","CCL21", "CCL19", "CCR7", "ACKR4", )
healthy_mean <- group_rowmeans(healthy_5pr@assays$RNA@data, healthy_5pr@meta.data$full_annotation)
healthy_detect <- group_rowmeans(healthy_5pr@assays$RNA@counts > 0, healthy_5pr@meta.data$full_annotation)
psc_mean <- group_rowmeans(obj@assays$RNA@data, obj@meta.data$cell_type)
psc_detect <- group_rowmeans(obj@assays$RNA@counts > 0, obj@meta.data$cell_type)

Cytokines <- c("CCL21", "CCL19", "CCR7", "ACKR4", "CXCL10", "CXCR3", "IFNG", "IGNGR1", "IFNGR2",
		"IL2", "IL2RA", "IL2RB", "IL2RG", "IL21", "IL21R", 
		"CXCL16", "CXCR6", "CCL3", "CCL4", "CCL5", "CCR5", "CCR1", "CCR4", "CCR3",
		"CX3CL1", "CX3CR1", "IL6", "IL6R", "IL6ST", "CXCL8", "CXCR1", "CXCR2")

png(paste(outname, "TLM_cytokines_psc.png", sep="_"), width=8, height=6,units="in", res=150)
DotPlot(obj, features=Cytokines, group.by="cell_type")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
png(paste(outname, "TLM_cytokines_healthy.png", sep="_"), width=8, height=6,units="in", res=150)
DotPlot(healthy_5pr, features=Cytokines, group.by="full_annotation")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

psc_detect[c("CXCL8", "CXCR1", "CXCR2"), c("Cholangiocyte", "Stellate", "ppLSEC", "cNK", "lrNK", "CD8+T", "CD4+T", "AntiB")]
healthy_detect[c("CXCL8", "CXCR1", "CXCR2"), c("Cholangiocyte", "Stellate", "ppLSEC", "cNK", "lrNK", "CD8+T", "CD4+T", "AntiB")]
these_colors <- c(cell_type_cols[["Cholangiocyte"]], cell_type_cols[["Stellate"]], 
cell_type_cols[["ppLSEC"]], cell_type_cols[["cNK"]], 
cell_type_cols[["lrNK"]], cell_type_cols[["CD8+T"]], 
cell_type_cols[["CD4+T"]], cell_type_cols[["AntiB"]])

these_types <- c("Cholangiocyte", "Stellate", "ppLSEC", "cNK", "lrNK", "CD8+T", "CD4+T", "AntiB")

tab1 <- rbind(psc_detect["CXCL8", these_types], 
		healthy_detect["CXCL8", these_types])

tab2 <- rbind(colSums(psc_detect[c("CXCR1", "CXCR2"),these_types]), 
		colSums(healthy_detect[c("CXCR1", "CXCR2"), these_types]))

png(paste(outname, "TLM_CXCL8.png", sep="_"), width=8*0.5, height=6*2*0.5,units="in", res=150)
par(mfrow=c(2,1))
barplot(tab1[1,], names=these_types, col=these_colors, ylab="Detection Rate", main="CXCL8")
barplot(tab2[1,], names=these_types, col=these_colors, ylab="Detection Rate", main="CXCR1/CXCR2")
dev.off()

NK_enrichments <- data.frame(term=c("NK-mediated cytotoxicity", "Antigen presentation", "Th1/T2 differentiation", "Chemokine signaling", "Endocytosis"), pvalue=c(1.947*10^-16, 1.072*10^-14, 2.334*10^-6, 2.536*10^-5, 5.798*10^-4)) # from gprofiler


tab3 <- rbind((psc_detect[c("IFNG"),these_types]), 
		(healthy_detect[c("IFNG"), these_types]))
# ========== #

## Cell-Type Frequncie - in silico Gating ##

#obj <- readRDS(paste(outname, "_Annotated.rds", sep=""))
x <- readRDS(paste(outname, "Objects_forCellCellInterations.rds", sep="_"))
healthy_5pr <- x$healthy_5pr
obj <- x$psc
rm(x)
source("~/scripts/LiverMap2.0/My_R_Scripts.R")
healthy_mean <- group_rowmeans(healthy_5pr@assays$RNA@data, factor(healthy_5pr@meta.data$full_annotation))
healthy_detect <- group_rowmeans(healthy_5pr@assays$RNA@counts > 0, factor(healthy_5pr@meta.data$full_annotation))
psc_mean <- group_rowmeans(obj@assays$RNA@data, factor(obj@meta.data$cell_type))
psc_detect <- group_rowmeans(obj@assays$RNA@counts > 0, factor(obj@meta.data$cell_type))

healthy_freqs <- table(factor(healthy_5pr@meta.data$full_annotation))
psc_freqs <- table(factor(obj@meta.data$cell_type))
healthy_freqs_by_sample <- table(factor(paste(healthy_5pr@meta.data$full_annotation, healthy_5pr@meta.data$sample)))
psc_freqs_by_sample <-  table(factor(paste(obj@meta.data$cell_type, obj@meta.data$sample)))

my_split_names <- function(x) {matrix(unlist(strsplit(names(x), " ")), nrow=2, ncol=length(x))}

## What % of PTPRC+ are MRC1+? ##
#png(paste(outname, bg_marker, "healthy_umap.png", sep="_"), width=6, height=4, units="in", res=150)
#print(FeaturePlot(healthy_5pr, bg_marker))
#dev.off()
#png(paste(outname, bg_marker, "psc_umap.png", sep="_"), width=6, height=4, units="in", res=150)
#print(FeaturePlot(obj, bg_marker))
#dev.off()
#png(paste(outname, fg_marker, "healthy_umap.png", sep="_"), width=6, height=4, units="in", res=150)
#print(FeaturePlot(healthy_5pr, fg_marker))
#dev.off()
#png(paste(outname, fg_marker, "psc_umap.png", sep="_"), width=6, height=4, units="in", res=150)
#print(FeaturePlot(obj, fg_marker))
#dev.off()

# Set up functions

get_sample_counts <- function(these_types, freqs_by_sample) {
        sample_labels <-  my_split_names(freqs_by_sample)
        tmp <- freqs_by_sample[sample_labels[1,] %in% these_types]
        lab <- sample_labels[2, sample_labels[1,] %in% these_types]
        d <- split(seq(length(tmp)), lab)
        return(sapply(d, function(X){sum(tmp[X])}))
}

# Healthy freqs
bg_types <- colnames(healthy_detect)[healthy_detect["PTPRC",] > 0.1]
fg_types <- colnames(healthy_detect)[healthy_detect["MRC1",] > 0.1]
fg_types <- fg_types[fg_types %in% bg_types]

n_fg_healthy <- sum(healthy_freqs[names(healthy_freqs) %in% fg_types])
n_bg_healthy <- sum(healthy_freqs[names(healthy_freqs) %in% bg_types])

n_fg_healthy_sample <- get_sample_counts(fg_types, healthy_freqs_by_sample)
n_bg_healthy_sample <- get_sample_counts(bg_types, healthy_freqs_by_sample)


# PSC freqs
bg_types <- colnames(psc_detect)[psc_detect["PTPRC",] > 0.1]
fg_types <- colnames(psc_detect)[psc_detect["MRC1",] > 0.2]
fg_types <- fg_types[fg_types %in% bg_types]

n_fg_psc <- sum(psc_freqs[names(psc_freqs) %in% fg_types])
n_bg_psc <- sum(psc_freqs[names(psc_freqs) %in% bg_types])

n_fg_psc_sample <- get_sample_counts(fg_types, psc_freqs_by_sample)
n_bg_psc_sample <- get_sample_counts(bg_types, psc_freqs_by_sample)

t.test( n_fg_psc_sample/n_bg_psc_sample, n_fg_healthy_sample/n_bg_healthy_sample)

pval <- fisher.test(rbind(c(n_fg_healthy, n_bg_healthy-n_fg_healthy), c(n_fg_psc, n_bg_psc-n_fg_psc)))
if (pval$p.value < 10^-10) {
	n_stars <- floor(log10(-1*log10(pval$p.value)))+1
	stars <- paste(rep("*", n_stars), collapse="")
} else if (pval$p.value < 0.05) {
	stars <- "*"
} else {
	stars <- ""
}
png(paste(outname, "MRC1_vs_PTPRC_barplot.png", sep="_"), width=2, height=4, units="in", res=300)
par(mar=c(8,4,1,1))
ymax = max(c(n_fg_psc/n_bg_psc,  n_fg_healthy/n_bg_healthy))
coords<-barplot(c(n_fg_psc/n_bg_psc, n_fg_healthy/n_bg_healthy), names=c("PSC(5'sc)", "Healthy(5'sc)"), col=c("black", "grey65"), ylab="Percent of MRC1+ Cells (%)", las=2, ylim= c(0, ymax*1.5))
lines(coords, y=rep(ymax*1.1, 2) )
text(mean(coords), ymax*1.1, pos=3, labels=stars)
dev.off()


## CD68+ PTPRC+ vs PTPRC+ ##
bg_types <-  colnames(healthy_detect)[healthy_detect["PTPRC",] > 0.1]
fg_types <- colnames(healthy_detect)[healthy_detect["CD68",] > 0.1]
fg_types <- fg_types[fg_types %in% bg_types]

n_fg_healthy <- sum(healthy_freqs[names(healthy_freqs) %in% fg_types])
n_bg_healthy <- sum(healthy_freqs[names(healthy_freqs) %in% bg_types])

n_fg_healthy_sample <- get_sample_counts(fg_types, healthy_freqs_by_sample)
n_bg_healthy_sample <- get_sample_counts(bg_types, healthy_freqs_by_sample)

bg_types <-  colnames(psc_detect)[psc_detect["PTPRC",] > 0.1]
fg_types <- colnames(psc_detect)[psc_detect["CD68",] > 0.6]
fg_types <- fg_types[fg_types %in% bg_types]

n_fg_psc <- sum(psc_freqs[names(psc_freqs) %in% fg_types])
n_bg_psc <- sum(psc_freqs[names(psc_freqs) %in% bg_types])

n_fg_psc_sample <- get_sample_counts(fg_types, psc_freqs_by_sample)
n_bg_psc_sample <- get_sample_counts(bg_types, psc_freqs_by_sample)

pval <- fisher.test(rbind(c(n_fg_healthy, n_bg_healthy-n_fg_healthy), c(n_fg_psc, n_bg_psc-n_fg_psc)))

pval <- t.test( n_fg_psc_sample/n_bg_psc_sample, n_fg_healthy_sample/n_bg_healthy_sample)

if (pval$p.value < 10^-10) {
        n_stars <- floor(log10(-1*log10(pval$p.value)))+1
        stars <- paste(rep("*", n_stars), collapse="")
} else if (pval$p.value < 0.05) {
        stars <- "*"
} else {
        stars <- ""
}
png(paste(outname, "CD68_vs_PTPRC_barplot.png", sep="_"), width=2, height=4, units="in", res=300)
par(mar=c(8,4,1,1))
ymax = max(c(n_fg_psc/n_bg_psc,  n_fg_healthy/n_bg_healthy))
coords<-barplot(c(n_fg_psc/n_bg_psc, n_fg_healthy/n_bg_healthy), names=c("PSC(5'sc)", "Healthy(5'sc)"), col=c("black", "grey65"), ylab="Percent of CD68+ Cells (%)", las=2, ylim= c(0, 0.8))
points(rep(coords[1,1], length(n_fg_psc_sample)), n_fg_psc_sample/n_bg_psc_sample, pch=16, col="black", cex=0.75)
points(rep(coords[2,1], length(n_fg_healthy_sample)), n_fg_healthy_sample/n_bg_healthy_sample, pch=16, col="black", cex=0.75)
lines(coords, y=rep(ymax*1.1, 2) )
text(mean(coords), ymax*1.1, pos=3, labels=stars)
dev.off()

## MCR1 vs CD68
bg_types <-  colnames(healthy_detect)[healthy_detect["PTPRC",] > 0.1 & healthy_detect["CD68",] > 0.1]
fg_types <- colnames(healthy_detect)[healthy_detect["MRC1",] > 0.1]
fg_types <- fg_types[fg_types %in% bg_types]

n_fg_healthy <- sum(healthy_freqs[names(healthy_freqs) %in% fg_types])
n_bg_healthy <- sum(healthy_freqs[names(healthy_freqs) %in% bg_types])

bg_types <-  colnames(psc_detect)[psc_detect["PTPRC",] > 0.1 & psc_detect["CD68",] > 0.6]
fg_types <- colnames(psc_detect)[psc_detect["MRC1",] > 0.2]
fg_types <- fg_types[fg_types %in% bg_types]

n_fg_psc <- sum(psc_freqs[names(psc_freqs) %in% fg_types])
n_bg_psc <- sum(psc_freqs[names(psc_freqs) %in% bg_types])

pval <- fisher.test(rbind(c(n_fg_healthy, n_bg_healthy-n_fg_healthy), c(n_fg_psc, n_bg_psc-n_fg_psc)))
if (pval$p.value < 10^-10) {
        n_stars <- floor(log10(-1*log10(pval$p.value)))+1
        stars <- paste(rep("*", n_stars), collapse="")
} else if (pval$p.value < 0.05) {
        stars <- "*"
} else {
        stars <- ""
}
png(paste(outname, "MRC1_vs_CD68_barplot.png", sep="_"), width=2, height=4, units="in", res=300)
par(mar=c(8,4,1,1))
ymax = max(c(n_fg_psc/n_bg_psc,  n_fg_healthy/n_bg_healthy))
coords<-barplot(c(n_fg_psc/n_bg_psc, n_fg_healthy/n_bg_healthy), names=c("PSC(5'sc)", "Healthy(5'sc)"), col=c("black", "grey65"), ylab="Percent of MRC1+ Cells (%)", las=2, ylim= c(0, ymax*1.5))
lines(coords, y=rep(ymax*1.1, 2) )
text(mean(coords), ymax*1.1, pos=3, labels=stars)
dev.off()

## LYZ vs CD68
bg_types <-  colnames(healthy_detect)[healthy_detect["PTPRC",] > is.on.dthreshold & healthy_detect["CD68",] > is.on.dthreshold]
fg_types <- colnames(healthy_mean)[healthy_mean["LYZ",] > 1]
fg_types <- fg_types[fg_types %in% bg_types]

n_fg_healthy <- sum(healthy_freqs[names(healthy_freqs) %in% fg_types])
n_bg_healthy <- sum(healthy_freqs[names(healthy_freqs) %in% bg_types])

bg_types <-  colnames(psc_detect)[psc_detect["PTPRC",] > is.on.dthreshold & psc_detect["CD68",] > 0.6]
fg_types <- colnames(psc_mean)[psc_mean["LYZ",] > 2]
fg_types <- fg_types[fg_types %in% bg_types]

n_fg_psc <- sum(psc_freqs[names(psc_freqs) %in% fg_types])
n_bg_psc <- sum(psc_freqs[names(psc_freqs) %in% bg_types])

pval <- fisher.test(rbind(c(n_fg_healthy, n_bg_healthy-n_fg_healthy), c(n_fg_psc, n_bg_psc-n_fg_psc)))
if (pval$p.value < 10^-10) {
        n_stars <- floor(log10(-1*log10(pval$p.value)))+1
        stars <- paste(rep("*", n_stars), collapse="")
} else if (pval$p.value < 0.05) {
        stars <- "*"
} else {
        stars <- ""
}
png(paste(outname, "LYZ_vs_CD68_barplot.png", sep="_"), width=2, height=4, units="in", res=300)
par(mar=c(8,4,1,1))
ymax = max(c(n_fg_psc/n_bg_psc,  n_fg_healthy/n_bg_healthy))
coords<-barplot(c(n_fg_psc/n_bg_psc, n_fg_healthy/n_bg_healthy), names=c("PSC(5'sc)", "Healthy(5'sc)"), col=c("black", "grey65"), ylab="Percent of LYZ+ Cells (%)", las=2, ylim= c(0, ymax*1.5))
lines(coords, y=rep(ymax*1.1, 2) )
text(mean(coords), ymax*1.1, pos=3, labels=stars)
dev.off()


## HLADR vs CD68+/CD14+ PTPRC+ CD3-
bg_types <- colnames(healthy_detect)[
		healthy_detect["PTPRC",] > is.on.dthreshold & 
		(healthy_detect["CD68",] > is.on.dthreshold | healthy_detect["CD14",] > 0.51) &
		(healthy_detect["CD3E",] < is.on.dthreshold & healthy_detect["CD3D",] < is.on.dthreshold)
		]
fg_types <- colnames(healthy_detect)[
                healthy_mean["HLA-DRA",] > 2 |  healthy_mean["HLA-DPA1",] > 2 |  healthy_mean["HLA-DQB1",] > 2 |  healthy_mean["HLA-DPB1",] > 2
                ]
fg_types <- fg_types[fg_types %in% bg_types]

n_fg_healthy <- sum(healthy_freqs[names(healthy_freqs) %in% fg_types])
n_bg_healthy <- sum(healthy_freqs[names(healthy_freqs) %in% bg_types])

bg_types <- colnames(psc_detect)[
                psc_detect["PTPRC",] > is.on.dthreshold &
                (psc_detect["CD68",] > 0.6 | psc_detect["CD14",] > 0.51) &
                (psc_detect["CD3E",] < is.on.dthreshold & psc_detect["CD3D",] < is.on.dthreshold)
                ]
fg_types <- colnames(psc_detect)[
                psc_mean["HLA-DRA",] > 2 |  psc_mean["HLA-DPA1",] > 2 |  psc_mean["HLA-DQB1",] > 2 |  psc_mean["HLA-DPB1",] > 2
                ]
fg_types <- fg_types[fg_types %in% bg_types]

n_fg_psc <- sum(psc_freqs[names(psc_freqs) %in% fg_types])
n_bg_psc <- sum(psc_freqs[names(psc_freqs) %in% bg_types])

pval <- fisher.test(rbind(c(n_fg_healthy, n_bg_healthy-n_fg_healthy), c(n_fg_psc, n_bg_psc-n_fg_psc)))

if (pval$p.value < 10^-10) {
        n_stars <- floor(log10(-1*log10(pval$p.value)))+1
        stars <- paste(rep("*", n_stars), collapse="")
} else if (pval$p.value < 0.05) {
        stars <- "*"
} else {
        stars <- ""
}
png(paste(outname, "Macrophage_HLADR_Prop_barplot.png", sep="_"), width=2, height=4, units="in", res=300)
par(mar=c(8,4,1,1))
ymax = max(c(n_fg_psc/n_bg_psc,  n_fg_healthy/n_bg_healthy))
coords<-barplot(c(n_fg_psc/n_bg_psc, n_fg_healthy/n_bg_healthy), names=c("PSC(5'sc)", "Healthy(5'sc)"), col=c("black", "grey65"), ylab="Percent of HLADR Meyloid cells (%)", las=2, ylim= c(0, ymax*1.5))
lines(coords, y=rep(ymax*1.1, 2) )
text(mean(coords), ymax*1.1, pos=3, labels=stars)
dev.off()




###########################

tab <- obj@assays$RNA@counts
tab <- tab[unique(rownames(tab)),]
met <- data.frame(Cell=obj@meta.data[,"cell_ID"], cell_type=obj@meta.data[,"cell_type"])
colnames(tab) <- as.character(met[,1])
write.table(tab, paste(outname, "psc_counts.txt", sep="_"), col.names=T, row.names=T)
write.table(met, paste(outname, "psc_meta.txt", sep="_"), col.names=T, row.names=F, quote=F)

tab <- healthy_5pr@assays$RNA@counts
tab <- tab[unique(rownames(tab)),]
met <- data.frame(Cell=healthy_5pr@meta.data[,"cell_ID"], cell_type=healthy_5pr@meta.data[,"full_annotation"])
colnames(tab) <- as.character(met[,1])
write.table(tab, paste(outname, "norm_counts.txt", sep="_"), col.names=T, row.names=T)
write.table(met, paste(outname, "norm_meta.txt", sep="_"), col.names=T)

# Proliferation
healthy_5pr <- Seurat::CellCycleScoring(healthy_5pr, g2m.features=cc.genes$g2m.genes, s.features=cc.genes$s.genes)
obj <- Seurat::CellCycleScoring(obj, g2m.features=cc.genes$g2m.genes, s.features=cc.genes$s.genes)
CC_healthy <- table(healthy_5pr@meta.data$full_annotation, healthy_5pr@meta.data$Phase)
CC_psc <- table(obj@meta.data$cell_type, obj@meta.data$Phase)
CC_healthy <- CC_healthy/rowSums(CC_healthy)
CC_psc <- CC_psc/rowSums(CC_psc)


common <- rownames(CC_psc)[rownames(CC_psc) %in% rownames(CC_healthy)]
Prop_G2M <- data.frame(PSC=CC_psc[match(common, rownames(CC_psc)),2],
                        Norm=CC_healthy[match(common, rownames(CC_healthy)),2])

Prop_S <- data.frame(PSC=CC_psc[match(common, rownames(CC_psc)),3],
                        Norm=CC_healthy[match(common, rownames(CC_healthy)),3])

Prop_G2M+Prop_S
mat <- as.matrix(Prop_G2M+Prop_S); mat <- mat[c(1,13,8,7,9,4,5,11,12,3,2,10,6),]

png(paste(outname, "Cycling_barplot.png", sep="_"), width=2*1.5, height=4*1.5, units="in", res=150)
par(mar=c(3,6.5,1,1))
barplot(t(as.matrix(mat))*100, beside=T, horiz=T, col=c("grey35", "grey65"), las=2, xlab="Cycling (%)")
dev.off()

# Sync datasets
healthy_5pr <- healthy_5pr[,healthy_5pr@meta.data$full_annotation %in% unique(obj@meta.data$cell_type)]

HEALTHY_5pr_pseudo <- get_pseudobulk(healthy_5pr@assays$RNA@counts, healthy_5pr@meta.data$full_annotation, healthy_5pr@meta.data$donor)
PSC_5pr_pseudo <- get_pseudobulk(obj@assays$RNA@counts, obj@meta.data$cell_type, obj@meta.data$donor)

saveRDS(list(healthy_meta=healthy_5pr@meta.data, PSC_meta=obj@meta.data, healthy=HEALTHY_5pr_pseudo, PSC=PSC_5pr_pseudo), paste(outname, "pseudobulks.rds", sep="_"))



require(ggplot2)
new_colour_scheme <- unlist(cell_type_cols)[order(names(cell_type_cols))]
names(new_colour_scheme) <- sort(names(cell_type_cols))
png(paste(outname, "manual_anno_umap_synced.png", sep="_"), width=6*1.5, height=6*1.2, units="in", res=300)
DimPlot(obj, reduction="umap", group.by="cell_type", pt.size=.1, label=TRUE) +
        scale_color_manual(values=new_colour_scheme)
dev.off()
new_colour_scheme <- unlist(cell_type_cols)[order(names(cell_type_cols))]
names(new_colour_scheme) <- sort(names(cell_type_cols))
new_colour_scheme<-new_colour_scheme[names(new_colour_scheme) %in% healthy_5pr@meta.data$full_annotation]

set.seed(1012)
healthy_5pr <- RunUMAP(healthy_5pr, dims=1:20)

png(paste(outname, "manual_anno_umap_syncedhealthy_newumap.png", sep="_"), width=6*1.5, height=6*1.2, units="in", res=300)
DimPlot(healthy_5pr, reduction="umap", group.by="full_annotation", pt.size=.1, label=TRUE) +
        scale_color_manual(values=new_colour_scheme)
dev.off()

#### DE - edgeR ####

require(edgeR)
pseudos <- readRDS(paste(outname, "pseudobulks.rds", sep="_"))

HEALTHY_5pr_pseudo <- pseudos$healthy
PSC_5pr_pseudo <- pseudos$PSC

## Remove Hepato genes ##

#hepatoDE <- read.table("PSC_SC_5pr_3_ClusterDE.csv", sep=",")
#hepatoDE_genes <- rownames(hepatoDE[hepatoDE$avg_logFC > 0.7,])

hepatoDE1 <- read.table("/cluster/projects/macparland/TA/LiverMap2.0/RawData/Analysis/EmptyDrops/DE_CentralHep_Coarse_Manual_wilcox.csv", header=T, sep=",")
hepatoDE2 <- read.table("/cluster/projects/macparland/TA/LiverMap2.0/RawData/Analysis/EmptyDrops/DE_PortalHep_Coarse_Manual_wilcox.csv", header=T, sep=",")
hepatoDE_genes <- unique(c(rownames(hepatoDE1[hepatoDE1$avg_logFC > 0 & hepatoDE1$pct.2 > 1/3,]), 
			   rownames(hepatoDE2[hepatoDE2$avg_logFC > 0 & hepatoDE2$pct.2 > 1/3,]))) # 18/08/2022

## EdgeR ##

common_genes <- sort(intersect(rownames(HEALTHY_5pr_pseudo), rownames(PSC_5pr_pseudo)))
#common_genes <- common_genes[!common_genes %in% hepatoDE_genes]

HEALTHY_5pr_pseudo <- HEALTHY_5pr_pseudo[match(common_genes, rownames(HEALTHY_5pr_pseudo)),]
PSC_5pr_pseudo <- PSC_5pr_pseudo[match(common_genes, rownames(PSC_5pr_pseudo)),]

types <- sort(unique(sapply(strsplit(colnames(HEALTHY_5pr_pseudo), "_"), "[", 1)))
types <- sub("\\+T","", types)
types <- sub("\\/DC","", types)

set.seed(39029)
for(type in types) {
	psc <- round(PSC_5pr_pseudo[,grep(type, colnames(PSC_5pr_pseudo))])
	norm <- HEALTHY_5pr_pseudo[,grep(type, colnames(HEALTHY_5pr_pseudo))]
	if (ncol(norm) == 0 | ncol(psc) == 0) {print(paste("Warning: skipping", type)); next;}
	pseudobulk_matrix <- cbind(psc,norm)
	groups <- c(rep("PSC", ncol(psc)), rep("Norm", ncol(norm)))

	require(edgeR)
	eRobj <- DGEList(counts=pseudobulk_matrix, group=groups)
	eRobj <- calcNormFactors(eRobj)
	eRobj1 <- estimateCommonDisp(eRobj, verbose=T)
	eRobj1 <- estimateTagwiseDisp(eRobj1)
	design.mat <- model.matrix(~ 0 + eRobj1$samples$group)
	de <- exactTest(eRobj1, pair=c("PSC", "Norm"))
	de_mat<- data.frame(topTags(de, nrow(pseudobulk_matrix)))
	de_mat <- de_mat[de_mat$FDR < 0.05,]
	de_mat <- de_mat[! rownames(de_mat) %in% hepatoDE_genes,]
	write.table(de_mat, paste(outname, type, "edgeR_pseudobulkDE_newhepatofilter.csv", sep="_"), sep=",")
}

## Comparison Plots ##

# Relative Frequencies of Cell-types
summary_data <- readRDS(paste(outname, "pseudobulks.rds", sep="_"))
# B-cells vs all lymphocytes
cell_type_ns_healthy <- table(summary_data$healthy_meta$full_annotation)
cell_type_ns_psc <- table(summary_data$PSC_meta$cell_type)

lympho <- c("CD8+T", "CD4+T", "Treg", "cNK", "lrNK", "MatB", "AntiB", "gdT","NKT")
myeloid <- c("LAM-like","MHCII/DC", "Kupffer", "Monocyte", "pDC", "Neutrophil") #exclude Neutrophils from sum?
nonPara <- c("cvLSEC", "ppLSEC", "Stellate","Fibroblast", "Cholangiocyte")

cell_type_cols[["MHCII/DC"]] <- cell_type_cols[["cDC"]]
cell_type_cols[["Monocyte"]] <- cell_type_cols[["Monocyte-derived"]]


#pie(cell_type_ns_healthy[names(cell_type_ns_healthy) %in% lympho])
#pie(cell_type_ns_psc[names(cell_type_ns_psc) %in% lympho])
#pie(cell_type_ns_healthy[names(cell_type_ns_healthy) %in% myeloid])
#pie(cell_type_ns_psc[names(cell_type_ns_psc) %in% myeloid])
#pie(cell_type_ns_healthy[names(cell_type_ns_healthy) %in% nonPara])
#pie(cell_type_ns_psc[names(cell_type_ns_psc) %in% nonPara])

tmp1 <- cell_type_ns_healthy[match(nonPara, names(cell_type_ns_healthy))]
names(tmp1) <- nonPara
tmp1[is.na(tmp1)] <- 0
tmp2 <- cell_type_ns_psc[match(nonPara, names(cell_type_ns_psc))]
png(paste(outname, "Relative_cellnumbers_nonPara.png", sep="_"), width=6*0.7, height=7*0.7, units="in", res=150)
barplot(cbind(tmp1/sum(tmp1), tmp2/sum(tmp2))*100, names=c("Control", "PSC"), col=convert_using_list(names(tmp1), cell_type_cols), legend.text=nonPara, ylab="Percent of Cells (%)", xlim=c(0,6), ylim=c(0, 125), las=2)
dev.off()

#### DE Plots ####

all_disease_de <- list()
files <- Sys.glob(paste(outname,"*edgeR_pseudobulkDE_noHep2.csv", sep=""))
for (f in files) {
	tab <- read.table(f, sep=",")
	tmp <- unlist(strsplit(f, "_"))
	type <- tmp[length(tmp)-3]
	all_disease_de[[type]] <- tab;
}

tmp<- all_disease_de[["P-Hepato"]]
sort(rownames(tmp)[tmp[,1] < 0])

# Example where same gene different directions in different cell-types?
up <- c()
down <- c()
for(t in all_disease_de){
	up <- c(up, rownames(t)[t$logFC >0])
	down <- c(down, rownames(t)[t$logFC <0])
}

sum(up %in% down)

of_interest <- unique(up[up%in%down])
mat <- matrix(0, nrow=length(of_interest), ncol=length(all_disease_de))
rownames(mat) <- of_interest
colnames(mat) <- names(all_disease_de)
for(t in names(all_disease_de)){
	FC <- all_disease_de[[t]]$logFC
	FC <- FC[match(rownames(mat), rownames(all_disease_de[[t]]))]
	FC[is.na(FC)] <- 0
	mat[,colnames(mat)==t] <- FC*-1
}
to_plot <- as.character(t(round(mat, digits=1)))
to_plot[to_plot =="0"] <- "";
png(paste(outname, "Opposite_genes.png", sep="_"), width=6*1.5, height=6, units="in", res=300)
require(RColorBrewer)
cols <- rev(brewer.pal(7,"RdBu"))
image(t(round(mat, digits=1)), col=cols, xaxt=NULL, yaxt=NULL, axes=FALSE)
axis(1, at=seq(from=0, to=1, length=ncol(mat)), labels=colnames(mat), las=2)
axis(2, at=seq(from=0, to=1, length=nrow(mat)), labels=rownames(mat), las=2)
dev.off()


image.scale <- function(z, zlim, col = heat.colors(12),
breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
 if(!missing(breaks)){
  if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
 }
 if(missing(breaks) & !missing(zlim)){
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
 }
 if(missing(breaks) & missing(zlim)){
  zlim <- range(z, na.rm=TRUE)
  zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
  zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
 }
 poly <- vector(mode="list", length(col))
 for(i in seq(poly)){
  poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
 }
 xaxt <- ifelse(horiz, "s", "n")
 yaxt <- ifelse(horiz, "n", "s")
 if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
 if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
 if(missing(xlim)) xlim=XLIM
 if(missing(ylim)) ylim=YLIM
 plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
 for(i in seq(poly)){
  if(horiz){
   polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
  }
  if(!horiz){
   polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
  }
 }
}

png(paste(outname, "Opposite_genes2.png", sep="_"), width=6*1.5, height=6, units="in", res=300)
require(RColorBrewer)
par(layout(rbind(c(1,2), c(1,2)), widths=c(6,1)))
cols <- rev(brewer.pal(7,"RdBu"))
image(t(round(mat, digits=1)), col=cols, xaxt=NULL, yaxt=NULL, axes=FALSE)
axis(1, at=seq(from=0, to=1, length=ncol(mat)), labels=colnames(mat), las=2)
axis(2, at=seq(from=0, to=1, length=nrow(mat)), labels=rownames(mat), las=2)
image.scale(as.vector(mat), col=cols, horiz=FALSE)
dev.off()

pseudobulks <- readRDS(paste(outname, "pseudobulks.rds", sep="_"))

png(paste(outname, "Edge_signature.png", sep="_"), width=10*1.2, height=10, units="in", res=300)
FeaturePlot(obj, c("SQSTM1", "IL32", "APOE", "MT1G", "HLA-A", "ANGPTL8", "AKR1B10", "CES1", "APOA2"))
dev.off()

# Number of DE genes vs cell-type
n_norm_down <- sapply(all_disease_de, function(x) {sum(x$logFC < 0)})
n_norm_up <- sapply(all_disease_de, function(x) {sum(x$logFC > 0)})

type_names <- names(all_disease_de);
my_xlim=c(max(c(n_norm_up,n_norm_down))*-1, max(c(n_norm_up,n_norm_down)))
png(paste(outname, "nDE_barplot.png", sep="_"), width=5, height=5, units="in", res=300)
par(mar=c(4,6.5,1,1))
xloc <- barplot(n_norm_down, names=type_names, xlim=my_xlim, horiz=TRUE, las=1, col="red")
barplot(-1*n_norm_up, names=type_names, xlim=my_xlim, horiz=TRUE, add=TRUE, las=1, col="blue",
		xlab="N DE genes in PSC vs Control")
dev.off()

# Pathway Analysis

source("~/scripts/LiverMap2.0/My_R_Scripts.R")
require(fgsea)
require(gProfileR)

immune_path <- gmtPathways("/cluster/projects/macparland/TA/ExternalData/MSigDb_immune_signatures_c7.all.v7.1.symbols.gmt")
Hallmark_path <- gmtPathways("/cluster/projects/macparland/TA/ExternalData/MSigDbHalmarkPathways.gmt")
MSigAll <- gmtPathways("/cluster/projects/macparland/TA/ExternalData/MSigDb_curated_c2.all.v7.1.symbols.gmt")
reactome <- gmtPathways("/cluster/projects/macparland/TA/ExternalData/ReactomePathways.gmt")
BaderMSig <- gmtPathways("/cluster/projects/macparland/TA/ExternalData/BaderLab25Aug2020/Human_MSigdb_August_01_2020_symbol.gmt.txt")
BaderReact <- gmtPathways("/cluster/projects/macparland/TA/ExternalData/BaderLab25Aug2020/Human_Reactome_August_01_2020_symbol.gmt.txt")


extract_pathway_dat <- function(res, this_pathway, this_name) {
        genes <- unlist(res$rich[unlist(res$rich[,1]) == this_pathway,8])
        names(genes) <- rep(this_name, length(genes));

        bar_point <- unlist(abs(log10(res$rich[unlist(res$rich[,1]) == this_pathway,"padj"]))); names(bar_point)<- this_name;
	direction <- sign(paths$rich$NES[unlist(res$rich[,1]) == this_pathway])
        return(list(genes=genes, point=bar_point, dir=direction))
}

type="Cholangiocyte"; set.seed(101);
scores <- all_disease_de[[type]][,1]*-1; names(scores) <- rownames(all_disease_de[[type]]);
paths <- do_fgsea(sort(scores), pathways=BaderMSig)

names(V(paths$graph))

this_names <- rev(c("IFNa Response", "IFNg Response", "Bile Metabolism", "EMT", "Xenobiotic Metabolism", "Coagulation", "ECM", "Matrisome Associated", "Matrisome")) # lrNK
this_names <- rev(c("IL12", "CD8 TCR", "IFNa Response", "Allograft Rejection", "TCR pathway", "Downstream CD8", "IFNg Response", "CXCR4 Pathway", "MYC targets", "Adipogenesis", "Myogenesis", "EMT", "Estrogen response", "Fatty Acid Metabolism", "Bile Metabolism", "Matrisome Associated", "Matrisome", "ECM", "Coagulation", "Xenobiotic Metabolism" )) # cNK
this_names <- c("Matrisome", "Coagulation", "Matrisome Associated", "ECM", "Xenobiotic Metablism", "Estrogen Response", "Bile Acid Metablism", "IFNg Response", "Allograft Rejection") # CD8+T
this_names <- c("ECM", "Coagulation", "Xenobiotic Metabolism","Matrisome associated", "Matrisome") # CD4+T
this_names <- c("Xenobiotic Metabolism", "ECM Regulators", "Coagulation", "Matrisome") # Monocyte
this_names <- c("Inflammatory Response", "Matrisome Associated", "Allograft Rejection", "Matrisome", "NABA Secreted Factors") # C-Hepato
this_names <- c("Matrisome", "Matrisome Associated", "IFNg Response", "Allograft Rejection", "IFNa Response") # P-Hepato
this_names <- c("Matrisome", "Matrisome Associated") # Kupffer
this_names <- c("ECM", "Coagulation", "Xenobiotic Metabolism", "Matrisome associated", "Matrisome", "TNFA Signalling") # cvLSEC

bars <- c();
for(i in 1:length(V(paths$graph))) {
	p = names(V(paths$graph))[i]
	n = this_names[i]
	summary <- extract_pathway_dat(paths, p, n)
	bars <- c(bars, summary$point*summary$dir)
}
png(paste(outname, type, "pathways.png", sep="_"), width=8, height=8, units="in", res=150)
V(paths$graph)$label <- this_names
plot(paths$graph, vertex.color=paths$vertex_col)
dev.off()
png(paste(outname, type, "pathways2.png", sep="_"), width=8, height=8, units="in", res=150)
par(mar=c(4,10,1,1))
barplot(bars, names=this_names, las=2, horiz=TRUE, xlab="log10(p-value)", col=paths$vertex_col)
dev.off()




#Kupffer - downregulation of Complement Cascade
#lrNK - up regulation of KLRC1/2/D1 - inhibitory molecules & HLA-I


tmp<- all_disease_de[["Kupffer"]]
sort(rownames(tmp)[tmp[,1] < 0])

### Cell-Coll Communication ###
# Version Conflicts
#library(CellChat)
#library(patchwork)
#options(stringsAsFactors = FALSE)

#objs <- readRDS(paste(outname, "Objects_forCellCellInterations.rds", sep="_"))

#healthy_cellchat <- createCellChat(object=objs$healthy_5pr, group.by="full_annotation")
#healthy_cellchat@DB <- CellChatDB.human
#healthy_cellchat <- subsetData(healthy_cellchat)
#psc_cellchat <- createCellChat(object=objs$psc, group.by="cell_type")
#psc_cellchat@DB <- CellChatDB.human
#psc_cellchat <- subsetData(psc_cellchat)

#cell_chat_pipeline <- function(cellchat_obj) {
#	cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
#	cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
#	cellchat_obj <- computeCommunProb(cellchat_obj, type="truncatedMean", trim=0.05)
#	cellchat_obj <- computeCommunProbPathway(cellchat_obj, type="truncatedMean", trim=0.05)
#	cellchat_obj <- aggregateNet(cellchat_obj)
#	return(cellchat_obj)
#}

#healthy_cellchat <- cell_chat_pipeline(healthy_cellchat)
#psc_cellchat <- cell_chat_pipeline(psc_cellchat)






