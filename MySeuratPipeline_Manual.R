script_dir = "/cluster/home/tandrews/scripts/LiverMap2.0"
auto_anno_dir = "/cluster/projects/macparland/TA/AutoAnnotation"
source(paste(script_dir, "Setup_autoannotation.R", sep="/"))
#source(paste(script_dir,"Setup_autoannotation.R", sep="/"))

#i <- as.numeric(as.character(commandArgs(trailingOnly=TRUE)))
i = 1

my_dataset_names <- c("C66", "C68", "C69", "C70", "C72")
seeds <- c(6228, 9013, 7738, 3172, 8224)
folders <- c("McGilvray_Sonya__C66_Caudate", 
		"MacParland_Sonya__C68_Total_liver", 
		"MacParland_Sonya__human_c69_liver_19years_male_3pr_v2",
		"C70_with_top_up_sequencing",
		"C72")
root_dir <-"/cluster/projects/macparland/TA/LiverMap2.0/RawData"

mt_filter <- 50
ng_filter <- 100
nc_filter <- 10
nkNN <- 20
npcs <- 20
nhvg <- 2000

set.seed(seeds[i])
name <- my_dataset_names[i]
folder <- paste(root_dir, folders[i], sep="/");
#res <- 5
res <- 3


require(dplyr)
require(Seurat)
require(Matrix)

# Read the data
mydata <- Read10X(data.dir = paste(folder, "filtered_gene_bc_matrices/GRCh38", sep="/"))
mydata <- mydata[,1:1000]

warning(paste(dim(mydata), "inpute dimensions of", name))

myseur <- CreateSeuratObject(counts = mydata, project = name, min.cells = nc_filter, min.features = ng_filter) ### <- ERROR! gene & cell filter are backwards. - Fixed on 20/04/2020

# Mitochondrial filter
myseur[["percent.mt"]] <- PercentageFeatureSet(myseur, pattern = "^MT-")

myseur <- subset(myseur, subset = nFeature_RNA > ng_filter & percent.mt < mt_filter)

# Add metadata
myseur@meta.data$cell_barcode <- colnames(myseur)
myseur@meta.data$donor <- rep(name, ncol(myseur));
myseur@meta.data$cell_ID <- paste(myseur@meta.data$donor, myseur@meta.data$cell_barcode, sep="_");

# sctransform normalization
#require("sctransform")
#norm <- sctransform::vst(Matrix(myseur@assays$RNA@counts), res_clip_range=c(-Inf, Inf), method="nb_fast");

#myseur@assays$RNA@data <- norm$y;
myseur <- NormalizeData(myseur);
# Scale
myseur <- ScaleData(myseur);
# HVG genes (n=2000)
myseur <- FindVariableFeatures(myseur, selection.method = "vst", nfeatures = 2000)
# PCA 
myseur <- RunPCA(myseur, features = VariableFeatures(object = myseur))
ElbowPlot(myseur)

# Clustering
myseur <- FindNeighbors(myseur, dims = 1:npcs)
myseur <- FindClusters(myseur, resolution = res, k.param=nkNN)
# Visualization with TSNE & UMAP
myseur <- RunTSNE(myseur, dims = 1:npcs)
myseur <- RunUMAP(myseur, dims = 1:npcs, parallel=FALSE)
png(paste(name, "_default_tsne.png", sep=""), width=6, height=6, units="in", res=100)
DimPlot(myseur, reduction = "tsne")
dev.off()
png(paste(name, "_default_umap.png", sep=""), width=6, height=6, units="in", res=100)
DimPlot(myseur, reduction = "umap")
dev.off()

#Cell-cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

myseur <- CellCycleScoring(myseur, s.features = s.genes, g2m.features=g2m.genes, set.ident=TRUE)

# AutoAnnotation with scmap
myseur <- run_scmap_seurat(myseur, scmap_ref=map1_ref);

simplify_annotations <- function(annotations) {
	simplified <- as.character(annotations)
	simplified[simplifed %in% c(
		"AntibodysecretingBcells", 
		"MatureBcells")] <- "Bcells"
	simplified[simplifed %in% c(
		"CD3abTcells", "gdTcells1", "gdTcells2")] <- "Tcells"
	simplified[simplifed %in% c(
		"PericentralHep", "UnidentifiedHep", "PeriportalHep",
		"interzonalHep",)] <- "Hepatocyte"
}


general_labs <- simplify_annotations(myseur@meta.data$scmap_cell_anno)
general_labs2 <- simplify_annotations(myseur@meta.data$scmap_cluster_anno)
general_labs[general_labs != general_labs2] <- "ambiguous"
general_labs[general_labs == "unassigned"] <- "ambiguous"
myseur@meta.data$general_labs <- general_labs;
inconsistent <- myseur@meta.data$scmap_cell_anno != myseur@meta.data$scmap_cluster_anno;
myseur@meta.data$consistent_labs <- myseur@meta.data$scmap_cell_anno;
myseur@meta.data$consistent_labs[inconsistent] <- general_labs[inconsistent]


saveRDS(myseur, paste(name,"Anno_SeurObj.rds", sep="_"));

