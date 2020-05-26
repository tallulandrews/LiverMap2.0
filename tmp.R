require(methods)
script_dir = "/cluster/home/tandrews/scripts/LiverMap2.0"
auto_anno_dir = "/cluster/projects/macparland/TA/AutoAnnotation"
source(paste(script_dir, "Setup_autoannotation.R", sep="/"))

#i <- as.numeric(as.character(commandArgs(trailingOnly=TRUE)))
#i = 31

Params <- read.table(paste(script_dir, "LiverMap_SampleProcessingParams.csv", sep="/"), sep=",", header=T)

mt_filter <- Params[i,"MTfilter"]
ng_filter <- Params[i,"nGenefilter"]
nc_filter <- Params[i,"nCellfilter"]
nkNN <- Params[i, "kNN"]
npcs <- Params[i, "nPCs"]
nhvg <- Params[i, "nHGV"]
this_seed <- Params[i, "Seed"]
max_cells <- 2000000
name <- paste(as.character(Params[i, "Name"], "nolimit", sep="_"))
folder <- as.character(Params[i, "Directory_UHN"])
rds <- as.character(Params[i, "RDS_file"])

set.seed(this_seed)
#res <- 5
res <- 3


require(dplyr)
require(Seurat)
require(Matrix)


	my_10xfolder <- paste(folder, "filtered_gene_bc_matrices/GRCh38", sep="/")

	# Read the data
	mydata <- Read10X(data.dir = paste(folder, "filtered_gene_bc_matrices/GRCh38", sep="/"))
	print(paste("Stats out:", dim(mydata), "input dimensions of", name))
	# Enforce a maximum number of cells captured based on expected number.
	mydata <- mydata[,Matrix::colSums(mydata) > ng_filter]

myseur <- CreateSeuratObject(counts = mydata, project = name, min.cells = nc_filter, min.features = ng_filter) ### <- ERROR! gene & cell filter are backwards. - Fixed on 20/04/2020

# Mitochondrial filter
myseur[["percent.mt"]] <- PercentageFeatureSet(myseur, pattern = "^MT-")

myseur <- subset(myseur, subset = nFeature_RNA > ng_filter & percent.mt < mt_filter)

print(paste("Stats out:", dim(myseur), "seurat dimensions of", name))
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

print(paste("Stats out:", length(unique(myseur@meta.data$seurat_clusters)), "seurat clusters in", name))

#Cell-cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

myseur <- CellCycleScoring(myseur, s.features = s.genes, g2m.features=g2m.genes, set.ident=TRUE)

# AutoAnnotation with scmap
myseur <- run_scmap_seurat(myseur, scmap_ref=map1_ref);

simplify_annotations <- function(annotations) {
	simplified <- as.character(annotations)
	simplified[simplified %in% c(
		"AntibodysecretingBcells", 
		"MatureBcells")] <- "Bcells"
	simplified[simplified %in% c(
		"CD3abTcells", "gdTcells1", "gdTcells2")] <- "Tcells"
	simplified[simplified %in% c(
		"PericentralHep", "UnidentifiedHep", "PeriportalHep",
		"interzonalHep")] <- "Hepatocyte"
	return(simplified);
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

