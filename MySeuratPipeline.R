my_metadata_table <- read.table("Metadata20LiverMapPlusParams.csv", sep=",", header=T, stringsAsFactors=FALSE);
source("Setup_autoannotation.R")

all_anno <- list();
for(dataset_row in 1:21) {

set.seed(my_metadata_table$Seed[dataset_row])
name <- my_metadata_table$Name[dataset_row]
folder <- my_metadata_table$Directory[dataset_row]
mt_filter <- my_metadata_table$MTfilter[dataset_row]
ng_filter <- my_metadata_table$nGenefilter[dataset_row]
nc_filter <- my_metadata_table$nCellfilter[dataset_row]
nhvg <- my_metadata_table$nHVG[dataset_row]
npcs <- my_metadata_table$nPCs[dataset_row]
nkNN <- my_metadata_table$kNN[dataset_row]
res <- 5


require(dplyr)
require(Seurat)
require(Matrix)

# Read the data
mydata <- Read10X(data.dir = paste(folder, "filtered_gene_bc_matrices/GRCh38", sep="/"))
myseur <- CreateSeuratObject(counts = mydata, project = name, min.cells = ng_filter, min.features = nc_filter) ### <- ERROR! gene & cell filter are backwards.

# Mitochondrial filter
myseur[["percent.mt"]] <- PercentageFeatureSet(myseur, pattern = "^MT-")

myseur <- subset(myseur, subset = nFeature_RNA > ng_filter & percent.mt < mt_filter)

# sctransform normalization
require("sctransform")
norm <- sctransform::vst(Matrix(myseur@assays$RNA@counts), res_clip_range=c(-Inf, Inf), method="nb_fast");

myseur@assays$RNA@data <- norm$y;
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

# AutoAnnotation with scmap
myseur <- run_scmap_seurat(myseur, scmap_ref=map1_ref);
saveRDS(myseur, paste(name,"Anno_SeurObj.rds", sep="_"));

anno_tab <- myseur@meta.data
anno_tab$cell_barcode <- colnames(myseur);

rownames(anno_tab) <- paste(anno_tab$orig.ident, anno_tab$cell_barcode, sep="_")

all_anno[[name]] <- anno_tab;

}
saveRDS(all_anno, "All20_automatedannotation.rds");



