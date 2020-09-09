require(methods)
script_dir = "/cluster/home/tandrews/scripts/LiverMap2.0"
auto_anno_dir = "/cluster/projects/macparland/TA/AutoAnnotation"

## Auto Annotation Stuff##
source(paste(script_dir, "Colour_Scheme.R", sep="/"))
simplify_annotations <- function(annotations) {
        simplified <- as.character(annotations)
        simplified[simplified %in% c(
                "AntibodysecretingBcells",
                "MatureBcells", "B_cells", "any_B_cell", "Anti_B")] <- "Bcells"
        simplified[simplified %in% c(
                "CD3abTcells", "gdTcells1", "gdTcells2", "gd_T_1", "gd_T_2", "CD3_T")] <- "Tcells"
        simplified[simplified %in% c(
                "PericentralHep", "UnidentifiedHep", "PeriportalHep",
                "interzonalHep", "Hep", "PortalHep2", "PortaHep1", "PortHep", "CentralHep1", "CentralHep")] <- "Hepatocyte"
        return(map_cell_types(simplified));
}

do_annotation <- function(myseur) {

        myseur <- run_scmap_seurat(myseur, scmap_ref=map1_ref);
        marker_anno <- Use_markers_for_anno(myseur@assays$RNA@data, myseur$seurat_clusters);



        general_labs <- simplify_annotations(myseur@meta.data$scmap_cell_anno)
        general_labs2 <- simplify_annotations(myseur@meta.data$scmap_cluster_anno)
        general_labs[general_labs != general_labs2] <- "ambiguous"
        general_labs[general_labs == "unassigned"] <- "ambiguous"
        myseur@meta.data$general_labs <- general_labs;
        inconsistent <- myseur@meta.data$scmap_cell_anno != myseur@meta.data$scmap_cluster_anno;
        myseur@meta.data$consistent_labs <- as.character(myseur@meta.data$scmap_cell_anno);
        myseur@meta.data$consistent_labs[inconsistent] <- as.character(general_labs[inconsistent]);
        myseur@meta.data$marker_labs <- as.character(marker_anno$cell_assign)
        myseur@meta.data$marker_general_labs <- simplify_annotations(as.character(marker_anno$cell_assign))
        return(myseur)
}

source(paste(script_dir, "Setup_autoannotation.R", sep="/"))

args <- as.character(commandArgs(trailingOnly=TRUE))

Params <- read.table(paste(script_dir, "LiverMap_SampleProcessingParams.csv", sep="/"), sep=",", header=T)

if (grepl("rds", args[1])) {

mt_filter <- 50
ng_filter <- 100
nc_filter <- 10
nkNN <- 20
npcs <- 20
nhvg <- 2000
this_seed <- 101
max_cells <- 20000

file = args[1]
f <- unlist(strsplit(file, "/"))
f <- f[length(f)]

f <- unlist(strsplit(f, "[_\\.]"))
name <- paste(f[1:(length(f)-3)], collapse="_");
rds <- file
out_tag = paste(c(f[(length(f)-2):(length(f)-1)], "SeurObj"), collapse="_")
folder <- as.character(Params[Params$Name==name,"Directory_UHN"])

} else {

i <- as.numeric(args[1])

mt_filter <- Params[i,"MTfilter"]
ng_filter <- Params[i,"nGenefilter"]
nc_filter <- Params[i,"nCellfilter"]
nkNN <- Params[i, "kNN"]
npcs <- Params[i, "nPCs"]
nhvg <- Params[i, "nHGV"]
this_seed <- Params[i, "Seed"]
max_cells <- 20000
name <- as.character(Params[i, "Name"])
folder <- as.character(Params[i, "Directory_UHN"])
rds <- as.character(Params[i, "RDS_file"])
out_tag = "Anno_SeurObj"
}

print(paste(name,"_", out_tag, ".rds", sep=""))

if (!file.exists( paste(name,"_", out_tag, ".rds", sep=""))) {
set.seed(this_seed)
#res <- 5
#res <- 3
# chanaged on 13-07-2020
res <- 1

print("Read Data")

require(dplyr)
require(Seurat)
require(Matrix)

if (file.exists(rds)) {
	mydata <- readRDS(rds)
} else {

	my_10xfolder <- paste(folder, "filtered_gene_bc_matrices/GRCh38", sep="/")

	# Read the data
	mydata <- Seurat::Read10X(data.dir = paste(folder, "filtered_gene_bc_matrices/GRCh38", sep="/"))
	print(paste("Stats out:", dim(mydata), "input dimensions of", name))
	# Enforce a maximum number of cells captured based on expected number.
	mydata <- mydata[,Matrix::colSums(mydata) > ng_filter]
	if (ncol(mydata) > 20000) {
		ng_detect <- Matrix::colSums(mydata);
		q_thresh <- quantile(ng_detect, 1-max_cells/length(ng_detect))
		mydata <-  mydata[,ng_detect > q_thresh]
	}
}
myseur <- Seurat::CreateSeuratObject(counts = mydata, project = name, min.cells = nc_filter, min.features = ng_filter) ### <- ERROR! gene & cell filter are backwards. - Fixed on 20/04/2020

# Mitochondrial filter
myseur[["percent.mt"]] <- Seurat::PercentageFeatureSet(myseur, pattern = "^MT-")

myseur <- subset(myseur, subset = nFeature_RNA > ng_filter & percent.mt < mt_filter)

print(paste("Stats out:", dim(myseur), "seurat dimensions of", name))
# Add metadata
myseur@meta.data$cell_barcode <- colnames(myseur)
myseur@meta.data$donor <- rep(name, ncol(myseur));
myseur@meta.data$cell_ID <- paste(myseur@meta.data$donor, myseur@meta.data$cell_barcode, sep="_");

#add corrected table
if (grepl("emptyDrops", rds)) {
	corrected_rds <- gsub("table", "correctedtable", rds)
	corrected_mat <- readRDS(corrected_rds)
	corrected_mat <- corrected_mat[,match(colnames(myseur), colnames(corrected_mat))]
	corrected_mat <- corrected_mat[match(rownames(myseur), rownames(corrected_mat)),]
	myseur@assays$My_Corrected <- corrected_mat;
	myseur@meta.data$corrected_percent.mt <- Matrix::colSums(corrected_mat[
							grepl("MT-",rownames(corrected_mat)),])/
						 Matrix::colSums(corrected_mat)
}


print("Seurat Pipeline")
# sctransform normalization
#require("sctransform")
#norm <- sctransform::vst(Matrix(myseur@assays$RNA@counts), res_clip_range=c(-Inf, Inf), method="nb_fast");

#myseur@assays$RNA@data <- norm$y;
myseur <- Seurat::NormalizeData(myseur);
# Scale
myseur <- Seurat::ScaleData(myseur);
# HVG genes (n=2000)
myseur <- Seurat::FindVariableFeatures(myseur, selection.method = "vst", nfeatures = 2000)
# PCA 
myseur <- Seurat::RunPCA(myseur, features = VariableFeatures(object = myseur))
ElbowPlot(myseur)

# Clustering
myseur <- Seurat::FindNeighbors(myseur, dims = 1:npcs)
myseur <- Seurat::FindClusters(myseur, resolution = res, k.param=nkNN)
# Visualization with TSNE & UMAP
myseur <- Seurat::RunTSNE(myseur, dims = 1:npcs)
myseur <- Seurat::RunUMAP(myseur, dims = 1:npcs, parallel=FALSE)

#Cell-cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

myseur <- Seurat::CellCycleScoring(myseur, s.features = s.genes, g2m.features=g2m.genes, set.ident=TRUE)


# AutoAnnotation with scmap
myseur <- do_annotation(myseur)

saveRDS(myseur, paste(name,"_",out_tag, ".rds", sep=""));
} else {
myseur <- readRDS( paste(name,"_",out_tag, ".rds", sep=""));
}

print(paste("Stats out:", length(unique(myseur@meta.data$seurat_clusters)), "seurat clusters in", name))

##### Make plots ####
agg_coord_by_cluster <- function(coords, clusters) {
	x <- split(seq(nrow(coords)), clusters)
	result <- sapply(x, function(a) apply(coords[a,],2,median))
	return(result)
}

umap_lab_pos <- agg_coord_by_cluster(myseur@reductions$umap@cell.embeddings, myseur@meta.data$seurat_clusters)

require("ggplot2")
require("Seurat")

# UMAP + Ref scmap anno
new_colour_scheme <- Cell_type_colours[order(Cell_type_colours[,1]),]
myseur@meta.data$consistent_labs <- map_cell_types(myseur@meta.data$consistent_labs)
new_colour_scheme <- new_colour_scheme[new_colour_scheme[,1] %in% myseur@meta.data$consistent_labs,]

png(paste(name, out_tag, "_refanno_umap.png", sep="_"), width=6, height=6, units="in", res=100)
DimPlot(myseur, reduction="umap", group.by="consistent_labs", pt.size=.1)+scale_color_manual(values=new_colour_scheme[,2])+annotate("text", x=umap_lab_pos[1,], y=umap_lab_pos[2,], label=colnames(umap_lab_pos), colour="grey35")
dev.off()


# UMAP + Marker scmap anno

new_colour_scheme <- Cell_type_colours[order(Cell_type_colours[,1]),]
myseur@meta.data$marker_labs <- map_cell_types(myseur@meta.data$marker_labs)
new_colour_scheme <- new_colour_scheme[new_colour_scheme[,1] %in% myseur@meta.data$marker_labs,]

png(paste(name, out_tag, "_markanno_umap.png", sep="_"), width=6, height=6, units="in", res=100)
DimPlot(myseur, reduction="umap", group.by="marker_labs", pt.size=.1)+scale_color_manual(values=new_colour_scheme[,2])+annotate("text", x=umap_lab_pos[1,], y=umap_lab_pos[2,], label=colnames(umap_lab_pos), colour="grey35")
dev.off()

png(paste(name, out_tag, "_default_umap.png", sep="_"), width=6, height=6, units="in", res=100)
Seurat::DimPlot(myseur, reduction = "umap")
dev.off()

png(paste(name, out_tag, "perMT.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::FeaturePlot(myseur, "percent.mt", reduction="umap")
dev.off()
png(paste(name, out_tag, "nFeature.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::FeaturePlot(myseur, "nFeature_RNA", reduction="umap")
dev.off()
png(paste(name, out_tag, "CCphase.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::DimPlot(myseur, group.by="Phase", reduction="umap")
dev.off()
png(paste(name, out_tag, "sample.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::DimPlot(myseur, group.by="donor", reduction="umap")
dev.off()

SoupX_outfile <- paste(paste("SoupX/", name, sep=""), out_tag, "SoupX.rds", sep="_")
if (!file.exists(SoupX_outfile)) {
# SoupX
## Need to make the EmptyDrops version of this!!!
require("SoupX")
require("Seurat")
set.seed(this_seed)
# Create SoupX object
raw <- Read10X(paste(folder, "raw_gene_bc_matrices/GRCh38", sep="/"))
raw <- raw[,Matrix::colSums(raw) > 10]
raw <- raw[rownames(raw) %in% rownames(myseur),]
keep_cells <- colnames(myseur);
tot_umi <- Matrix::colSums(raw);
tot_is_cell <- min(tot_umi[colnames(raw) %in% colnames(myseur)])
raw_is_cell <- raw[,keep_cells]
myseur <- myseur[,match(colnames(raw_is_cell), colnames(myseur))]

sc <- SoupChannel(raw, raw_is_cell, metaData=myseur@meta.data, soupRange=c(10, tot_is_cell*0.9))

anno_cluster <- cell_anno_to_cluster_anno(as.character(unlist(myseur@meta.data$consistent_labs)), myseur@meta.data$seurat_clusters)
anno_cluster <- anno_cluster$lab[match(myseur@meta.data$seurat_clusters, anno_cluster$cluster)]

sc <- setClusters(sc, anno_cluster);
sc <- setDR(sc, cbind(myseur@reductions$umap@cell.embeddings[,1], myseur@reductions$umap@cell.embeddings[,2]));

# Use known genes to estimate contamination fraction per cell
non_expressed_gene_list = list(HB =c("HBB", "HBA1", "HBA"), 
				TCR=c("TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2"), 
				Drug=c("CYP2E1", "CYP3A7", "CYP2A7", "CYP2A6"), 
				BCR=c("IGKC", "IGHE", "IGHM", "IGLC1", "IGLC3", "IGLC2"), 
				Clot=c("F10", "F5", "F9", "F7", "FGB"), 
				Bile=c("SLC10A1", "SLC01B1", "SLCO1B3", "SLC22A1", "SLC22A7", 
					"ABCB1", "ABCB11", "ABCB4", "ABCC3", "ABCC6"));
non_expressed_gene_list <- lapply(non_expressed_gene_list, function(x){return(x[x %in% rownames(myseur)])})

useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = non_expressed_gene_list, clusters=NULL)
sc = calculateContaminationFraction(sc, non_expressed_gene_list, useToEst = useToEst, cellSpecificEstimates=TRUE)
quantile(sc$metaData$rho)

# Correcting the expression matrix.
out = adjustCounts(sc)
out <- out[,match(colnames(myseur), colnames(out))]

myseur@assays$SoupCorrected <- out;
saveRDS(myseur, SoupX_outfile);
} else {
myseur <- readRDS(SoupX_outfile)
}


#myseur <- do_annotation(myseur)
#orig_labels <- myseur@meta.data

labels <- list()
agreement <- c()
myseur@assays$Orig_Counts <- myseur@assays$RNA@counts

# Compare autoannotation across corrections.
if (length(myseur@assays) >2) {
	for (i in 2:length(myseur@assays)) {
		this_corr_mat <- myseur@assays[[i]]
		myseur<- myseur[rownames(myseur) %in% rownames(this_corr_mat),]
		myseur<- myseur[, colnames(myseur) %in% colnames(this_corr_mat)]
		this_corr_mat <- this_corr_mat[match(rownames(myseur), rownames(this_corr_mat)),]
		this_corr_mat <- this_corr_mat[,match(colnames(myseur), colnames(this_corr_mat))]
		
		myseur@assays$RNA@counts <- this_corr_mat
		myseur <- Seurat::NormalizeData(myseur);
		myseur <- Seurat::ScaleData(myseur);
		myseur <- do_annotation(myseur);
		corr_name <- names(myseur@assays)[i];
		labels[[corr_name]] <- myseur@meta.data
		# plot
		png(paste(name, out_tag, corr_name, "expression.png", sep="_"), width=6, height=6, units="in", res=100)
		Seurat::FeaturePlot(myseur, features=c("ALB","APOB", "SLC22A1","HBB", "HBB", "APOA1", "MARCO", "TNF"))
		dev.off()

		myseur[["percent.mt"]] <- Seurat::PercentageFeatureSet(myseur, pattern = "^MT-")
		myseur <- subset(myseur, subset = nFeature_RNA > ng_filter & percent.mt < mt_filter)
		png(paste(name, out_tag, corr_name, "MT.png", sep="_"), width=6, height=6, units="in", res=100)
		Seurat::FeaturePlot(myseur, features=c("percent.mt"))
		dev.off()

		#agreement of autoannotation labels
		is.lab <- myseur@meta.data[,"general_labs"] != "ambiguous" & 
			  myseur@meta.data[,"marker_general_labs"] != "ambiguous" & 
			  myseur@meta.data[,"marker_general_labs"] != "Hepatocyte"
		score <- sum(myseur@meta.data[is.lab,"general_labs"] == myseur@meta.data[is.lab,"marker_general_labs"])/sum(is.lab)
		agreement[corr_name] <- score;
		
	}
}

png(paste(name, out_tag, corr_name, "Score.png", sep="_"), width=8, height=4, units="in", res=100)
par(mar=c(4,10,1,1))
barplot(agreement, horiz=TRUE, xlab="Autoannotation Agreement", col="grey65")
dev.off()

print(agreement)
