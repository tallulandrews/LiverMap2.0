##### NOTE ####
# Requires sample name and folder of Seurat raw matrix as input.
# Pipeline:
# EmptyDrops - select cells
# Seurat - quick standard pipeline & clustering
# scmap - automatic annotation using map 1.0
# SoupX - remove background RNA
# Last Updated: 28/07/2020 by: Tallulah Andrews
# Change log:
# Created by: Tallulah Andrews (28/07/2020)
##############

###### Load Requirements ######
a1 = suppressPackageStartupMessages(require(methods))
a2 = suppressPackageStartupMessages(require(Matrix))
a3 = suppressPackageStartupMessages(require("DropletUtils"))
a4 = suppressPackageStartupMessages(require(dplyr))
a5 = suppressPackageStartupMessages(require(Seurat))
a6 = suppressPackageStartupMessages(require(optparse))
a7 = suppressPackageStartupMessages(require("SoupX"))
a8 = suppressPackageStartupMessages(require("ggplot2"))
if (! (a1 & a2 & a3 & a4 & a5 & a6 & a7 & a8) ) {
	print("Cannot load all required packages: methods, Matrix, DropletUtils, dplyr, Seurat, optparse, SoupX", "ggplot2")
	exit()
}

script_dir = "/cluster/home/tandrews/scripts/LiverMap2.0"
source(paste(script_dir, "My_R_Scripts.R", sep="/"));
source("/cluster/home/tandrews/R-Scripts/Ensembl_Stuff.R"); #Map gene names
## Auto Annotation Stuff ##
# Colour Scheme #
source(paste(script_dir, "Colour_Scheme.R", sep="/"))

source(paste(script_dir, "Setup_autoannotation.R", sep="/"))
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

## ARGUMENTS & Parameters ##


option_list <- list(
    make_option(c("-i", "--input_dir"),
        help="Directory of raw unfiltered cell ranger UMI counts (required)"),
    make_option(c("-o", "--out_prefix"), default="Cleaned_output",
        help="prefix for output files [default %default]"),
    make_option(c("--organism"), default="Human",
        help="Organism, currently supported: [Human, Mouse, Rat]"),
    make_option(c("-M", "--max_cells"), type="integer", default=20000,
        help="Maximum number of plausible cells [default %default]",
        metavar="number"),
    make_option(c("-m", "--min_cells"), type="integer", default=100,
        help="Minimum number of plausible cells [default %default]",
        metavar="number"),
    make_option(c("--FDR"), default=0.01,
        help="FDR threshold for EmptyDrops [default %default]",
        metavar="false discovery rate"),
    make_option(c("--trim"), default=10,
        help="threshold below which droplets are ignored completely [default %default]",
        metavar="false discovery rate"),
    make_option(c("--mt"), default=50,
        help="threshold for mitochondrial expression percentage [default %default]",
        metavar="number"),
    make_option(c("-g", "--genes"), default=100,
        help="Threshold for minimum genes per cell [default %default]",
        metavar="number"),
    make_option(c("-c", "--cells"), default=10,
        help="Threshold for minimum cells a gene is detected in [default %default]"),
    make_option(c("--npcs"), default=20,
        help="number of principal components to use [default %default]"),
    make_option(c("--kNN"), default=20,
        help="k used for the nearest neighbour network [default %default]"),
    make_option(c("--res"), default=1,
        help="Resolution used for Seurat clustering [default %default]")
    )

OPTS <- parse_args(OptionParser(option_list=option_list))
print(OPTS)


rawdata <- Read10X(data.dir = OPTS$input_dir)
# Remove rows & columns that are completely zero
rawdata <- rawdata[,Matrix::colSums(rawdata) > OPTS$trim]
rawdata <- rawdata[Matrix::rowSums(rawdata) > 0,]
print(paste("Raw Data:", dim(rawdata), "input",c("gene","cells"), "of", OPTS$out_prefix))

# fixed edgecase when fewer droplets pass trim threshold than 
# specified for minimum plausible cells (20Aug2020)
if( ncol(rawdata) < OPTS$max_cells) {
	OPTS$max_cells <- ncol(rawdata) - 2;
}

# Rank droplet barcodes by total UMIs
br.out <- barcodeRanks(rawdata)

## Get total UMIs that correspond to the min & max thresholds.
plausible_cell_threshold <- max(br.out$total[br.out$rank > OPTS$max_cells]);
mandatory_cell_threshold <- max(br.out$total[br.out$rank < OPTS$min_cells]);

# My calling threshold
n_umi_sorted <- br.out$total[order(br.out$total, decreasing = T)]
rank_sorted <- 1:length(n_umi_sorted);

slope <- diff(log(n_umi_sorted))/diff(log(rank_sorted))
spar = 0.5
inflection = OPTS$min_cells
while(inflection == OPTS$min_cells) { 
	smooth_slope <- smooth.spline(slope, spar=spar) # changed to 0.4 from 0.5 on 19Aug202
	inflection <- which(smooth_slope$y == min(smooth_slope$y[OPTS$min_cells:OPTS$max_cells]))
	spar=spar*0.8
}
spar = spar/0.8;
my_inflection <- n_umi_sorted[inflection];

is_empty_threshold = max(plausible_cell_threshold, 
			br.out$knee, 
			br.out$inflection, my_inflection)
# fixed the comparison below on 19 Aug2020
if (is_empty_threshold > n_umi_sorted[OPTS$min_cells*10]){
	if (is.finite(plausible_cell_threshold)) {
		is_empty_threshold = plausible_cell_threshold
	} else {
		is_empty_threshold = min(br.out$knee, br.out$inflection)
	}
}

# Threshold Plot
png(paste(OPTS$out_prefix, "inflection_points.png", sep="_"), width=4, height=4, units="in", res=50)
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")
abline(h=br.out$knee, col="darkorange", lty=2, lwd=2)
abline(h=br.out$inflection, col="firebrick", lty=2, lwd=2)
abline(h=my_inflection, col="darkgoldenrod2", lty=2, lwd=2)
abline(h=is_empty_threshold, col="dodgerblue", lty=2, lwd=2)
legend("bottomleft", lty=2, lwd=2, c("knee", "inflection", "log-inflection", "empty_threshold"), col=c("darkorange", "firebrick", "darkgoldenrod2", "dodgerblue"), bty="n");
dev.off()



########### Begin if init seurat output file exists
if (!file.exists(paste(OPTS$out_prefix, "EmptyOnly.rds", sep="_"))) {
# Run EmptyDrops
set.seed(100)
e.out <- emptyDrops(rawdata, lower=is_empty_threshold, niters=100000, ignore=OPTS$trim, retain=mandatory_cell_threshold)

# Clean up results
e.out <- e.out[!is.na(e.out$PValue),] # Remove NAs (too few UMIs to estimate PValue)
e.out$q.value <- e.out$PValue;
e.out$q.value[e.out$Limited] <- 0; # Set those with p < 1/n_simulations to p = 0 to ensure they are kept after multiple testing correction
e.out$q.value <- p.adjust(e.out$q.value, method="fdr"); # apply FDR multiple testing correction.

# Significance threshold
is.cell <- e.out$q.value <= OPTS$FDR


# Plot of selected cells
png(paste(OPTS$out_prefix, "emptyDrops.png", sep="_"), width=6, height=6, units="in", res=150)
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
    xlab="Total UMI count", ylab="-Log Probability", pch=18)
dev.off()

# Subset the matrix to the selected droplets and save it to a file.
# Get number of detected genes/droplet
emptymat <- rawdata[,match(rownames(e.out),colnames(rawdata))]
emptymat <- emptymat[,is.cell]

print(paste("EmptyDrops :", dim(emptymat), "output",c("gene","cells"), "of", OPTS$out_prefix))

############### Quick Seurat Pipeline ##############
set.seed(8478)

myseur <- Seurat::CreateSeuratObject(counts = emptymat, project = OPTS$out_prefix, min.cells = OPTS$cells, min.features = OPTS$genes) 

# Mitochondrial filter
myseur[["percent.mt"]] <- Seurat::PercentageFeatureSet(myseur, pattern = "^MT-", )
if (sum(myseur[["percent.mt"]]) == 0) {
	myseur[["percent.mt"]] <- Seurat::PercentageFeatureSet(myseur, pattern = "^mt-", ) # Rat
}
if (sum(myseur[["percent.mt"]]) == 0) {
	myseur[["percent.mt"]] <- Seurat::PercentageFeatureSet(myseur, pattern = "^Mt-", ) # Mouse
}

myseur <- subset(myseur, subset = nFeature_RNA > OPTS$genes & percent.mt < OPTS$mt)

print(paste("Seurat :", dim(myseur),  "output",c("gene","cells"), "of", OPTS$out_prefix))

# Remap Genes using Orthologs
rename_genes <- function(myseur, new_names, old_names) {
	new_names[new_names==""] <- old_names[new_names==""]
	new_names[duplicated(new_names)] <- old_names[duplicated(new_names)]
	rownames(myseur@assays$RNA@counts) <- new_names
	if (nrow(myseur@assays$RNA@data)==length(new_names)) {
		rownames(myseur@assays$RNA@data) <- new_names
	}
	if (nrow(myseur@assays$RNA@scale.data)==length(new_names)) {
		rownames(myseur@assays$RNA@scale.data) <- new_names
	}
	return(myseur)
}
	
genes <- rownames(myseur)
myseur[["RNA"]][["origID"]] <- genes
#myseur@misc[["orig_gene_ID"]] <- genes
raw_genes <- rownames(rawdata) # fix gene ID to remove hypens like Seurat does.
raw_genes <- sub("_","-", raw_genes)
rownames(rawdata) <- raw_genes
if (OPTS$organism == "Mouse") {
	hgenes <- General_Map(genes, in.org="Mmus", out.org="Hsap", in.name="symbol", out.name="symbol")
	myseur <- rename_genes(myseur, hgenes, genes)
	hgenes <- General_Map(rownames(rawdata), in.org="Mmus", out.org="Hsap", in.name="symbol", out.name="symbol")
	tofix <- hgenes=="" | duplicated(hgenes)
	hgenes[tofix] <- rownames(rawdata)[tofix]
	rownames(rawdata) <- hgenes
}
if (OPTS$organism == "Rat") {
	hgenes <- General_Map(genes, in.org="Rat", out.org="Hsap", in.name="symbol", out.name="symbol")
	myseur <- rename_genes(myseur, hgenes, genes)
	hgenes <- General_Map(rownames(rawdata), in.org="Rat", out.org="Hsap", in.name="symbol", out.name="symbol")
	tofix <- hgenes=="" | duplicated(hgenes)
	hgenes[tofix] <- rownames(rawdata)[tofix]
	rownames(rawdata) <- hgenes
}


# Add metadata
myseur@meta.data$cell_barcode <- colnames(myseur)
myseur@meta.data$donor <- rep(OPTS$out_prefix, ncol(myseur));
myseur@meta.data$cell_ID <- paste(myseur@meta.data$donor, myseur@meta.data$cell_barcode, sep="_");

orig.meta.data <- myseur@meta.data;
}
###### End if init seurat doesn't exist


print("Seurat Pipeline")
run_seurat_pipeline <- function(myseur, out_tag) {
	# Normalize
	myseur <- Seurat::NormalizeData(myseur);
	# Scale
	myseur <- Seurat::ScaleData(myseur);
	# HVG genes (n=2000)
	myseur <- Seurat::FindVariableFeatures(myseur, selection.method = "vst", nfeatures = 2000)
	hvgs <- VariableFeatures(myseur); hvgs <- hvgs[!grepl("^MT-", hvgs)]; # added 22Sept2020
	hvgs <- VariableFeatures(myseur); hvgs <- hvgs[!grepl("^Mt-", hvgs)]; # added 22Sept2020
	hvgs <- VariableFeatures(myseur); hvgs <- hvgs[!grepl("^mt-", hvgs)]; # added 22Sept2020
	VariableFeatures(myseur) <- hvgs;
	# PCA 
	myseur <- Seurat::RunPCA(myseur, features = VariableFeatures(object = myseur))

	# Clustering
	myseur <- Seurat::FindNeighbors(myseur, dims = 1:OPTS$npcs)
	myseur <- Seurat::FindClusters(myseur, resolution = OPTS$res, k.param=OPTS$kNN)
	# Visualization with TSNE & UMAP
	myseur <- Seurat::RunUMAP(myseur, dims = 1:OPTS$npcs, parallel=FALSE)

	#Cell-cycle
	s.genes <- cc.genes$s.genes
	g2m.genes <- cc.genes$g2m.genes

	myseur <- Seurat::CellCycleScoring(myseur, s.features = s.genes, g2m.features=g2m.genes, set.ident=TRUE)

	# AutoAnnotation with scmap
	myseur <- do_annotation(myseur)


	print(paste("Clusters:", length(unique(myseur@meta.data$seurat_clusters)), "in", OPTS$out_prefix))

	############### PLOTTING ###########

	##### Make plots ####
	agg_coord_by_cluster <- function(coords, clusters) {
		x <- split(seq(nrow(coords)), clusters)
		result <- sapply(x, function(a) apply(coords[a,],2,median))
		return(result)
	}

	umap_lab_pos <- agg_coord_by_cluster(myseur@reductions$umap@cell.embeddings, myseur@meta.data$seurat_clusters)

	# UMAP + Ref scmap anno
	new_colour_scheme <- Cell_type_colours[order(Cell_type_colours[,1]),]
	myseur@meta.data$consistent_labs <- map_cell_types(myseur@meta.data$consistent_labs)
	new_colour_scheme <- new_colour_scheme[new_colour_scheme[,1] %in% myseur@meta.data$consistent_labs,]

	png(paste(OPTS$out_prefix, out_tag, "_refanno_umap.png", sep="_"), width=6, height=6, units="in", res=100)
	DimPlot(myseur, reduction="umap", group.by="consistent_labs", pt.size=.1)+scale_color_manual(values=new_colour_scheme[,2])+annotate("text", x=umap_lab_pos[1,], y=umap_lab_pos[2,], label=colnames(umap_lab_pos), colour="grey35")
	dev.off()

	# UMAP + Marker scmap anno

	new_colour_scheme <- Cell_type_colours[order(Cell_type_colours[,1]),]
	myseur@meta.data$marker_labs <- map_cell_types(myseur@meta.data$marker_labs)
	new_colour_scheme <- new_colour_scheme[new_colour_scheme[,1] %in% myseur@meta.data$marker_labs,]
	
	png(paste(OPTS$out_prefix, out_tag, "_markanno_umap.png", sep="_"), width=6, height=6, units="in", res=100)
	print(DimPlot(myseur, reduction="umap", group.by="marker_labs", pt.size=.1)+scale_color_manual(values=new_colour_scheme[,2])+annotate("text", x=umap_lab_pos[1,], y=umap_lab_pos[2,], label=colnames(umap_lab_pos), colour="grey35"))
	dev.off()

	# General QC Plots
	png(paste(OPTS$out_prefix, out_tag, "_default_umap.png", sep="_"), width=6, height=6, units="in", res=100)
	print(Seurat::DimPlot(myseur, reduction = "umap"))
	dev.off()

	png(paste(OPTS$out_prefix, out_tag, "perMT.png", sep="_"), width=6, height=6, units="in", res=150)
	print(Seurat::FeaturePlot(myseur, "percent.mt", reduction="umap")+theme(
		axis.line=element_blank(),axis.text.x=element_blank(),
	          axis.text.y=element_blank(),axis.ticks=element_blank(),
        	  axis.title.x=element_blank(),  axis.title.y=element_blank()))
	dev.off()
	png(paste(OPTS$out_prefix, out_tag, "nFeature.png", sep="_"), width=6, height=6, units="in", res=150)
	print(Seurat::FeaturePlot(myseur, "nFeature_RNA", reduction="umap")+theme(
                axis.line=element_blank(),axis.text.x=element_blank(),
                  axis.text.y=element_blank(),axis.ticks=element_blank(),
                  axis.title.x=element_blank(),  axis.title.y=element_blank()))
	dev.off()
	png(paste(OPTS$out_prefix, out_tag, "CCphase.png", sep="_"), width=6, height=6, units="in", res=150)
	print(Seurat::DimPlot(myseur, group.by="Phase", reduction="umap"))
	dev.off()

	png(paste(OPTS$out_prefix, out_tag, "MarkerGenes.png", sep="_"), width=12, height=12, units="in", res=150)
	print(Seurat::FeaturePlot(myseur, features=c("ALB","CYP3A4", "SCD", "MARCO", "LYZ", "CD68", "IGKC", "TRAC", "PTPRC"))+theme(
                axis.line=element_blank(),axis.text.x=element_blank(),
                  axis.text.y=element_blank(),axis.ticks=element_blank(),
                  axis.title.x=element_blank(),  axis.title.y=element_blank()))
	dev.off()

	return(myseur)
}

if (!file.exists(paste(OPTS$out_prefix, "EmptyOnly.rds", sep="_"))) {
	myseur <- run_seurat_pipeline(myseur, "initSeurat");

	saveRDS(myseur, paste(OPTS$out_prefix, "EmptyOnly.rds", sep="_"));

} else {
	myseur <- readRDS(paste(OPTS$out_prefix, "EmptyOnly.rds", sep="_"));
	orig.meta.data <- myseur@meta.data;
}
############## SoupX ##############

SoupX_outfile <- paste(OPTS$out_prefix, "SoupX.rds", sep="_")
# SoupX
require("Seurat")
set.seed(4671)
# Create SoupX object
#my_seur_genes <- unlist(myseur[["RNA"]][["origID"]])
#my_seur_genes <- unlist(myseur@misc[["orig_gene_ID"]])
#raw_data_genes <- rownames(rawdata)
rawdata <- rawdata[rownames(rawdata) %in% rownames(myseur),]
myseur <- myseur[rownames(myseur) %in% rownames(rawdata),]
#rawdata <- rawdata[match(unlist(myseur[["RNA"]][["origID"]]), rownames(rawdata)),]
#rawdata <- rawdata[match(unlist(myseur@misc[["orig_gene_ID"]]), rownames(rawdata)),]
rawdata <- rawdata[match(rownames(myseur), rownames(rawdata)),]
#rownames(rawdata) <- rownames(myseur)
keep_cells <- colnames(myseur);
tot_umi <- Matrix::colSums(rawdata);
tot_is_cell <- min(tot_umi[colnames(rawdata) %in% colnames(myseur)])
raw_is_cell <- rawdata[,keep_cells]
myseur <- myseur[,match(colnames(raw_is_cell), colnames(myseur))]

sc <- SoupChannel(rawdata, raw_is_cell, metaData=myseur@meta.data, soupRange=c(10, tot_is_cell*0.9))

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
#Jawaria - 7Aug2020
#expanded_expressed_gene_list = list(HB =c("HBB", "HBA1", "HBA"), 
#				TCR=c("TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2"), 
#				Drug=c("CYP2E1", "CYP3A7", "CYP2A7", "CYP2A6", "CYP2B6", "CYP2C8"), 
#				BCR=c("IGKC", "IGHE", "IGHM", "IGLC1", "IGLC3", "IGLC2"), 
#				Clot=c("F10", "F5", "F9", "F7", "FGB", "FGG", "FGL1"), 
#				Bile=c("SLC10A1", "SLC01B1", "SLCO1B3", "SLC22A1", "SLC22A7", 
#					"ABCB1", "ABCB11", "ABCB4", "ABCC3", "ABCC6"),
#				Lipid=c("APOC1", "APOC3", "APOA2", "APOA1", "APOE", "APOH"));

# For snRNAseq
expanded_non_genes = list(
		AntiB=c("IGKC","JCHAIN","IGHA1","IGLC1","IGLC2","IGLC3"),
		MatB=c("CD22","CD37","CD79B","FCRL1","LTB","DERL3","IGHG4"),
		CD3T=c("CD8A","CD8B","CD3D","CD3G","TRAC","IL32","TRBC1","TRBC2"),
		Hep=c("CYP1A2","CYP2E1","CYP3A4","GLUL","DCXR","FTL","GPX2","GSTA1","CYP2A7","FABP1","HAL","AGT","ALDOB","SDS"),
		LSEC=c("FCN2","CLEC1B","CLEC4G","PVALB","S100A13","GJA5","SPARCL1","CLEC14A","PLVAP","EGR3"),
		Eryth=c("HBB","HBA1","HBA2"),
		NKT=c("CSTW","IL7R","GZMB","GZMH","TBX21","HOPX","PRF1","S100B","TRDC","TRGC1","TRGC2","IL2RB","KLRB1","NCR1","NKG7","NCAM1","XCL2","XCL1","CD160","KLRC1"),
		Mac=c("VCAN","S100A8","MNDA","LYZ","FCN1","CXCL8","VCAN","VCAM1","TTYH3","TIMD4","SLC40A1","RAB31","MARCO","HMOX1","C1QC"),
		Chol=c("PROM1","SOX9","KRT7","KRT19","CFTR","EPCAM","CLDN4","CLDN7","ANXA4","TACSTD2"),
		Stel=c("ACTA2","COL1A1","RBP1","TAGLN","ADAMTSL2","GEM","LOXL1","LUM"),
		Endo=c("PECAM1","TAGLN","VWF","FLT1","MMRN1","RSPO3","LYPD2","LTC4S","TSHZ2","IL1R1")
		)
		

non_expressed_gene_list <- lapply(non_expressed_gene_list, function(x){return(x[x %in% rownames(myseur)])})
size <- unlist(lapply(non_expressed_gene_list, length))
non_expressed_gene_list<-non_expressed_gene_list[size > 3]

expanded_non_genes <- lapply(expanded_non_genes, function(x){return(x[x %in% rownames(myseur)])})
size <- unlist(lapply(expanded_non_genes, length))
expanded_non_genes<-expanded_non_genes[size > 3]

if (length(non_expressed_gene_list) < 4) {
	non_expressed_gene_list <- expanded_non_genes
}

useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = non_expressed_gene_list, clusters=myseur$seurat_clusters)
sc = calculateContaminationFraction(sc, non_expressed_gene_list, useToEst = useToEst, cellSpecificEstimates=TRUE)
quantile(sc$metaData$rho)

# Correcting the expression matrix.
out = adjustCounts(sc)
out <- out[,match(colnames(myseur), colnames(out))]

print(identical(colnames(out), rownames(myseur@meta.data)))

soup_seurat <- CreateSeuratObject(counts=out, meta.data=orig.meta.data, project=OPTS$out_prefix, min.cells = OPTS$cells, min.features = OPTS$genes);

saveRDS(soup_seurat, paste(OPTS$out_prefix, "SoupX_prepipeline.rds", sep="_"))

soup_seurat <- run_seurat_pipeline(soup_seurat, "SoupSeurat");
saveRDS(soup_seurat, paste(OPTS$out_prefix, "SoupX.rds", sep="_"))
