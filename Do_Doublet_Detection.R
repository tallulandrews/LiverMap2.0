my_metadata_table <- read.table("Metadata20LiverMapPlusParams.csv", sep=",", header=T, stringsAsFactors=FALSE);

source("../AutoAnnotation/Setup_autoannotation.R")

all_doublet <- list();


for(dataset_row in 1:20) {

set.seed(my_metadata_table$Seed[dataset_row])
name <- my_metadata_table$Name[dataset_row]
npcs <- my_metadata_table$nPCs[dataset_row]
nkNN <- my_metadata_table$kNN[dataset_row]
res <- 5


print(name);
require(dplyr)
require(Seurat)
require(Matrix)
require(DoubletFinder)

if (file.exists(paste(name,"Anno_SeurObj.rds", sep="_"))) {
myseur <- readRDS(paste(name,"Anno_SeurObj.rds", sep="_"));

clus_lab <- myseur@meta.data$seurat_clusters
type_lab <- myseur@meta.data$marker_anno #fewer labels so better estimation!


sweep.res.list_liver <- paramSweep_v3(myseur, PCs = 1:10, sct = FALSE)
sweep.stats_liver <- summarizeSweep(sweep.res.list_liver, GT = FALSE)
bcmvn_liver <- find.pK(sweep.stats_liver)
homotypic.prop <- modelHomotypic(type_lab) #Could also use clus_lab       
rate <- 0.9/100*n_cells/1000 #from:http://cgs.hku.hk/portal/files/GRC/Events/Seminars/2017/20170904/chromium_single_cell.pdf
nExp_poi <- round(rate*n_cells)  
n_cells <- ncol(myseur);
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
myseur <- doubletFinder_v3(myseur, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
myseur <- doubletFinder_v3(myseur, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

anno_tab <- myseur@meta.data
anno_tab$cell_barcode <- colnames(myseur);

rownames(anno_tab) <- paste(anno_tab$orig.ident, anno_tab$cell_barcode, sep="_")

all_doublet[[name]] <- anno_tab;
}
}

saveRDS(all_doublet, "All20_automatedannotation.rds");





