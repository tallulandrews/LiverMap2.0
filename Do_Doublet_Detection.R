my_metadata_table <- read.table("~/scripts/LiverMap2.0/Metadata20LiverMapPlusParams.csv", sep=",", header=T, stringsAsFactors=FALSE);

require(dplyr)
require(Seurat)
require(Matrix)
require(DoubletFinder)

all_doublet <- list();


for(dataset_row in 1:nrow(my_metadata_table)) {

set.seed(my_metadata_table$Seed[dataset_row])
name <- my_metadata_table$Name[dataset_row]
npcs <- my_metadata_table$nPCs[dataset_row]
nkNN <- my_metadata_table$kNN[dataset_row]
res <- 5


print(name);

if (file.exists(paste(name,"Anno_SeurObj2.rds", sep="_"))) {
myseur <- readRDS(paste(name,"Anno_SeurObj2.rds", sep="_"));
} else {
print("Anno Object file missing!")
}
if (!file.exists(paste(name, "Anno_SeurObj.rds", sep="_"))) {next;}

clus_lab <- myseur@meta.data$seurat_clusters
type_lab <- myseur@meta.data$marker_anno #fewer labels so better estimation!

# optimize parameters
sweep.res.list_liver <- paramSweep_v3(myseur, PCs = 1:npcs, sct = TRUE)
sweep.stats_liver <- summarizeSweep(sweep.res.list_liver, GT = FALSE)
pN <- as.numeric(as.character((sweep.stats_liver[sweep.stats_liver[,3] == max(sweep.stats_liver[,3]),1])))
bcmvn_liver <- find.pK(sweep.stats_liver)
pK <- bcmvn_liver[bcmvn_liver$BCmetric == max(bcmvn_liver$BCmetric),2]
pK <- as.numeric(as.character(pK))

# estimate homotypics
n_cells <- ncol(myseur);
homotypic.prop <- modelHomotypic(type_lab) #Could also use clus_lab       
rate <- 0.9/100*n_cells/1000 #from:http://cgs.hku.hk/portal/files/GRC/Events/Seminars/2017/20170904/chromium_single_cell.pdf
nExp_poi <- round(rate*n_cells)  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#find doublets -- Buggy!!
myseur <- doubletFinder_v3(myseur, PCs = 1:npcs, pN = pN, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = T) 
myseur <- doubletFinder_v3(myseur, PCs = 1:npcs, pN = pN, pK = pK, nExp = nExp_poi.adj, reuse.pANN = paste("pANN", pN, pK, nExp_poi, sep="_"), sct = T)

anno_tab <- myseur@meta.data
anno_tab$cell_barcode <- colnames(myseur);

rownames(anno_tab) <- paste(anno_tab$orig.ident, anno_tab$cell_barcode, sep="_")

all_doublet[[name]] <- anno_tab;
}

saveRDS(all_doublet, "All20_doubletDetection.rds");





