my_metadata_table <- read.table("~/scripts/LiverMap2.0/LiverMap_SampleProcessingParams.csv", sep=",", header=T, stringsAsFactors=FALSE);
source("~/scripts/LiverMap2.0/Colour_Scheme.R")
source("~/scripts/LiverMap2.0/Setup_autoannotation.R")

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

tmp_file <- paste(name, "doublet_tmp.rds", sep="_")


print(name);

if(!file.exists(tmp_file)) {

if (file.exists(paste(name,"Anno_SeurObj.rds", sep="_"))) {
myseur <- readRDS(paste(name,"Anno_SeurObj.rds", sep="_"));
} else {
print("Anno Object file missing!")
}

clus_lab <- myseur@meta.data$seurat_clusters

type_lab <- cell_anno_to_cluster_anno(simplify_annotations(myseur@meta.data$general_labs, c("Mac", "T", "Hep")), myseur@meta.data$seurat_clusters)
type_lab <- type_lab[as.numeric(clus_lab),2]

# optimize parameters
sweep.res.list_liver <- paramSweep_v3(myseur, PCs = 1:npcs, sct = F)
sweep.stats_liver <- summarizeSweep(sweep.res.list_liver, GT = FALSE)
best <- which(sweep.stats_liver$BCreal == max(sweep.stats_liver$BCreal));

pN <- as.numeric(as.character((sweep.stats_liver[best,1])))

bcmvn_liver <- find.pK(sweep.stats_liver)
pK <- bcmvn_liver[bcmvn_liver$MeanBC == max(bcmvn_liver$MeanBC),2]
pK <- as.numeric(as.character(pK))

# estimate homotypics
n_cells <- ncol(myseur);
homotypic.prop <- modelHomotypic(type_lab) #Could also use clus_lab       
rate <- 0.9/100*n_cells/1000 #from:http://cgs.hku.hk/portal/files/GRC/Events/Seminars/2017/20170904/chromium_single_cell.pdf
nExp_poi <- round(rate*n_cells)  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#find doublets
pN <- min(250/n_cells, pN);
pN <- max(20/n_cells, pN);

#myseur <- doubletFinder_v3(myseur, PCs = 1:npcs, pN = pN, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = F) 
myseur <- doubletFinder_v3(myseur, PCs = 1:npcs, pN = pN, pK = pK, nExp = nExp_poi.adj, sct = F)

tmp <- ncol(myseur@meta.data)
myseur@meta.data$doublet_score <- myseur@meta.data[,tmp-1]
myseur@meta.data$DF_assignment <- myseur@meta.data[,tmp]

myseur@meta.data$cell_barcode <- colnames(myseur);
anno_tab <- myseur@meta.data

rownames(anno_tab) <- paste(anno_tab$orig.ident, anno_tab$cell_barcode, sep="_")
saveRDS(anno_tab, tmp_file)
all_doublet[[name]] <- anno_tab;
} else {
anno_tab <- readRDS(tmp_file)


all_doublet[[name]] <- anno_tab;

}
}

saveRDS(all_doublet, "All20_doubletDetection2.rds");





