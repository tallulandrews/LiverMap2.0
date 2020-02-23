my_metadata_table <- read.table("Metadata20LiverMapPlusParams.csv", sep=",", header=T, stringsAsFactors=FALSE);

source("../AutoAnnotation/Setup_autoannotation.R")

all_anno <- list();


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

#scmap based
if (file.exists(paste(name,"Anno_SeurObj.rds", sep="_"))) {
myseur <- readRDS(paste(name,"Anno_SeurObj.rds", sep="_"));

orig_dims <- dim(myseur)

myseur <- run_scmap_seurat(myseur, scmap_ref=map1_ref);

new_dims <- dim(myseur)
print(name)
print(identical(orig_dims, new_dims))
# cluster based
norm <- myseur@assays$RNA@data
clus_lab <- myseur@meta.data$seurat_clusters

res <- Use_markers_for_anno(norm, clus_lab)
myseur@meta.data$marker_anno <- res$cell_assign;

saveRDS(myseur, paste(name,"Anno_SeurObj2.rds", sep="_"));

anno_tab <- myseur@meta.data
anno_tab$cell_barcode <- colnames(myseur);

rownames(anno_tab) <- paste(anno_tab$orig.ident, anno_tab$cell_barcode, sep="_")

all_anno[[name]] <- anno_tab;
}
}

saveRDS(all_anno, "All20_automatedannotation.rds");





