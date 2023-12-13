require(dplyr)
require(Seurat)
require(Matrix)
require(DoubletFinder)

args <- commandArgs(trailingOnly=TRUE)
source("/cluster/home/tandrews/scripts/LiverMap2.0/Colour_Scheme.R")

myseur <- readRDS(args[1])
tmp <- unlist(strsplit(args[1], "/"))
tmp <- tmp[length(tmp)]
name <- sub(".rds", "", tmp)

print(name);

out_file <- paste(name, "doublets.rds", sep="_")

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
#pN <- min(250/n_cells, pN);
#pN <- max(20/n_cells, pN);

#myseur <- doubletFinder_v3(myseur, PCs = 1:npcs, pN = pN, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = F) 
myseur <- doubletFinder_v3(myseur, PCs = 1:npcs, pN = pN, pK = pK, nExp = nExp_poi.adj, sct = F)

tmp <- ncol(myseur@meta.data)
myseur@meta.data$doublet_score <- myseur@meta.data[,tmp-1]
myseur@meta.data$DF_assignment <- myseur@meta.data[,tmp]

myseur@meta.data$cell_barcode <- colnames(myseur);
anno_tab <- myseur@meta.data

rownames(anno_tab) <- paste(anno_tab$orig.ident, anno_tab$cell_barcode, sep="_")
saveRDS(anno_tab, tmp_file)
