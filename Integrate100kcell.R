require("Seurat")

seurfiles <- c("C37_SeurObj.rds", "C46_SeurObj.rds", "C52_SeurObj.rds", "C59_SeurObj.rds", "C39_SeurObj.rds", "C48_SeurObj.rds", "C53_SeurObj.rds", "C61_SeurObj.rds", "C41_SeurObj.rds", "C49_SeurObj.rds", "C54_SeurObj.rds", "C62_SeurObj.rds", "C42_SeurObj.rds", "C50_SeurObj.rds", "C56_SeurObj.rds", "C63_SeurObj.rds", "C43_SeurObj.rds", "C51_SeurObj.rds", "C58_SeurObj.rds", "C64_SeurObj.rds");

samp_names <- unlist(lapply(strsplit(seurfiles, "_"), function(x){x[[1]]}))
obj_list <- list()
common_genes <- c();
for (i in 1:length(seurfiles)) {
	n <- samp_names[i];
	obj <- readRDS(seurfiles[i]);
	if (length(common_genes) == 0) {
		common_genes <- rownames(obj);
	} else {
		common_genes <- intersect(common_genes, rownames(obj));
	}
	obj@meta.data$donor <- rep(n, ncol(obj));
	obj_list[[n]] <- obj
}

common_genes <- sort(common_genes);
common_genes <- common_genes[!grepl("^MT-", common_genes)]
for (i in 1:length(obj_list)) {
	obj_list[[i]] <- obj_list[[i]][match(common_genes, rownames(obj_list[[i]])),]
}


anchors <- FindIntegrationAnchors(object.list=obj_list, dims=1:20)

liver.integrated <- IntegrateData(achorset = anchors, dims=1:20)

saveRDS(liver.integrated, "basicintegration.rds")
