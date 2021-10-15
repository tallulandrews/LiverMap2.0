args <- commandArgs(trailingOnly=TRUE)
#args = object, cluster_col, projname
obj <- readRDS(args[1])
cluster_col <- args[2]
proj_name <- args[3]

files <- Sys.glob("*Subcluste*.rds")

p_names <- c("AntiB_harmony", "Cholangiocyte_harmony", 
		"Endo_harmony", "Hepatocyte1_harmony",
		"Hepatocyte2_harmony", "Macrophage_harmony",
		"NKT_harmony", "Stellate_harmony")

for (i in 1:length(files)){

obj <- readRDS(files[i]); proj_name=p_names[i]

clusters <- obj@meta.data[,cluster_col]
confounders <- cbind(obj@meta.data$assay_type, obj@meta.data$sample)
#groups <- apply(confounders, 1, paste, collapse="_")

source("~/scripts/LiverMap2.0/My_R_Scripts.R")
pseudobulk <- get_pseudobulk_means(obj@assays$RNA@scale.data, 
		clusters, obj@meta.data$sample)
pseudobulk_sample <- unlist(lapply(strsplit(colnames(pseudobulk), "_"), function(x){x <- unlist(x)[-1]; paste(x, collapse="_")}))

metadata <- obj@meta.data[,c("sample", "donor", "assay_type", "donor_age", "donor_age_group", "donor_sex", "donor_bmi", "donor_bmi_group")]
metadata <- unique(metadata)
metadata <- metadata[match(pseudobulk_sample, metadata$sample),]

saveRDS(list(pseudobulk_mat=pseudobulk, metadata=metadata),paste(proj_name,  "subcluster_pseudobulks.rds", sep="_"))

}



