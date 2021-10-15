obj <- readRDS("Merged_EmptyOnly_obj_Map2.2_ImportedClusters.rds")
core_manual <- c("PortalHep", "CentralHep", "CD3abTcells", "OtherHep", "Macrophage", "LSECs", "PortalHep", "Stellate", "gdTcells", "AntiBcell", "Cholangiocyte", "ErythHep", "Eryth")
coarse_manual <- c("PortalHep", "PortalHep", "PortalHep", "CD3abTcells", "CentralHep", "CentralHep", "cvLSECs", "InfMac", "NonInfMac", "CentralHep", "PortalEndo", "gdTcells", "NKcells", "Stellate", "AntiBcell", "Cholangiocyte", "MatBcell", "ErythHep", "Eryth", "CentralHep")
fine_manual <- c("PortalHep", "OtherHep", "PortalHep", "CentralHep", "PortalHep", "CentralHep", "NKcells", "CD3abTcells", "cvLSECs", "CentralHep", "Macrophage", "CentralHep", "PortalEndo", "InfMac", "NonInfMac", "gdTcells", "PortalHep", "OtherHep", "MatBcell", "cvLSECs", "CentralHep", "Macrophage", "NKcells", "Stellate", "OtherHep", "AntiBcell", "Cholangiocyte", "PortalHep", "ErythHep", "Ambiguous", "CD3abTcells", "Ertyh", "PortalHep", "NonInfMac", "CentralHep", "PortalHep")
obj@meta.data$Core_Manual_Anno <- core_manual[obj@meta.data$Core_clusters]
obj@meta.data$Coarse_Manual_Anno <- coarse_manual[obj@meta.data$Coarse_clusters]
obj@meta.data$Fine_Manual_Anno <- fine_manual[obj@meta.data$Fine_clusters]
#saveRDS(obj, "Merged_EmptyOnly_obj_Map2.2_ImportedClusters_ManualAnno.rds")
source("~/scripts/LiverMap2.0/My_R_Scripts.R")
obj@meta.data$donor_age[obj@meta.data$donor == "C68"] <- 61
obj@meta.data$donor_age[obj@meta.data$donor == "C72"] <- 23
obj@meta.data$donor_sex[obj@meta.data$donor == "C68"] <- "F"
obj@meta.data$donor_sex[obj@meta.data$donor == "C72"] <- "F"
obj@meta.data$donor_bmi[obj@meta.data$donor == "C72"] <- 18.4
obj@meta.data$donor_age_group[obj@meta.data$donor == "C68"] <- "elderly"
obj@meta.data$donor_age_group[obj@meta.data$donor == "C72"] <- "young"
obj@meta.data$donor_bmi_group[obj@meta.data$donor == "C72"] <- "normal"
#saveRDS(obj, "Merged_EmptyOnly_obj_Map2.2_ImportedClusters_ManualAnno.rds")
pseudobulk_3pr <- get_pseudobulk(obj@assays$RNA@counts[,obj@meta.data$assay_type=="3pr"], 
		obj@meta.data$Coarse_Manual_Anno[obj@meta.data$assay_type=="3pr"], 
		obj@meta.data$donor[obj@meta.data$assay_type=="3pr"])
#saveRDS(pseudobulk_3pr, "Merged_EmptyOnly_obj_Map2.2_allgenes_3pr_pseudobulks.rds")
pseudobulk_5pr <- get_pseudobulk(obj@assays$RNA@counts[,obj@meta.data$assay_type=="5pr"], 
		obj@meta.data$Coarse_Manual_Anno[obj@meta.data$assay_type=="5pr"], 
		obj@meta.data$donor[obj@meta.data$assay_type=="5pr"])
#saveRDS(pseudobulk_5pr, "Merged_EmptyOnly_obj_Map2.2_allgenes_5pr_pseudobulks.rds")
source("~/scripts/LiverMap2.0/Map2.2_add_metadata.R")
donor <- sapply(strsplit(colnames(pseudobulk_3pr), "_"), function(x){x[[2]]})
type <- sapply(strsplit(colnames(pseudobulk_3pr), "_"), function(x){x[[1]]})
meta <- get_metadata(donor)
meta_3pr <- cbind(meta, donor, type)
donor <- sapply(strsplit(colnames(pseudobulk_5pr), "_"), function(x){x[[2]]})
type <- sapply(strsplit(colnames(pseudobulk_5pr), "_"), function(x){x[[1]]})
meta <- get_metadata(donor)
meta_5pr <- cbind(meta, donor, type)
saveRDS(list(psuedobulk=pseudobulk_3pr, metadata=meta_3pr), "Merged_EmptyOnly_obj_Map2.2_allgenes_3pr_pseudobulks.rds")
saveRDS(list(psuedobulk=pseudobulk_5pr, metadata=meta_5pr), "Merged_EmptyOnly_obj_Map2.2_allgenes_5pr_pseudobulks.rds")
