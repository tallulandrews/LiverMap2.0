#### ----- Set Up ----- ####

require("Seurat")
require("RColorBrewer")
require("ggplot2")
source("~/scripts/LiverMap2.0/Colour_Scheme.R")
source("~/scripts/LiverMap2.0/My_R_Scripts.R")

## ColourScheme ##
Healthy_cols = brewer.pal(9, "Blues")[5]
PSC_cols = brewer.pal(9, "Reds")[5]
PBC_cols = brewer.pal(9, "Purples")[5]

dir="/cluster/projects/macparland/TA/PostReview_Disease_vs_Healthy_Map"


##### Input #####
outname <- "SN_Integrated_Map";

print(outname)

## Filters: ##
n_features = 750
n_umi = 1000
scale_size = 2000



###### Transitioning Hepatocytes ######
### Cell-type Frequencies ###
merged_obj = readRDS(paste(outname, "_Annotated.rds", sep=""))

hepato_clusters <- c("P-Hepato2", "C-Hepato2", "P-Hepato", "C-Hepato", "I-Hepato")
chol_clusters <- c("CholMucus", "Chol")
stellate_clusters <- c("Stellate", "aStellate", "Fibroblast")

table(merged_obj@meta.data$cell_type, merged_obj@meta.data$Phenotype)


require(Seurat)
this_obj <- merged_obj[,merged_obj@meta.data$cell_type %in% hepato_clusters]
Idents(this_obj) <- this_obj@meta.data$cell_type
hepato_markers <- FindAllMarkers(this_obj)

this_obj <- merged_obj[,merged_obj@meta.data$cell_type %in% c(hepato_clusters, chol_clusters,stellate_clusters)]
Idents(this_obj) <- this_obj@meta.data$cell_type
parenchymal_markers <- FindAllMarkers(this_obj)

saveRDS(list(hepato=hepato_markers, paranchymal=parenchymal_markers), paste(outname, "parenchymal_markers.rds", sep="_"))

hepato_markers <- hepato_markers[hepato_markers$p_val_adj < 0.05,]
parenchymal_markers <- parenchymal_markers[parenchymal_markers$p_val_adj < 0.05,]
# Are any Hepato clusters more chol-like?
hep_and_chol_markers <- hepato_markers$gene %in% parenchymal_markers$gene[parenchymal_markers$cluster %in% chol_clusters]
n_chol_in_hep <- table(hepato_markers$cluster[hep_and_chol_markers])

tmp <- table(merged_obj@meta.data$cell_type, merged_obj@meta.data$Phenotype); tmp2 <- t(tmp)/colSums(tmp)
props_hep <- tmp2[,match(names(n_chol_in_hep),colnames(tmp2))]

png("Figure_Disease_Hep_Chol.png", width=6*2, height=6, units="in", res=300)
par(mfrow=c(1,2))
par(mar=c(4,6,1,1))
loc <- barplot(props_hep, col=c(Healthy_cols, PBC_cols, PSC_cols), beside=TRUE, horiz=T, las=2, xlab="frequency", ylab="")
text(apply(props_hep, 2, max)+0.025, loc[2,], c("","","","*","*"))
barplot(n_chol_in_hep, horiz=T, las=2, xlab="Cholangiocyte Markers (n)", ylab="")
dev.off()

# I-Hepato are disease-specific, are they cholangiocyte like?
chol_markers <- parenchymal_markers[parenchymal_markers$cluster %in% chol_clusters,]; chol_markers[order(chol_markers$p_val_adj, decreasing=FALSE),]

top_chol_in_hep <- c(head(parenchymal_markers[parenchymal_markers$cluster %in% chol_clusters & parenchymal_markers$gene %in% hepato_markers$gene,], 20)$gene)
top_hep_not_chol <- c(head(parenchymal_markers[parenchymal_markers$cluster %in% hepato_clusters & !parenchymal_markers$gene %in% chol_markers ,], 20)$gene)

merged_obj@meta.data$Phenotype2 <- as.character(merged_obj@meta.data$Phenotype)
merged_obj@meta.data$Phenotype2[merged_obj@meta.data$Phenotype2=="Healthy"] <- "NDD"
merged_obj@meta.data$Phenotype2 <- factor(merged_obj@meta.data$Phenotype2, levels=c("PSC", "PBC", "NDD"))
this_obj <- merged_obj[,merged_obj@meta.data$cell_type %in% c(hepato_clusters, chol_clusters)]
this_obj@meta.data$cell_type <- factor(this_obj@meta.data$cell_type)
Idents(this_obj) <- this_obj@meta.data$cell_type

require(ggplot2)
genes <- c("HAL", "SDS", "PDE7B", "CPS1","CYP3A4", "ABCC2", "GBE1", "TPRG1", "ALB", "APOA1","BICC1", "FGF13", "DCDC2", "ANKS1B",  "ANXA4", "ACSM3", "CTNND2", "RALYL", "CFTR", "EFNA5", "MACC1", "TMC5", "NRG3", "MUC5B", "KRT19", "KRT7")
Idents(this_obj) <- factor(Idents(this_obj), levels=c("P-Hepato", "P-Hepato2", "C-Hepato", "C-Hepato2","I-Hepato","Chol","CholMucus"))
png("SN_Chol_vs_Hep.png", width=10, height=5, units="in", res=100)
DotPlot(this_obj, features=genes, cols=c(PSC_cols, PBC_cols, Healthy_cols), split.by="Phenotype2")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
















## Disease DE ##

all_disease_de_PSC_vs_Healthy <- list()
files <- Sys.glob(paste(outname,"*edgeR_pseudobulkDE.csv", sep=""))
for (f in files) {
        tab <- read.table(f, sep=",")
        tmp <- unlist(strsplit(f, "_"))
        type <- tmp[length(tmp)-2]
        if (tab[1,1] != "logFC") {
                tab$cell_type <- type
                all_disease_de_PSC_vs_Healthy[[type]] <- tab;
        }
}

all_disease_de_PBC_vs_Healthy <- list()
files <- Sys.glob(paste(outname,"*edgeR_pseudobulkDE_PBC-Healthy.csv", sep=""))
for (f in files) {
        tab <- read.table(f, sep=",")
        tmp <- unlist(strsplit(f, "_"))
        type <- tmp[length(tmp)-3]
        if (tab[1,1] != "logFC") {
                tab$cell_type <- type
                all_disease_de_PBC_vs_Healthy[[type]] <- tab;
        }
}

all_disease_de_PSC_vs_PBC <- list()
files <- Sys.glob(paste(outname,"*edgeR_pseudobulkDE_PSC-PBC.csv", sep=""))
for (f in files) {
        tab <- read.table(f, sep=",")
        tmp <- unlist(strsplit(f, "_"))
        type <- tmp[length(tmp)-3]
        if (tab[1,1] != "logFC") {
                tab$cell_type <- type
                all_disease_de_PSC_vs_PBC[[type]] <- tab;
        }
}

# Which disease DE genes for hepatocytes are also cholangiocyte markers? #

all_disease_de_PSC_vs_PBC[[hepato_clusters]]


