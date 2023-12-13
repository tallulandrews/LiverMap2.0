#### ----- Set Up ----- ####

require("Seurat")
require("ggplot2")
source("~/scripts/LiverMap2.0/Colour_Scheme.R")

dir="/cluster/projects/macparland/TA/PostReview_Disease_vs_Healthy_Map"

args <- commandArgs(trailingOnly=TRUE)

##### Input #####
if (args[1] == "SC") {
	outname <- "SC_Integrated_Map";
} else if (args[1] == "SN") {
	outname <- "SN_Integrated_Map";
} else {
	print("Error:Must provide either SC or SN as argument!")
	exit()
}

#### ----- Cytokines ----- ####


require("Seurat")
source("~/scripts/LiverMap2.0/My_R_Scripts.R")
obj <- readRDS(paste(outname, "_Annotated.rds", sep=""))

require("RColorBrewer")
Healthy_cols = brewer.pal(9, "Blues")
PSC_cols = brewer.pal(9, "Reds")
PBC_cols = brewer.pal(9, "Purples")

## Custom_thing - TLM ##
# Manual survey of cytokine interactions (from KEGG):
# CCL21, CCL19 -> CCR7, ACKR4
# CXCL10 -> CXCR3
# IFNG -> IFNGR1, IFNGR2
# IL2 -> IL2RA, IL2RB, IL2RG, IL2RB, IL2RG
# IL21R -> IL21
# CXCR6 -> CXCL16
# CCL3 -> CCR5, CCR1
# CCL4 -> CCR5
# CCL5 -> CCR3, CCR5, CCR4, CCR1
# CX3CR1 -> CX3CL1
# IL32 ->
# IL6 -> IL6R, IL6ST
# CXCL8 -> CXCR1, CXCR2

# Pathway IL2 -> NK/T -> IL32 -> Macs -> TNFalpha
genes <- c("IFNG", "IFNA1", "TNF", "CXCL8", "CXCL10","CCL21", "CCL19", "CCR7", "ACKR4" )

phenogroup_expr <- get_pseudobulk(obj@assays$RNA@data, obj@meta.data$cell_type, obj@meta.data$Phenotype, method="mean")
phenogroup_detect <- get_pseudobulk(obj@assays$RNA@counts > 0, obj@meta.data$cell_type, obj@meta.data$Phenotype, method="mean")
phenogroup_pheno <- sapply(strsplit("_",colnames(phenogroup_detect)))
phenogroup_celltype <- sapply(strsplit("_", colnames(phenogroup_detect)))

Cytokines <- c("CCL21", "CCL19", "CCR7", "ACKR4", "CXCL10", "CXCR3", "IFNG", "IGNGR1", "IFNGR2",
		"IL2", "IL2RA", "IL2RB", "IL2RG", "IL21", "IL21R", 
		"CXCL16", "CXCR6", "CCL3", "CCL4", "CCL5", "CCR5", "CCR1", "CCR4", "CCR3",
		"CX3CL1", "CX3CR1", "IL6", "IL6R", "IL6ST", "CXCL8", "CXCR1", "CXCR2")

png(paste(outname, "cytokine_dotplot.png", sep="_"), width=8, height=6,units="in", res=150)
DotPlot(obj, features=Cytokines, group.by="cell_type", split.by="Phenotype", cols=c(Healthy_cols[5], PBC_cols[5], PSC_cols[5]))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

these_types <- c("Chol", "Stellate", "ppLSEC", "cNK", "lrNK", "CD8T", "CD4T", "AntiB", "MatB","Monocyte","Kupffer", "LAM-like" )

psc_detect <- phenotgroup_detect[,phenogroup_pheno=="PSC"]
pbc_detect <- phenotgroup_detect[,phenogroup_pheno=="PBC"]
healthy_detect <- phenotgroup_detect[,phenogroup_pheno=="Healthy"]

psc_detect[c("CXCL8", "CXCR1", "CXCR2"), these_types]
healthy_detect[c("CXCL8", "CXCR1", "CXCR2"), these_types]
these_colors <- c(cell_type_cols[["Cholangiocyte"]], cell_type_cols[["Stellate"]], 
cell_type_cols[["ppLSEC"]], cell_type_cols[["cNK"]], 
cell_type_cols[["lrNK"]], cell_type_cols[["CD8+T"]], 
cell_type_cols[["CD4+T"]], cell_type_cols[["AntiB"]])


PSC_detect <- phenogroup_detect[,phenogroup_pheno == "PSC"];colnames(phenotgroup_detect) <- phenogroup_celltype[phenogroup_pheno == "PSC"]
PBC_detect <- phenogroup_detect[,phenogroup_pheno == "PBC"];colnames(phenotgroup_detect) <- phenogroup_celltype[phenogroup_pheno == "PBC"]
Healthy_detect <- phenogroup_detect[,phenogroup_pheno == "Healthy"];colnames(phenotgroup_detect) <- phenogroup_celltype[phenogroup_pheno == "Healthy"]

png(paste(outname, "_barplots_CXCL8.png", sep="_"), width=8*0.5,*3 height=6*2*0.5,units="in", res=150)
par(mfcol=c(2,3))
barplot(PSC_detect["CXCL8",these_types], names=these_types, col=these_colors, ylab="Detection Rate", main="PSC\nCXCL8")
barplot(colSums(PSC_detect[c("CXCR1", "CXCR2"),these_types], names=these_types, col=these_colors, ylab="Detection Rate", main="CXCR1/CXCR2")

barplot(Healthy_detect["CXCL8",these_types], names=these_types, col=these_colors, ylab="Detection Rate", main="Healthy\nCXCL8")
barplot(colSums(Healthy_detect[c("CXCR1", "CXCR2"),these_types], names=these_types, col=these_colors, ylab="Detection Rate", main="CXCR1/CXCR2")

barplot(PBC_detect["CXCL8",these_types], names=these_types, col=these_colors, ylab="Detection Rate", main="PBC\nCXCL8")
barplot(colSums(PBC_detect[c("CXCR1", "CXCR2"),these_types], names=these_types, col=these_colors, ylab="Detection Rate", main="CXCR1/CXCR2")
dev.off()

#NK_enrichments <- data.frame(term=c("NK-mediated cytotoxicity", "Antigen presentation", "Th1/T2 differentiation", "Chemokine signaling", "Endocytosis"), pvalue=c(1.947*10^-16, 1.072*10^-14, 2.334*10^-6, 2.536*10^-5, 5.798*10^-4)) # from gprofiler


