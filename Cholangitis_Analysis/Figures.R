## Cell-Type Colour Scheme
require("RColorBrewer")
require("ggplot2")
require(Seurat)
source("~/scripts/LiverMap2.0/Colour_Scheme.R")
source("~/scripts/LiverMap2.0/My_R_Scripts.R")

## ColourScheme ##
Healthy_cols = brewer.pal(9, "Blues")[5]
PSC_cols = brewer.pal(9, "Reds") [5]
PBC_cols = brewer.pal(9, "Purples") [5]

blues <- brewer.pal(5, "Blues")
greens <- brewer.pal(5, "Greens")
oranges <- brewer.pal(5, "Oranges")
reds <- brewer.pal(5, "Reds")


cell_type_cols <- list(
                "CD8+T"=Cell_type_colours[Cell_type_colours[,1] == "gdTcells2",2],
                "CD8T"=Cell_type_colours[Cell_type_colours[,1] == "gdTcells2",2],
                "CD4+T"=Cell_type_colours[Cell_type_colours[,1] == "Tcell",2],
                "CD4T"=Cell_type_colours[Cell_type_colours[,1] == "Tcell",2],
                "NKT"="purple",
                "CD8T-cNK"="pink",
                "CD4T-lrNK"="teal",
                "NKTcell"="purple",
                "cNK-CD8T"=greens[3],
                "lrNK-CD4T"="magenta",
                "CD4T-lrNK"="magenta",
                "CD4-lrNK"="magenta",
                "CD3T-lrNK"="magenta",
                "Treg"=greens[5],
                "MatB"="salmon",
                "MatB/CD4+"="salmon",
                "AntiB"=Cell_type_colours[Cell_type_colours[,1] == "AntiBcell",2],
                "Cellcycle"=Cell_type_colours[Cell_type_colours[,1] == "gdTcells1",2],
                "gdT"=Cell_type_colours[Cell_type_colours[,1] == "gdTcells1",2],
                "Tcell"=Cell_type_colours[Cell_type_colours[,1] == "gdTcells1",2],
                "lrNK"="#DB8E00",
                "cNK"=Cell_type_colours[Cell_type_colours[,1] == "NKcells",2],
                "miscNK"=Cell_type_colours[Cell_type_colours[,1] == "NKcells",2],
                "MAST"="green",
                "Kupffer"=Cell_type_colours[Cell_type_colours[,1] == "NonInfMac",2],
                "Monocyte"=Cell_type_colours[Cell_type_colours[,1] == "InfMac",2],
                "ActMac"="seagreen",
                "MHCII"="orchid",
                "cDC2"="blue",
                "cDC"="blue",
                "pDC"=blues[3],
                "LAM-like"=blues[5],
                "Neutrophil"="black",
                "RBC"=Cell_type_colours[Cell_type_colours[,1] == "Eryth",2],
                "MatB--RBC"=Cell_type_colours[Cell_type_colours[,1] == "Eryth",2],
                "Prolif"="navy",
                "Prolif-Mac"="navy",
                "Prolif-RBC"=Cell_type_colours[Cell_type_colours[,1] == "Eryth",2],
                "Hepato"=reds[3],
                "I-Hepato"=reds[2],
                "InterHep"=reds[2],
                "UnkHepato"=reds[3],
                "Unk-Hep"=reds[3],
                "Hepato--Mac"="purple",
                "cvLSEC--Mac"="purple",
                "C-Hepato"="red",
                "C-Hepato2"=Cell_type_colours[Cell_type_colours[,1] == "CentralHep",2],
                "P-Hepato"=Cell_type_colours[Cell_type_colours[,1] == "Eryth",2],
                "P-Hepato2"="navy",
                "Hep-Chol"="navy",
                "Stellate"=Cell_type_colours[Cell_type_colours[,1] == "Stellate",2],
                "a-Stellate"="lightgreen",
                "aStellate"="lightgreen",
                "q-Stellate"=oranges[5],
                "Stellate"=oranges[3],
                "Fibroblast"="grey50",
               "VSMC"="purple",
                "Doublet"="black",
                "Cholangiocyte"=Cell_type_colours[Cell_type_colours[,1] == "Cholangiocyte",2],
                "Chol"=Cell_type_colours[Cell_type_colours[,1] == "Cholangiocyte",2],
                "Chol-Mucosal"="goldenrod1",
                "CholMucus"="goldenrod1",
                "cvLSEC"=Cell_type_colours[Cell_type_colours[,1] == "cvLSECs",2],
                "ppLSEC"=Cell_type_colours[Cell_type_colours[,1] == "PortalLSECs",2],
                "cvEndo"=Cell_type_colours[Cell_type_colours[,1] == "Portalendo",2],
                "Arterial"="violet")


### Read Data ###
merged_obj = readRDS("SC_Integrated_Map_Annotated.rds")
outname = "Figures_SC"
prefix="SC_Integrated_Map"
#merged_obj = readRDS("SN_Integrated_Map_Annotated.rds")
#outname = "Figures_SN"
merged_obj@meta.data$Phenotype_plot <- as.character(merged_obj@meta.data$Phenotype)
merged_obj@meta.data$Phenotype_plot[merged_obj@meta.data$Phenotype_plot == "Healthy"] <- "NDD"


### Figure: UMAPs with PSC, PBC & Healthy overlay. ###

merged_obj <- merged_obj[,!grepl("Doublet", merged_obj@meta.data$cell_type)]
merged_obj@meta.data$cell_type <- factor(merged_obj@meta.data$cell_type, levels = rev(names(sort(table(merged_obj@meta.data$cell_type)))))

new_colour_scheme <- unlist(cell_type_cols)[match(levels(merged_obj@meta.data$cell_type), names(cell_type_cols))]
names(new_colour_scheme) <- levels(merged_obj@meta.data$cell_type)

png(paste(outname, "CellType_UMAP_with_labels.png", sep="_"), width=12, height=4, units="in", res=300)
DimPlot(merged_obj, split.by="Phenotype_plot", group.by="cell_type", label=TRUE)+scale_color_manual(values=new_colour_scheme)
dev.off()

png(paste(outname, "CellType_UMAP_no_labels.png", sep="_"), width=12, height=4, units="in", res=300)
DimPlot(merged_obj, split.by="Phenotype_plot", group.by="cell_type")+scale_color_manual(values=new_colour_scheme)
dev.off()

### Figure : Cell-type Frequency ###
freqs <-  table(merged_obj@meta.data$cell_type, merged_obj@meta.data$orig.ident)
freqs <- apply(freqs,2,function(x){x+1/sum(x)})

cell_type_freq_plot <- function(props) {
        biogroup <- rep("Healthy", ncol(props))
        biogroup[grepl("PSC", colnames(props))] <- "PSC"
        biogroup[grepl("PBC", colnames(props))] <- "PBC"

        source("~/scripts/LiverMap2.0/My_R_Scripts.R")
        to_plot <- props
        barplot_toplot <-group_rowmeans(to_plot, biogroup)

        tmp <- table(biogroup); exclude <- names(tmp[tmp == 1])
        if (length(exclude) > 0) {
                props <- props[,!biogroup %in% exclude]
                biogroup <- biogroup[!biogroup %in% exclude]
        }

        sig <- apply(props, 1, function(x){oneway.test(x~biogroup)$p.value})
        sig2 <- apply(props, 1, function(x){wilcox.test(x[biogroup=="PSC"],x[biogroup=="Healthy"])$p.value})
        comb_sig <- sig2
        comb_sig[rowSums(props == 0) > 1] <- sig[rowSums(props == 0) > 1]
        comb_sig[is.na(comb_sig)] <- sig[is.na(comb_sig)]
        sig_stars <- rep("", length(comb_sig))
        sig_stars[comb_sig < 0.05] <- "*"
        sig_stars[comb_sig < 0.01] <- "**"
        sig_stars[comb_sig < 0.001] <- "***"

        par(mar=c(4,12,1,1))
        loc <- barplot(t(barplot_toplot), col=c(Healthy_cols, PBC_cols, PSC_cols), beside=TRUE, horiz=T, las=2, xlab="frequency", ylab="", xlim=c(0, max(barplot_toplot)*1.2))
        text( apply(t(barplot_toplot),2,max), colMeans(loc), sig_stars, pos=4)
        legend("topright", fill=c(Healthy_cols, PBC_cols, PSC_cols), c("NDD", "PBC", "PSC"), bty="n")
}

props <- t(t(freqs)/colSums(freqs))*100
png(paste(outname, "type_byPhenotype_all.png", sep="_"), width=8, height=12, units="in", res=300)
cell_type_freq_plot(props)
dev.off()
png(paste(outname, "type_byPhenotype_major.png", sep="_"), width=8/1.5, height=12/2, units="in", res=300)
cell_type_freq_plot(props[rowMeans(props) > 1,])
dev.off()
png(paste(outname, "type_byPhenotype_minor.png", sep="_"), width=8/2, height=12/3, units="in", res=300)
cell_type_freq_plot(props[rowMeans(props) <= 1,])
dev.off()

### Figure nDE ###

all_disease_de <- list()
files <- Sys.glob(paste("SC_DE_files/",prefix,"*edgeR_pseudobulkDE_PBC-Healthy.csv", sep=""))
for (f in files) {
        tab <- read.table(f, sep=",")
        tmp <- unlist(strsplit(f, "_"))
        type <- tmp[length(tmp)-3]
        all_disease_de[[type]] <- tab;
}

# Number of DE genes vs cell-type
n_norm_down <- sapply(all_disease_de, function(x) {sum(x$logFC < 0)})
n_norm_up <- sapply(all_disease_de, function(x) {sum(x$logFC > 0)})

type_names <- names(all_disease_de);
my_xlim=c(max(c(n_norm_up,n_norm_down))*-1, max(c(n_norm_up,n_norm_down)))
png(paste(outname, "nDE_barplot.png", sep="_"), width=5, height=5, units="in", res=300)
par(mar=c(4.5,8.5,1,1))
xloc <- barplot(n_norm_down, names=type_names, xlim=my_xlim, horiz=TRUE, las=1, col="red")
barplot(-1*n_norm_up, names=type_names, xlim=my_xlim, horiz=TRUE, add=TRUE, las=1, col="blue",
                xlab="N DE genes in PSC vs Control")
dev.off()


### Figure DE: PSC vs PBC vs Healthy ###
files1 <- Sys.glob(paste(prefix, "*edgeR_pseudobulkDE_PSC-PBC.csv", sep=""))
all_disease_de_PSC_vs_Healthy <- list()
all_disease_de_PSC_vs_PBC <- list()
all_disease_de_PBC_vs_Healthy <- list()
for (f in files1) {
        tab <- read.table(f, sep=",")
        tmp <- unlist(strsplit(f, "_"))
        type <- tmp[length(tmp)-3]
        all_disease_de_PSC_vs_PBC[[type]] <- tab;
	file_tmp <- paste(prefix, type, "edgeR_pseudobulkDE_rerun.csv", sep="_")
	if (file.exists(file_tmp)) {
	        tab <- read.table(file_tmp, sep=",")
	        all_disease_de_PSC_vs_Healthy[[type]] <- tab;
	}
	file_tmp <- paste(prefix, type, "edgeR_pseudobulkDE_PBC-Healthy.csv", sep="_")
	if (file.exists(file_tmp)) {
        	tab <- read.table(file_tmp, sep=",")
        	all_disease_de_PBC_vs_Healthy[[type]] <- tab;
	}
}

core_types <- intersect(names(all_disease_de_PSC_vs_Healthy), intersect(names(all_disease_de_PSC_vs_PBC), names(all_disease_de_PBC_vs_Healthy)))
core_types <- core_types[!grepl("Doublet", core_types)]

total_DE <- list(PSC_vs_healthy=c(), PBC_vs_healthy=c(), PSC_and_PBC_vs_healthy=c())
venn_lists <- list(PSC_vs_healthy=c(), PBC_vs_healthy=c())

for(type in core_types) {
	psc_vs_health = all_disease_de_PSC_vs_Healthy[[type]]; 
	if (!is.factor(psc_vs_health[,1])) {
		psc_vs_health_id <- paste(sign(psc_vs_health[,1]), rownames(psc_vs_health))
	} else {
		# no DE genes
		psc_Vs_health_id <- c()
	}
	pbc_vs_health = all_disease_de_PBC_vs_Healthy[[type]]; 
	if (!is.factor(pbc_vs_health[,1])) {
		pbc_vs_health_id <- paste(sign(pbc_vs_health[,1]), rownames(pbc_vs_health))
	} else {
		# no DE genes
		pbc_Vs_health_id <- c()
	}
	common <- intersect(psc_vs_health_id, pbc_vs_health_id)
	total_DE$PSC_vs_healthy <- c(total_DE$PSC_vs_healthy, paste(psc_vs_health_id[!psc_vs_health_id %in% common], type))
	total_DE$PBC_vs_healthy <- c(total_DE$PBC_vs_healthy, paste(pbc_vs_health_id[!pbc_vs_health_id %in% common], type))
	total_DE$PSC_and_PBC_vs_healthy <- c(total_DE$PSC_and_PBC_vs_healthy, paste(common, type))
	venn_lists$PSC_vs_healthy <- c(venn_lists$PSC_vs_healthy, paste(psc_vs_health_id, type))
	venn_lists$PBC_vs_healthy <- c(venn_lists$PBC_vs_healthy, paste(pbc_vs_health_id, type))
}

require(ggvenn)

png(paste(outname, "PSC-Healthy_vs_PBC-Healthy_venn.png", sep="_"), width=6*0.6, height=4*0.6, units="in", res=300)
ggvenn::ggvenn(venn_lists, fill_color=c(PSC_cols, PBC_cols))
dev.off()

saveRDS(venn_lists, paste(outname, "PSC-Healthy_vs_PBC-Healthy_venn.rds", sep="_"))
#require(gplots)
#out <- venn(venn_lists)

### Figure Edge Signature ###
png(paste(outname, "Edge_signature_DotPlot.png", sep="_"), width=10*0.9, height=10, units="in", res=300)
DotPlot(merged_obj, features=c("SQSTM1", "IL32", "APOE", "MT1G", "HLA-A", "ANGPTL8", "AKR1B10", "CES1", "APOA2"), group.by="cell_type") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


png(paste(outname, "Edge_signature_DotPlot2.png", sep="_"), width=10*0.9, height=10, units="in", res=300)
DotPlot(merged_obj, features=c("SQSTM1", "HKDC1", "SULT1C2", "DAPK2", "AFF3", "PLA2G4C", 
				"KRT23", "AKR1B10", "CREB5", "CXCL8","BIRC3", "LAMA3", 
				"ENAH", "GALNT18", "RASEF", "THSD4", "FGF13", "DCDC2", "BICC1",
				"KRT7", "KRT19"), group.by="cell_type") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


png(paste(outname, "Edge_signature_DotPlot3.png", sep="_"), width=10*0.9*0.70, height=10*0.70, units="in", res=300)
DotPlot(merged_obj[,merged_obj@meta.data$Phenotype == "PSC"], features=c("AKR1B10", "LCN2", "SQSTM1", "TXN", "GPNMB", "COL4A1", "COL4A2", "COL1A1", "THY1", "TIMP1", "FASN", "MT1G", "PLA2G2A", "M1TH",  "IL32", "HLA-A"), group.by="cell_type") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

## Figure Cytokines ###
require("Seurat")
source("~/scripts/LiverMap2.0/My_R_Scripts.R")

table(merged_obj@meta.data$cell_type)/ncol(merged_obj)

major_celltypes <- c("P-Hepato", "P-Hepato2", "I-Hepato", "C-Hepato", "Chol", "Stellate", "Fibroblast", "cvLSEC", "ppLSEC", "MatB", "AntiB", "CD8T", "CD4T", "cNK", "lrNK", "CD8T-cNK", "CD3T-lrNK", "Tcell", "Monocyte", "Kupffer", "LAM-like", "ActMac", "Neutrophil"); # SC
major_celltypes <- c("P-Hepato", "P-Hepato2", "I-Hepato", "C-Hepato", "C-Hepato2","Chol", "Stellate", "Fibroblast", "cvLSEC", "ppLSEC", "AntiB", "CD4T", "lrNK", "Monocyte", "Kupffer"); # SN

phenogroup_expr <- get_pseudobulk(merged_obj@assays$RNA@data, merged_obj@meta.data$cell_type, merged_obj@meta.data$Phenotype, method="mean")
phenogroup_detect <- get_pseudobulk(merged_obj@assays$RNA@counts > 0, merged_obj@meta.data$cell_type, merged_obj@meta.data$Phenotype, method="mean")
phenogroup_pheno <- sapply(strsplit(colnames(phenogroup_detect), "_"), function(x){unlist(x[[2]])} )
phenogroup_celltype <- sapply(strsplit(colnames(phenogroup_detect), "_"), function(x){unlist(x[[1]])})

Cytokines <- c("CCL21", "CCL19", "CCR7", "ACKR4", "CXCL10", "CXCR3", "IFNG", "IGNGR1", "IFNGR2",
                "IL2", "IL2RA", "IL2RB", "IL2RG", "IL21", "IL21R",
                "CXCL16", "CXCR6", "CCL3", "CCL4", "CCL5", "CCR5", "CCR1", "CCR4", "CCR3",
                "CX3CL1", "CX3CR1", "IL6", "IL6R", "IL6ST", "CXCL8", "CXCR1", "CXCR2")

png(paste(outname, "cytokine_dotplot.png", sep="_"), width=8, height=12,units="in", res=150)
DotPlot(merged_obj[, merged_obj@meta.data$cell_type %in% major_celltypes], features=Cytokines, group.by="cell_type", split.by="Phenotype_plot", cols=c(Healthy_cols, PBC_cols, PSC_cols))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

### ----- SC additional immune celltypes ----- ###

label_genes <- function(x, this_name) { names(x) <- rep(this_name, length(x)); return(x)}

marker_genes <- read.csv("SupplementaryTable2MarkerGeneLists.csv")
lymphocytes_genes <- marker_genes[marker_genes[,3] == "Lymphocyte",]
major_celltypes <- c("P-Hepato", "P-Hepato2", "I-Hepato", "C-Hepato", "Chol", "Stellate", "Fibroblast", "cvLSEC", "ppLSEC", "MatB", "AntiB", "CD8T", "CD4T", "cNK", "lrNK", "CD8T-cNK", "CD3T-lrNK", "Tcell", "Monocyte", "Kupffer", "LAM-like", "ActMac", "Neutrophil"); # SC
lympho_celltypes <- c("CD8T", "CD4T", "cNK", "lrNK", "CD8T-cNK", "CD3T-lrNK", "Tcell"); # SC

png(paste(outname, "lymphocyte_dotplot.png", sep="_"), width=9, height=4,units="in", res=300)
DotPlot(merged_obj[, merged_obj@meta.data$cell_type %in% lympho_celltypes], features=c(as.character(lymphocytes_genes[,1]), "KLRB1", "KLRG1", "IL2RA", "CTLA4", "FOXP3", "TIGIT"), group.by="cell_type")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

doublets <- readRDS(paste(prefix, "Doublets.rds", sep="_"))
doublets <- doublets[match(merged_obj@meta.data$cell_ID, doublets$cell_ID),]

tmp <- table(merged_obj@meta.data$cell_type, doublets$DF_assignment)
tmp2 <- round(tmp/rowSums(tmp), digits=2)*100

png(paste(outname, "lymphocyte_doublets.png", sep="_"), width=4, height=3,units="in", res=300)
par(mar=c(6,4,1,1))
barplot(tmp2[rownames(tmp2) %in% lympho_celltypes,1], xlab="", ylab="Doublets (%)", las=2,
	col=new_colour_scheme[match(rownames(tmp2), names(new_colour_scheme))])
dev.off()


######### Cell Type Marker Genes ############
key_markers <- read.csv("/cluster/home/tandrews/scripts/PostReview_Disease_vs_Healthy_Map/Celltype_Marker_Genes_from_PSC.csv")

subset <- merged_obj[rownames(merged_obj) %in% key_markers[,1],]
Idents(subset) <- subset@meta.data$cell_type
res <- FindAllMarkers(subset, logfc.threshold=0, min.pct=0, only.pos=TRUE)
out <- c();
for (g in key_markers[,1]) {
        if (sum(res$gene == g) < 1) {next;}
        de <- res[res$gene == g  & res$cluster != "Hepatocyte2_Female" & res$cluster != "NKT_Flush",]
        l2fc <- de[,2]
        out <- rbind(out,de[which(l2fc == max(l2fc)),])
}
saveRDS(res, paste(outname,"DE_stats_for_Supplementary_Table_diseasemap.res", sep="_"))
write.table(out,paste(outname, "DE_stats_for_Supplementary_Table_diseasemap.csv", sep="_"), sep=",")


######### Cluster Tables ##########


tab <- table(merged_obj@meta.data$cell_type, merged_obj@meta.data$Phenotype)
write.table(tab,paste(outname, "cell_counts_per_phenotype.csv", sep="_"), sep=",")

tab <- table(merged_obj@meta.data$cell_type, merged_obj@meta.data$orig.ident)
write.table(tab,paste(outname, "cell_counts_per_sample.csv", sep="_"), sep=",")


### Macrophages ###
#
# THIS IS NOW A SEPARATE SCRIPT
#



all_disease_de <- list()
all_disease_de_towrite <- data.frame()
files <- Sys.glob(paste(outname,"*edgeR_pseudobulkDE_PSC-PBC.csv", sep=""))
for (f in files) {
        tab <- read.table(f, sep=",")
        tmp <- unlist(strsplit(f, "_"))
        type <- tmp[length(tmp)-3]
	if (tab[1,1] != "logFC") {
		tab$cell_type <- type
        	all_disease_de[[type]] <- tab;
		all_disease_de_towrite <- rbind(all_disease_de_towrite, tab)
	}
}

write.table(all_disease_de_towrite, file=paste(outname, "PSC-PBC_summary.csv", sep="_"), sep=",", row.names=T, col.names=T)



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

Mac_types <- c("Monocyte", "Kupffer", "LAM-like", "Neutrophil", "ActMac", "Hepato--Mac")
for (t in Mac_types) {
	all <- union(rownames(all_disease_de_PSC_vs_Healthy[[t]]), rownames(all_disease_de_PBC_vs_Healthy[[t]]))
	out <- data.frame(gene=all, 
			logFC_PSC=rep(0,length(all)), FDR_PSC=rep(1, length(all)), 
			logFC_PBC=rep(0, length(all)), FDR_PBC=rep(1, length(all)))
	psc <- all_disease_de_PSC_vs_Healthy[[t]]
	psc <- psc[match(all, rownames(psc)),]


	pbc <- all_disease_de_PBC_vs_Healthy[[t]]
	pbc <- pbc[match(all, rownames(pbc)),]
	out$logFC_PSC <- psc$logFC; out$FDR_PSC <- psc$FDR
	out$logFC_PBC <- pbc$logFC; out$FDR_PBC <- pbc$FDR
	write.table(out, paste(outname, t, "PSC_vs_PBC_vs_healthy_DE.csv", sep="_"), sep=",")
}

	
### NKT cell DE ###
cd8tNK_marks1 <- FindMarkers(merged_obj, ident.1="CD8T-cNK", ident.2="cNK")
cd8tNK_marks2 <- FindMarkers(merged_obj, ident.1="CD8T-cNK", ident.2="CD8T")
write.table(cd8tNK_marks1, file=paste(outname, "CD8T-cNK_vs_cNK.csv", sep="_"), sep=",")
write.table(cd8tNK_marks2, file=paste(outname, "CD8T-cNK_vs_CD8T.csv", sep="_"), sep=",")
comb <- intersect(rownames(cd8tNK_marks1), rownames(cd8tNK_marks2))
cd8tNK_marks1 <- cd8tNK_marks1[match(comb, rownames(cd8tNK_marks1)),]
cd8tNK_marks2 <- cd8tNK_marks2[match(comb, rownames(cd8tNK_marks2)),]
geom_mean <- cd8tNK_marks1[,2]*cd8tNK_marks2[,2]



cd3tNK_marks1 <- FindMarkers(merged_obj, ident.1="CD3T-lrNK", ident.2="lrNK")
cd3tNK_marks2 <- FindMarkers(merged_obj, ident.1="CD3T-lrNK", ident.2="Tcell")
write.table(cd3tNK_marks1, file=paste(outname, "CD3T-lrNK_vs_lrNK.csv", sep="_"), sep=",")
write.table(cd3tNK_marks2, file=paste(outname, "CD3T-lrNK_vs_Tcell.csv", sep="_"), sep=",")
comb <- intersect(rownames(cd3tNK_marks1), rownames(cd3tNK_marks2))
cd3tNK_marks1 <- cd3tNK_marks1[match(comb, rownames(cd3tNK_marks1)),]
cd3tNK_marks2 <- cd3tNK_marks2[match(comb, rownames(cd3tNK_marks2)),]
geom_mean <- cd3tNK_marks1[,2]*cd3tNK_marks2[,2]

