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
	PSC = c(
	"/cluster/projects/macparland/TA/PSC/Processed/PSC018_5pr_caudate_EmptyOnly.rds",
	"/cluster/projects/macparland/TA/PSC/Processed/PSC014X_5pr_caudate_EmptyOnly.rds",
	"/cluster/projects/macparland/TA/PSC/Processed/PSC019_Caudate_5pr_EmptyOnly.rds",
	"/cluster/projects/macparland/TA/PSC/Processed/PSC024_Caudate_5pr_EmptyOnly.rds",
	"/cluster/projects/macparland/TA/PSC/Processed/PSC012_SC_5pr_EmptyOnly.rds",
	"/cluster/projects/macparland/TA/PSC/Processed/PSC016_SC_5pr_EmptyOnly.rds",
	"/cluster/projects/macparland/TA/PSC/Processed/PSC033_SC_5prV3_EmptyOnly.rds",
	"/cluster/projects/macparland/TA/PSC/Processed/PSC005_SC_Frozen_5prV2_EmptyOnly.rds"
	)

	PBC = c("/cluster/projects/macparland/TA/PBC/Processed/PBC005_Frozen_5prV2_EmptyOnly.rds",
	"/cluster/projects/macparland/TA/PBC/Processed/PBC001_Frozen_5pr_V2_EmptyOnly.rds")

	Healthy = c("/cluster/projects/macparland/TA/LiverMap2.0/Cleaned/C58_5pr_EmptyOnly.rds",
	    "/cluster/projects/macparland/TA/LiverMap2.0/Cleaned/C59_5pr_EmptyOnly.rds",
	    "/cluster/projects/macparland/TA/LiverMap2.0/Cleaned/C61_5pr_EmptyOnly.rds",
	    "/cluster/projects/macparland/TA/LiverMap2.0/Cleaned/C63_5pr_reseq_EmptyOnly.rds",
	    "/cluster/projects/macparland/TA/LiverMap2.0/Cleaned/C64_5pr_EmptyOnly.rds",
	    "/cluster/projects/macparland/TA/LiverMap2.0/Cleaned/C70_5pr_reseq_EmptyOnly.rds")

	all_samples <- c(PSC, PBC, Healthy)


	outname <- "SC_Integrated_Map";

print(outname)

sample_names <- strsplit(all_samples, "/");
sample_names <- sapply(sample_names, function(x){x[length(x)]})
sample_names <- sub("_EmptyOnly.rds", "", sample_names)

## Filters: ##
n_features = 750
n_umi = 1000
scale_size = 2000

## in silico gating ##

merged_obj = readRDS(paste(outname, "_Annotated.rds", sep=""))
PTPRC_pos <- merged_obj@assays$RNA@data["PTPRC",] > 0.75
CD68_pos <-  merged_obj@assays$RNA@data["CD68",] > 0.75
double_pos <- PTPRC_pos & CD68_pos
NDD <- sum(double_pos[merged_obj@meta.data$Phenotype =="Healthy"])/sum(PTPRC_pos[merged_obj@meta.data$Phenotype =="Healthy"])
PSC <- sum(double_pos[merged_obj@meta.data$Phenotype =="PSC"])/sum(PTPRC_pos[merged_obj@meta.data$Phenotype =="PSC"])
PBC <- sum(double_pos[merged_obj@meta.data$Phenotype =="PBC"])/sum(PTPRC_pos[merged_obj@meta.data$Phenotype =="PBC"])
nPSC <- sum(merged_obj@meta.data$Phenotype =="PSC" & PTPRC_pos)
nPBC <- sum(merged_obj@meta.data$Phenotype =="PBC" & PTPRC_pos)
nNDD <- sum(merged_obj@meta.data$Phenotype =="Healthy" & PTPRC_pos)

n = 326
(NDD-PBC)/sqrt(NDD*(1-NDD)/nNDD + PBC*(1-PBC)/nPBC)
 pnorm((NDD-PBC)/sqrt(NDD*(1-NDD)/nNDD + PBC*(1-PBC)/nPBC), lower.tail=FALSE)

png("in_silico_gating.png", width=3, height=5, units="in", res=300)
bar_loc <- barplot(c(NDD, PBC, PSC), xlab="", ylab="CD45+CD68+ / CD45+", col=c(Healthy_cols, PBC_cols, PSC_cols), names=c("NDD", "PBC", "PSC"), ylim=c(0, max(c(NDD, PBC, PSC))+0.1))
lines(bar_loc[c(1,2),1], rep(max(NDD, PBC)+0.025,2))
text(mean(bar_loc[1:2,1]), max(NDD, PBC)+0.025, "**", pos=3)
lines(bar_loc[c(1,3),1], rep(max(NDD, PSC)+0.05,2))
text(mean(bar_loc[1:3,1]), max(NDD, PSC)+0.05, "***",  pos=3)
lines(bar_loc[c(2,3),1], rep(max(PBC, PSC)+0.05,2))
text(mean(bar_loc[2:3,1]), max(PBC, PSC)+0.05, "*", pos=3)
dev.off()

#stars: * = p < 0.05, ** = p < 10^-5, *** = p < 10^-10




#### ----- Subcluster ----- ####

subcluster <- function(obj, clusters=c(0), cluster_col="seurat_clusters", res=2) {
	set.seed(101)
	subsample_cells <- obj[,obj@meta.data[,cluster_col] %in% clusters]
	subsample_cells <- RunPCA(subsample_cells)
	subsample_cells <- FindVariableFeatures(subsample_cells, method="vst", nfeatures=1000)
	subsample_cells <- RunPCA(subsample_cells, features=VariableFeatures(subsample_cells), ndims=1:20)
	subsample_cells <- FindNeighbors(subsample_cells, dims=1:12)
	subsample_cells <- FindClusters(subsample_cells, resolution=res)
	subsample_cells <- RunUMAP(subsample_cells, dims=1:20)
	png(paste("SubSubcluster_", outname, "_", paste(clusters, collapse="_"), "_umap.png", sep=""), width=8, height=8, units="in", res=300)
	print(DimPlot(subsample_cells))
	dev.off();
	return(subsample_cells)
}

###### Case vs Control ######
### Cell-type Frequencies ###
merged_obj = readRDS(paste(outname, "_Annotated.rds", sep=""))
tmp <- table(merged_obj@meta.data$cell_type, merged_obj@meta.data$Phenotype)
tmp <- tmp[!grepl("Hepato", rownames(tmp)),]
tmp2 <- t(tmp)/colSums(tmp)
tmp2[, order(tmp2[3,])]

subclustered_macs <- subcluster(merged_obj, clusters=c("Monocyte", "Kupffer", 
			"LAM-like", "Neutrophil", "ActMac"), 
			cluster_col="cell_type")
subclustered_macs@meta.data$subclusters <- subclustered_macs@meta.data$seurat_clusters
tmp <- table(subclustered_macs@meta.data$seurat_clusters, subclustered_macs@meta.data$Phenotype)
tmp <- table(subclustered_macs@meta.data$seurat_clusters, subclustered_macs@meta.data$orig.ident)


png(paste(outname, "MacSubcluster.png", sep="_"), width=12, height=4, units="in", res=300)
DimPlot(subclustered_macs, split.by="Phenotype", group.by="subclusters", label=TRUE)
dev.off()

png(paste(outname, "MacSubcluster_globaltypes.png", sep="_"), width=12, height=4, units="in", res=300)
DimPlot(subclustered_macs, split.by="Phenotype", group.by="cell_type", label=TRUE)
dev.off()

library(pals)
new_color_scheme = pals::alphabet2(n=length(unique(subclustered_macs@meta.data$orig.ident)))
names(new_color_scheme) <- levels(factor(subclustered_macs@meta.data$orig.ident))
png(paste(outname, "MacSubcluster_sample.png", sep="_"), width=12, height=4, units="in", res=300)
DimPlot(subclustered_macs, group.by="orig.ident", label=FALSE, split.by="Phenotype")  +
        scale_color_manual(values=new_color_scheme)
dev.off()

saveRDS(subclustered_macs, paste(outname, "MacSubcluster.rds", sep="_"))


subclustered_macs_markers <- FindAllMarkers(subclustered_macs, group.by="seurat_clusters")
write.table(subclustered_macs_markers, paste(outname, "MacSubcluster_markers.csv", sep="_"), sep=",")



tmp <-table(subclustered_macs@meta.data$seurat_clusters, subclustered_macs@meta.data$Phenotype)
tmp2 <- t(tmp)/colSums(tmp)

# Type specific: 
tmp2[3,] > 10*tmp2[2,] & tmp2[3,] > 10*tmp2[1,]


### Figures ###
subclustered_macs <- readRDS(paste(outname, "MacSubcluster.rds", sep="_"))

# Barplots

tmp <-table(subclustered_macs@meta.data$seurat_clusters, subclustered_macs@meta.data$Phenotype)
pheno_freq <- tmp/rowSums(tmp)
tmp <-table(subclustered_macs@meta.data$seurat_clusters, subclustered_macs@meta.data$orig.ident)
sample_freq <- tmp/rowSums(tmp)
colnames(sample_freq) <- sapply(strsplit(colnames(sample_freq), "_"), function(x){unlist(x[[1]])})

png("Figure_SC_MacSubcluster_freqBarplots.png", width=8*0.8, height=6*0.8, units="in", res=300)
layout(rbind(c(1,3), c(2,3)), widths=c(4,1), heights=c(1,1))
par(mar=c(4,4,1,1))
barplot(t(pheno_freq), col=c(Healthy_cols, PBC_cols, PSC_cols), xlab="", ylab="Phenotype", las=2)

sample_cols <- pals::alphabet2(ncol(sample_freq))
par(mar=c(4,4,1,1))
barplot(t(sample_freq), col=sample_cols, xlab="", ylab="Phenotype", las=2)
par(mar=c(0,0,0,0))
plot.new()
legend("topleft", fill=c(Healthy_cols, PBC_cols, PSC_cols), c("Healthy", "PBC", "PSC"), bty="n")
#par(mar=c(0,0,0,0))
#plot.new()
legend("bottomleft", fill=sample_cols, colnames(sample_freq), bty="n")
dev.off()


### Marker Plot ###
source("~/My_R_Packages/CellTypeProfiles/R/General.R")
source("~/My_R_Packages/CellTypeProfiles/R/Markers.R")

markers <- complex_markers(subclustered_macs@assays$RNA@data, subclustered_macs@meta.data$seurat_clusters)

markers2 <- read.csv(paste(outname, "MacSubcluster_markers.csv", sep="_"), sep=",")

table(factor(markers2[markers2$avg_logFC > 1,"gene"]))
toplot <- c("C1QA", "C1QB", "C1QC", "CD5L", "MARCO", "VCAM1", "VSIG4",
		"HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", 
		"LYZ", "S100A8", "S100A9", "VCAN", "LGMN", "GPNMB", "TREM2", "FCGR3B",
		"IL1B", "IL1R2", "IL18", "IL18BP", "IL10RA", "IL13RA1", "AREG", 
		"CCL18", "CCL4", "CCL20", "CCL3", "CXCL2", "CXCL3", "CXCL8", "CXCL10", "CXCL16", "BCL2L1", 
		"TNF", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TGFB1", "TGFBI", "NFKB1", "NFKBIA", "NFKBIZ")
		#"IFNG", "IFNK", "IFNLR1", "IFNGR1", "IFNL1", "IFNAR2", "IFNGR2")

cluster_theme <- c("PSC", "PSC", "Dis", "Dis", "Dis", 
			"Dis","Dis","Dis","Dis","PSC",
			"PSC", "PSC", "Dis","Dis","Dis",
			"Dis", "PSC", "Healthy", "PSC", "Dis",
			"PSC", "Dis", "Healthy", "Healthy", "Healthy",
			"Healthy", "Healthy", "Healthy", "PSC")
Idents(subclustered_macs) <- factor(paste(cluster_theme[as.numeric(subclustered_macs@meta.data$seurat_clusters)], subclustered_macs@meta.data$seurat_clusters))

png(paste(outname, "MacSubcluster_DotPlot.png", sep="_"), width=12, height=7, units="in", res=300)
DotPlot(subclustered_macs, features=toplot) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#Celltypes
cluster2type <- c("0"="Monocyte", "1"="Monocyte", "2"="Mac", "3"="HLA", "4"="Kupffer",
			"5"="LAM-like", "6"="Kupffer", "7"="Mac", "8"="Mac","9"="ActMono"
			"10"=, "11"="ActKupffer", "12"="Mac", "13"="Kupffer", "14"="Kupffer",
			"15"="Mac","16"="Mono","17"="Mono", "18"="Neutrophil", "19"="ActMono", 
			"20"="Mac", "21"="Mono", "22"="LAM-like", "23"="Kupffer", "24"="Mac",
			"25"="Kupffer", "26"="Mac", "27"="HLA", "28"="BCL2L1"
		)

## Pathway Plots ##

source("~/scripts/LiverMap2.0/Colour_Scheme.R")
source("~/scripts/LiverMap2.0/My_R_Scripts.R")

require(fgsea)
require(gProfileR)

immune_path <- gmtPathways("/cluster/projects/macparland/TA/ExternalData/MSigDb_immune_signatures_c7.all.v7.1.symbols.gmt")
Hallmark_path <- gmtPathways("/cluster/projects/macparland/TA/ExternalData/MSigDbHalmarkPathways.gmt")
MSigAll <- gmtPathways("/cluster/projects/macparland/TA/ExternalData/MSigDb_curated_c2.all.v7.1.symbols.gmt")
reactome <- gmtPathways("/cluster/projects/macparland/TA/ExternalData/ReactomePathways.gmt")
BaderMSig <- gmtPathways("/cluster/projects/macparland/TA/ExternalData/BaderLab25Aug2020/Human_MSigdb_August_01_2020_symbol.gmt.txt")
BaderReact <- gmtPathways("/cluster/projects/macparland/TA/ExternalData/BaderLab25Aug2020/Human_Reactome_August_01_2020_symbol.gmt.txt")



all_fsea_out <- list()
for (cluster in unique(markers2$cluster)) {
        set.seed(101);
	scores <- markers2$avg_logFC[markers2$cluster==cluster]
	names(scores) <- markers2$gene[markers2$cluster==cluster]
        paths <- do_fgsea(sort(scores), pathways=reactome)
        if (!is.null(paths)) {
                all_fsea_out[[as.character(cluster)]] <- paths
        }
}

all_fsea_out2 <- list()
for (cluster in unique(markers2$cluster)[22:28]) {
        set.seed(101);
	scores <- markers2$avg_logFC[markers2$cluster==cluster]
	names(scores) <- markers2$gene[markers2$cluster==cluster]
        paths <- do_fgsea(sort(scores), pathways=immune_path)
        if (!is.null(paths)) {
                all_fsea_out2[[as.character(cluster)]] <- paths
        }
}

saveRDS(list(all_fsea_out=all_fsea_out, all_fsea_out2=all_fsea_out2), file=paste(outname, "MacSubcluster_enrichments.rds", sep="_"))

richments <- readRDS(paste(outname, "MacSubcluster_enrichments.rds", sep="_"))
all_fsea_out=richments$all_fsea_out;
all_fsea_out2=richments$all_fsea_out2;

cluster_cols <- c(PSC_cols, PSC_cols, PBC_cols, PBC_cols, PBC_cols, 
			PBC_cols, PBC_cols, PBC_cols, PBC_cols, PSC_cols,
			PSC_cols, PSC_cols, PBC_cols, PBC_cols, PBC_cols,
			PBC_cols, PSC_cols, Healthy_cols, PSC_cols, PBC_cols,
			PSC_cols, PBC_cols, Healthy_cols, Healthy_cols, Healthy_cols,
			Healthy_cols, Healthy_cols, Healthy_cols, PSC_cols)
names(cluster_cols) <- 1:29 -1

all_paths <- sapply(all_fsea_out, function(x){unlist(as.character(x$rich$pathway))})
all_paths2 <- sapply(all_fsea_out[cluster_theme == "PSC"], function(x){unlist(as.character(x$rich$pathway))})
all_paths <- sapply(all_fsea_out2, function(x){unlist(as.character(x$rich$pathway))})
all_paths2 <- sapply(all_fsea_out2[cluster_theme[-21]=="PSC"], function(x){unlist(as.character(x$rich$pathway))})

make_path_plot <- function(this_path, fsea_out) {
	NES=c()
	for (cluster in names(fsea_out)) {
	        if (sum(fsea_out[[cluster]]$rich$pathway == this_path)) {
	                path_res = fsea_out[[cluster]]$rich
	                NES = c(NES, unlist(path_res[which(path_res$pathway == this_path)[1],"NES"]))
	        } else {
	                NES = c(NES, 0)
	        }
	}
	names(NES) = names(fsea_out)
	barplot(NES, main=this_path, col=cluster_cols[names(cluster_cols) %in% names(NES)], ylab="NES")
}
	

png(paste(outname, "pathway_barplots1.png", sep="_"), width=8*0.75, height=5*3*0.75, units="in", res=300)
par(mfrow=c(3,1))
make_path_plot("Interferon Signaling", all_fsea_out)
make_path_plot("MHC class II antigen presentation", all_fsea_out)
make_path_plot("Cellular responses to stress", all_fsea_out)
dev.off()


png(paste(outname, "pathway_barplots2.png", sep="_"), width=8*0.75, height=5*2*0.75, units="in", res=300)
par(mfrow=c(2,1))
make_path_plot("GSE9988_LPS_VS_CTRL_TREATED_MONOCYTE_UP", all_fsea_out2)
make_path_plot("GSE7509_DC_VS_MONOCYTE_UP", all_fsea_out2)
dev.off()


pdf(paste(outname, "pathway_barplots.pdf", sep="_"), width=8, height=5)
this_path = "Interferon Signaling"
this_path = "MHC class II antigen presentation"
this_path = "Cellular responses to stress"
NES=c()
for (cluster in names(all_fsea_out)) {
	if (sum(all_fsea_out[[cluster]]$rich$pathway == this_path)) {
		path_res = all_fsea_out[[cluster]]$rich
		NES = c(NES, unlist(path_res[which(path_res$pathway == this_path)[1],"NES"]))
	} else {
		NES = c(NES, 0)
	}		
}
names(NES) = names(all_fsea_out)
barplot(NES, main=this_path, col=cluster_cols[names(cluster_cols) %in% names(NES)])


this_path = "GSE9988_LPS_VS_CTRL_TREATED_MONOCYTE_UP"
this_path = "GSE3720_UNSTIM_VS_PMA_STIM_VD2_GAMMADELTA_TCELL_UP"
this_path = "GSE7509_DC_VS_MONOCYTE_UP"
NES=c()
for (cluster in names(all_fsea_out2)) {
        if (sum(all_fsea_out2[[cluster]]$rich$pathway == this_path)) {
                path_res = all_fsea_out2[[cluster]]$rich
                NES = c(NES, unlist(path_res[which(path_res$pathway == this_path)[1],"NES"]))
        } else {
                NES = c(NES, 0)
        }
}
names(NES) = names(all_fsea_out2)
barplot(NES, main=this_path, col=cluster_cols[names(cluster_cols) %in% names(NES)])

dev.off()


paths <- c()
for (c in c("0", "1", "9", "10", "11", "16", "18", "20", "28")) {
	paths <- c(paths, all_fsea_out2[[c]]$rich$pathway)
}
paths_pbc <- c()
for (c in c("2","3","4","5","6","8","12","13","14","15","21")) {
	paths_pbc <- c(paths_pbc, all_fsea_out2[[c]]$rich$pathway)
}

### Pathway Heatmap ###

pathways <- c("GSE41087_WT_VS_FOXP3_MUT_ANTI_CD3_CD28_STIM_CD4_TCELL_UP", "GSE3720_UNSTIM_VS_LPS_STIM_VD2_GAMMADELTA_TCELL_UP","GSE24634_IL4_VS_CTRL_TREATED_NAIVE_CD4_TCELL_DAY5_UP", "GSE16385_ROSIGLITAZONE_IFNG_TNF_VS_IL4_STIM_MACROPHAGE_UP", "GSE9988_LPS_VS_CTRL_TREATED_MONOCYTE_UP", "GSE36888_UNTREATED_VS_IL2_TREATED_TCELL_17H_UP", "GSE1925_CTRL_VS_IFNG_PRIMED_MACROPHAGE_3H_IFNG_STIM_UP", "GSE2405_0H_VS_24H_A_PHAGOCYTOPHILUM_STIM_NEUTROPHIL_UP", "GSE46606_UNSTIM_VS_CD40L_IL2_IL5_1DAY_STIMULATED_IRF4HIGH_SORTED_BCELL_DN","Interferon Signaling", "MHC class II antigen presentation", "Cellular responses to stress", "Cytokine Signaling in Immune system", "Neutrophil degranulation", "Signaling by Interleukins", "Complement cascade", "TCR signaling")
names(pathways) <- c("KO-FOXP3 CD4 TCells", "Stim gdTCells", "IL4 CD4TCells", "IFNG+TNF Mac", "LPS Monocytes", "IL2 TCells", "IFNG Mac", "Stim Neutrophil", "CD40L BCells dn", "Interferon Signaling", "MHCII presentation", "Cellular stress", "Cytokine Signaling", "Neutrophil degraulation", "Interleukin Signaling", "Complement cascade", "TCR signaling")
path_anno=data.frame(Type=c(rep("Immune",9), rep("KEGG",8)))
rownames(path_anno) <- names(pathways)

require("RColorBrewer")
markers2 <- read.csv(paste(outname, "MacSubcluster_markers.csv", sep="_"), sep=",")
Idents(subclustered_macs) <- factor(paste(cluster_theme[as.numeric(subclustered_macs@meta.data$seurat_clusters)], subclustered_macs@meta.data$seurat_clusters))
richments <- readRDS(paste(outname, "MacSubcluster_enrichments.rds", sep="_"))
all_fsea_out=richments$all_fsea_out;
all_fsea_out2=richments$all_fsea_out2;

unique_clusters <- sort(unique(Idents(subclustered_macs)))
unique_clusters_pheno <- data.frame(Phenotype=sapply(strsplit(as.character(unique_clusters), " "), function(x){x[[1]]}));
rownames(unique_clusters_pheno) <- unique_clusters
unique_cluster_num <- sapply(strsplit(as.character(unique_clusters), " "), function(x){x[[2]]})

out_mat <- matrix(NA, nrow=length(unique_clusters), ncol=length(pathways))
for (i in 1:length(unique_clusters)) {
        this_type <- as.character(unique_cluster_num[i])
        rich1 <- all_fsea_out[[this_type]]$rich
        rich2 <- all_fsea_out2[[this_type]]$rich
	rich <- rbind(rich1, rich2)
        out_mat[i,] <- unlist(rich[match(pathways, rich$pathway),"NES"])
}

rownames(out_mat) <- unique_clusters
colnames(out_mat) <- names(pathways)

my_heatmap_colours <- colorRampPalette(c(brewer.pal(n=9, "RdPu")[6:9], "black", rev(brewer.pal(n = 7, name = "RdYlBu")[1:4])))(200)

out_mat[is.na(out_mat)] <- 0
require(pheatmap)
legend_thing <- seq(from=floor(min(out_mat)), to=floor(max(out_mat)), length=5)
breaks = seq(from=-1*(max(abs(out_mat))), to=max(abs(out_mat)), length=length(my_heatmap_colours)+1)
png(paste(outname, "MacSubcluster_pathway_enrichment_heatmap.png", sep="_"), width=8, height=6*0.8, units="in", res=300)
pheatmap(t(out_mat), legend_breaks = c(legend_thing, max(out_mat)), main = "", legend_labels = c(legend_thing, "NES\n"), scale="none",
                col=my_heatmap_colours, breaks=breaks, annotation_col = unique_clusters_pheno, 
		annotation_colors=list(Phenotype=c(Dis=PBC_cols, Healthy=Healthy_cols, PSC=PSC_cols)),
		annotation_row=path_anno)
dev.off()

