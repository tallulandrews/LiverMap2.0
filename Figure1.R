
require(Seurat)
source("/cluster/home/tandrews/scripts/LiverMap2.0/Colour_Scheme.R")
source("/cluster/home/tandrews/scripts/LiverMap2.0/SubColour_Scheme.R")
colours_tech <- c("darkviolet", "orchid1")
source("/cluster/home/tandrews/scripts/LiverMap2.0/My_R_Scripts.R")

dir = "/cluster/projects/macparland/TA/LiverMap2.0/RawData/Analysis/EmptyDrops/Figures"

metaData <- readRDS("../Merged_EmptyOnly_obj_Map2.2_ImportedClusters_ManualAnno_MetaDataOnly.rds")

cluster_2_type <- table(metaData$Coarse_clusters, metaData$Coarse_Manual_Anno)
cluster_2_type <- colnames(cluster_2_type)[ apply(cluster_2_type, 1, function(x){which(x == max(x))}) ]
cluster_2_colour <- type_2_colour(cluster_2_type)

#mergedobj <- readRDS("../Map2.2_merged_integrated_harmony_plus_clusters.rds")
mergedobj <- readRDS("../Map_Empty_Integrated_harmony_plus_analysis.rds")
#mergedobj <- readRDS("../Merged_EmptyOnly_obj_Map2.2_ImportedClusters.rds")
#mergedobj <- readRDS("../Merged_integrated_with_metadata.rds")

# Central UMAP Plot
 require(ggplot2)
agg_coord_by_cluster <- function(coords, clusters) {
        x <- split(seq(nrow(coords)), clusters)
        result <- sapply(x, function(a) apply(coords[a,],2,median))
        return(result)
}

umap_lab_pos <- agg_coord_by_cluster(Reductions(mergedobj, "umap")@cell.embeddings, mergedobj@meta.data[,"Coarse_clusters"])
# Type colourscheme
new_colour_scheme <- Cell_type_colours[order(Cell_type_colours[,1]),]
mergedobj@meta.data[,"Cell Type"] <- map_cell_types(metaData[,"Coarse_Manual_Anno"])
new_colour_scheme <- new_colour_scheme[new_colour_scheme[,1] %in% mergedobj@meta.data[,"Cell Type"],]

# Change Macrophage annotations
new_colour_scheme[12,1] <- "Kupffer"
new_colour_scheme[9,1] <- "Monocyte-like"
mergedobj@meta.data[mergedobj@meta.data[,"Cell Type"]=="NonInfMac","Cell Type"] <- "Kupffer"
mergedobj@meta.data[mergedobj@meta.data[,"Cell Type"]=="InfMac","Cell Type"] <- "Monocyte-like"
metaData[metaData[,"Coarse_Manual_Anno"]=="NonInfMac","Coarse_Manual_Anno"] <- "Kupffer"
metaData[metaData[,"Coarse_Manual_Anno"]=="InfMac","Coarse_Manual_Anno"] <- "Monocyte-like"
new_colour_scheme <- new_colour_scheme[order(new_colour_scheme[,1]),]

#pdf("Figure1_UMAP_Cluster_and_Type.pdf", width=8, height=8)
png("Figure1_UMAP_Cluster_and_Type_relabelled.png", width=6, height=6, units="in", res=300)
print( DimPlot(mergedobj, reduction="umap", group.by="Cell Type", pt.size=.1)+scale_color_manual(values=new_colour_scheme[,2])+annotate("text", x=umap_lab_pos[1,], y=umap_lab_pos[2,], label=colnames(umap_lab_pos), colour="black") )
dev.off()
png("TLM_Figure1_UMAP_Cluster_and_Type_noLab.png", width=8, height=6, units="in", res=300)
print( DimPlot(mergedobj, reduction="umap", group.by="Cell Type", pt.size=.1)+scale_color_manual(values=new_colour_scheme[,2]))
dev.off()


# Aside: Flush vs Map - Supplementary #
is.Flush <- grepl("Flush", mergedobj@meta.data$sample)
tab <- table(is.Flush, mergedobj@meta.data[,"Cell Type"])
Prop_Flush <- tab[2,]/sum(tab[2,])
Prop_Caudate <- tab[1,]/sum(tab[1,])
diff <- Prop_Flush - Prop_Caudate
ordering <- order(Prop_Caudate, decreasing=T)
png("Suppl_Fig_Flush_vs_Caudate_WholeMap.png", width=6, height=4, units="in", res=300)
par(mar=c(10,4,1,1))
barplot(log2(Prop_Flush/Prop_Caudate)[ordering], col=new_colour_scheme[ordering,2], las=2, ylab="Flush vs Caudate (log2FC)")
dev.off()

# UMAP by demograph data

mergedobj@meta.data$donor_sex <- factor(metaData$donor_sex)

#pdf("Figure1_UMAP_Sex.pdf", width=8, height=8)
png("Figure1_UMAP_Sex.png", width=5, height=5, units="in", res=300)
print( DimPlot(mergedobj, reduction="umap", group.by="donor_sex", pt.size=.1)+scale_color_manual(values=colours_sex) )
dev.off()

#pdf("Figure1_UMAP_Age.pdf", width=8, height=8)
png("Figure1_UMAP_Age.png", width=5, height=5, units="in", res=300)
print( DimPlot(mergedobj, reduction="umap", group.by="donor_age_group", pt.size=.1)+scale_color_manual(values=colours_age) )
dev.off()

#pdf("Figure1_UMAP_Tech.pdf", width=8, height=8)
png("Figure1_UMAP_Tech.png", width=5, height=5, units="in", res=300)
print( DimPlot(mergedobj, reduction="umap", group.by="assay_type", pt.size=.1)+scale_color_manual(values=colours_tech) )
dev.off()

# Barplot by demographic data

bar_dat <- table(metaData$Coarse_Manual_Anno, factor(metaData$donor_sex, levels=c("F", "M")))
bar_dat <- bar_dat/rowSums(bar_dat)
#pdf("Figure1_Barplot_Sex.pdf", width=6, height=6)
png("Figure1_Barplot_Sex_renamed.png", width=4, height=4, units="in", res=300)
par(mar=c(8, 4, 1,1))
barplot(t(bar_dat)*100, col=colours_sex, las=2, main="", xlab="", ylab="Proportion (%)")
dev.off()

#bar_dat <- table(metaData$Coarse_Manual_Anno, mergedobj@meta.data$assay_type)
bar_dat <- table(metaData$Coarse_Manual_Anno, metaData$assay_type)
bar_dat <- bar_dat/rowSums(bar_dat)
#pdf("Figure1_Barplot_Tech.pdf", width=6, height=6)
png("Figure1_Barplot_Tech_renamed.png", width=4, height=4, units="in", res=300)
par(mar=c(8, 4, 1,1))
barplot(t(bar_dat)*100, col=colours_tech, las=2, main="", xlab="", ylab="Proportion (%)")
dev.off()

#bar_dat <- table(metaData$Coarse_Manual_Anno, mergedobj@meta.data$donor_age_group)
bar_dat <- table(metaData$Coarse_Manual_Anno, metaData$donor_age_group)
bar_dat <- bar_dat/rowSums(bar_dat)
reorder <- c(3,1,2)
#pdf("Figure1_Barplot_Age.pdf", width=6, height=6)
png("Figure1_Barplot_Age_renamed.png", width=4, height=4, units="in", res=300)
par(mar=c(8, 4, 1,1))
barplot(t(bar_dat[,reorder])*100, col=colours_age[reorder], las=2, main="", xlab="", ylab="Proportion (%)")
dev.off()



# Demographic data - Dot plot by age, coloured by sex
metadata <- read.delim("/cluster/home/tandrews/scripts/LiverMap2.0/Metadata20LiverMapPlusParams.csv", sep=",", header=T)
source("/cluster/home/tandrews/scripts/LiverMap2.0/Map2.2_add_metadata.R")
source("/cluster/home/tandrews/scripts/LiverMap2.0/SubColour_Scheme.R")

thing <- unique(metadata[,"Sample"])
metadata <- get_metadata(thing)
reject <- get_transplant_outcome(thing)

rownames(metadata) <- thing
metadata <- metadata[rownames(metadata) != "C62",]

#metadata <- read.delim("Metadata20LiverMapPlusParams.csv", sep=",", header=T)
#metadata <- metadata[1:(nrow(metadata)-2),]
# Manual stripchart

#pdf("Figure1_Demographics_Dots.pdf", width=4, height=8)
png("Figure1_Demographics_Dots.png", width=4*0.75, height=8*0.75,units="in", res=300)
par(mar=c(0.5,4,0.5,0.5))
xes=rep(1, nrow(metadata))
yes=metadata$age
xes[duplicated(yes)] <- 1.1
plot(xes, yes, pch=19, main="", ylab="Donor Age (years)", yaxt="n", xaxt="n", frame=FALSE, xlab="",
	col=colours_sex[metadata$sex], xlim=c(1,2))
axis(2, at=seq(from=0, to=70, by=10))
abline(h=c(35, 60), col="grey35", lty=3)
legend(1.2, (min(metadata$age, na.rm=T)+35-1)/2, bty="n", c("Young"), yjust=0.5)
legend(1.2, (35+60)/2, bty="n", c("Middle"), yjust=0.5)
legend(1.2, (60+max(metadata$age, na.rm=T)+1)/2, bty="n", c("Elderly"), yjust=0.5)
dev.off()

#par(mar=c(4,1,0,1))
#stripchart(metadata$age, method = "stack", offset = .5, at = .15, pch = 19,main = "", xlab = "Donor Age (years)", frame.plot=FALSE, xaxt="n")
#axis(1, at=seq(from=0, to=70, by=10))
#abline(v=c(35, 60), col="grey35", lty=3)
#legend("topleft", bty="n", c("Young"))
#legend("top", bty="n", c("Middle"))
#legend("topright", bty="n", c("Old"))

metaData <- readRDS("../Merged_EmptyOnly_obj_Map2.2_ImportedClusters_ManualAnno_MetaDataOnly.rds")
metaData[metaData[,"Coarse_Manual_Anno"]=="NonInfMac","Coarse_Manual_Anno"] <- "Kupffer"
metaData[metaData[,"Coarse_Manual_Anno"]=="InfMac","Coarse_Manual_Anno"] <- "Monocyte-like"

assay <- table(metaData$assay_type)
tmp <- table(metaData$assay_type, metaData$sample)
tmp <- rowSums(tmp > 0)
#pdf("Figure1_technology_bars.pdf", width=4, height=4)
png("Figure1_technology_bars_renamed.png", width=4.5, height=3, units="in", res=300)
par(mar=c(4, 2, 0.5, 1))
# cell count by technology
bars <- barplot(assay, horiz=T, xlab="N cells", ylab="", main="", xlim=c(0, 79000), col=colours_tech)
text(assay, bars[,1], paste("(",tmp,")", sep=""), pos=4)
dev.off()


## ---- Supplementary Figures ---- ##

mergedobj <- readRDS("../Merged_EmptyOnly_obj_Map2.2_ImportedClusters_ManualAnno_dimreduce.rds")

# Subcluster groups (from: Map2.2_IndiScaled_Subcluster.R)

cluster_sets <- list(
	Hepatocyte1 = c(0, 1, 4, 5, 9, 19), 
	Hepatocyte2=c(2, 17, 18), 
	NKT=c(3, 11, 12, 16), 
	Stellate=c(13), 
	Cholangiocyte=c(15), 
	Macrophage=c(7, 8), 
	Endo=c(6,10), 
	AntiB=c(14))

orig.clusters = "Coarse_clusters"

mergedobj@meta.data$Subcluster_Group <- as.character(mergedobj@meta.data[,orig.clusters])
for (c_i in 1:length(cluster_sets)){
	id <- names(cluster_sets)[c_i]
	mergedobj@meta.data$Subcluster_Group[mergedobj@meta.data$Subcluster_Group %in% as.character(cluster_sets[[c_i]])] <- id
}

png("SupplFigure1_Subcluster_Groups.png", width=8, height=6, units="in", res=300)
DimPlot(mergedobj, group.by="Subcluster_Group", label=TRUE)
dev.off()

#de <- FindMarkers(mergedobj, group.by="Subcluster_Group", ident.1="Hepatocyte1", ident.2="Hepatocyte2")
# Ribosomal vs Mitochondrial

rps_pct <- PercentageFeatureSet(mergedobj, "^RPS")
rpl_pct <- PercentageFeatureSet(mergedobj, "^RPL")
mt_pct <- PercentageFeatureSet(mergedobj, "^MT")
mergedobj@meta.data$Ribo.pct <- rps_pct+rpl_pct
mergedobj@meta.data$MT.pct <- mt_pct

png("SupplFigure1_MT_umap.png", width=6.5, height=6, units="in", res=300)
FeaturePlot(mergedobj, features="MT.pct", label=FALSE)
dev.off()

d <- split(seq(ncol(mergedobj)), mergedobj@meta.data[,"Cell Type"])
mt_by_type <- sapply(d, function(group) {return(c(
				mean(unlist(mergedobj@meta.data[unlist(group),"MT.pct"])), 
				sd(unlist(mergedobj@meta.data[unlist(group),"MT.pct"])), 
				length(unlist(group))))})

png("SupplFigure1_MT_barplot.png", width=6.5, height=6, units="in", res=300)
par(mar=c(7,4,1,1))
x_loc <- barplot(mt_by_type[1,], names=colnames(mt_by_type), col=new_colour_scheme[,2], xlab="", ylab="MT(%)", las=2, ylim=c(0, 35))
arrows(x_loc, mt_by_type[1,], x_loc, mt_by_type[1,]+2*mt_by_type[2,]/sqrt(mt_by_type[3,]), len=0)
legend("topleft", bty="n", lty=1, c("95% CI"))
dev.off()

png("SupplFigure1_ribo_umap.png", width=6.5, height=6, units="in", res=300)
FeaturePlot(mergedobj, features="Ribo.pct", label=FALSE)
dev.off()

png("ForMikeBoffa_lipoprotein_umaps.png", width=12.5, height=12, units="in", res=300)
FeaturePlot(mergedobj, features=c("APOC1", "APOA2", "LPA", "APOB"), label=FALSE)
dev.off()

# Zonation

Spatial_Portal <- c("ALDOB", "APOA1", "SDS", "CYP2A7", "FABP1", "TTR", "HAL", "AGT", "GLS2", "FGB", "FGA")
Spatial_Central <- c("CYP3A4", "CYP2E1", "ADH1B", "ADH4", "GLUL", "ADH1A", "CES1", "DCXR", "SLCO1B3", "GTPX2", "ADH1C", "GSTA1")

Halpern <- c("GLUL", "CYP1A2", "CYP2E1", "CYP27A1", "HAMP", "IGFBP2", "ARG1", "APOA1", "PCK1", "G6PC", "IGF1", "PIGR", "ACLY")

require(Seurat)
require(ggplot2)

Spatial_genes <- c(Spatial_Central, rev(Spatial_Portal))

png("Figure_Coarse_clusters_spatialgenes.png", width=8, heigh=5, units="in", res=300)
DotPlot(mergedobj, features=Spatial_genes, group.by = "Coarse_clusters") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
png("Figure_Coarse_ManualAnno_spatialgenes.png", width=8.5, heigh=4.5, units="in", res=300)
DotPlot(mergedobj, features=Spatial_genes, group.by = "Coarse_Manual_Anno") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
png("Figure_Haplern_spatialgenes.png", width=6.5, heigh=5, units="in", res=300)
DotPlot(mergedobj, features=Halpern, group.by = "Coarse_Manual_Anno") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

png("Figure_umap_GLUL.png", width=5, height=5, units="in", res=300)
FeaturePlot(mergedobj, features=c("GLUL"))
dev.off()
png("Figure_umap_HAMP.png", width=5, height=5, units="in", res=300)
FeaturePlot(mergedobj, features=c("HAMP"))
dev.off()


# Marker Genes


# Key markers
label_genes <- function(x, this_name) { names(x) <- rep(this_name, length(x)); return(x)}
genes <- c(label_genes(c("VWF", "ENG", "PECAM1", "RAMP3"), "Endo"),
           label_genes(c("DNASE1L3","CLEC1B", "CLEC4G", "CLEC4M"),"LSEC"),
           label_genes(c("CD3D", "CD3E", "TRAC", "CD8A", "IL7R","IL32", "PTPRC"),"CD3Tcells"),
           label_genes(c("CD79A", "CD79B", "JCHAIN", "IGKC", "IGHG4"),"B cell"),
           label_genes(c("GZMB", "GZMK", "FCGR3A", "GNLY", "NKG7"),"NKcell"),
           label_genes(c("HBB", "HBA1", "HBA2"),"Eryth"),
           label_genes(c("MYL9", "ACTA2", "COL6A1", "IGFBP7", "TAGLN"),"Stellate"),
           label_genes(c("MARCO", "CD5L", "C1QC", "HMOX1"),"NonInfMac"),
           label_genes(c("VCAN", "S100A12", "LYZ", "S100A8"),"InfMac"),
           label_genes(c("ANXA4", "KRT7", "KRT8", "KRT18"),"Chol"),
	   label_genes(c("CYP1A2", "CYP2E1", "CYP3A4", "GLUL"),"CentralHep"),
	   label_genes(c("IGF1", "ALB", "APOA1", "CYP2A7"),"PortalHep")
        )

require(ggplot2)
png("SupplFigure1_AnnoType_Dotplot.png", width=11, height=5.5, units="in", res=300)
DotPlot(mergedobj, group.by="Cell Type", features=genes) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

png("SupplFigure1_AnnoCluster_Dotplot.png", width=10, height=5.5, units="in", res=300)
DotPlot(mergedobj, group.by="Coarse_clusters", features=genes) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


# Halpern Comparison
Halpern <- readRDS("/cluster/projects/macparland/TA/ExternalData/HalpernSpatialLocHep/Human_ortho_sig_profiles.rds")
clean_halpern <- Halpern[! (grepl("^RPS", rownames(Halpern)) | grepl("^RPL", rownames(Halpern)) | grepl("^MT-", rownames(Halpern))),]


pseudobulks <- group_rowmeans(mergedobj, mergedobj@meta.data$Coarse_clusters)

hep_clusters <- c(0,1,2,4,5,9,17,18,19)
hep_clusters1 <- c(0,1,4,5,9,19)
hep_clusters2 <- c(2,17,18)

hep_bulks <- pseudobulks[(rownames(pseudobulks) %in% rownames(clean_halpern)),colnames(pseudobulks) %in% hep_clusters]
hep_bulks <- hep_bulks[rowSums(hep_bulks) > 0,]
bulks_range <- apply(hep_bulks, 1, function(x){max(x)-min(x)})
head(hep_bulks[order(bulks_range, decreasing=T),], 50)

hep_bulks <- hep_bulks[bulks_range/rowMeans(hep_bulks) > 0.35 & rowMeans(hep_bulks) > 0.01,]

sync_profiles <- function(to_anno_means, ref_profiles) {
        common_genes <- sort(rownames(to_anno_means)[rownames(to_anno_means) %in% rownames(ref_profiles)])
        to_anno_means <- to_anno_means[match(common_genes, rownames(to_anno_means)),]
        ref_profiles <- ref_profiles[match(common_genes, rownames(ref_profiles)),]

        tmp <- colnames(to_anno_means)
        to_anno_means <- t(apply(to_anno_means, 1, scale));
        colnames(to_anno_means) <- tmp

        tmp <- colnames(ref_profiles)
        ref_profiles <- t(apply(ref_profiles, 1, scale));
        colnames(ref_profiles) <- tmp
        return(list(anno=to_anno_means, ref=ref_profiles));
}

synced <- sync_profiles(hep_bulks, clean_halpern)

require(proxy)
cors <- simil(t(synced$anno), t(synced$ref))
cos <- simil(t(synced$anno), t(synced$ref), method="cosine")



require(gplots)
require(RColorBrewer)
heat_col <- rev(c(rev(brewer.pal(4, "Reds")), "white", brewer.pal(4,"Blues")))

png("SupplHep_vs_Halpern_cosine.png", width=4, height=4, units="in", res=300)
heatmap.2(cos, trace="none", Rowv=TRUE, Colv=FALSE, dendrogram=c("none"), symm=TRUE, scale="none", symbreaks=TRUE, col=heat_col, key.xlab="cosine similarity", key.title="")
dev.off()
png("SupplHep_vs_Halpern_cors.png", width=4, height=4, units="in", res=300)
heatmap.2(cors, trace="none", Rowv=TRUE, Colv=FALSE, dendrogram=c("none"), symm=TRUE, scale="none", symbreaks=TRUE, col=heat_col, key.xlab="Pearson Correlation", key.title="")
dev.off()


## Spatial Zonation score #
#spatial <- read.table("/cluster/home/tandrews/scripts/LiverMap2.0/C73_C1_RESEQ_Varimax1_zonation_geneloadings.txt")
#spatial_zone <- spatial[rownames(spatial) %in% rownames(hep_bulks),]; names(spatial_zone) <- rownames(spatial)[rownames(spatial) %in% rownames(hep_bulks)]
#hep_bulks <- hep_bulks[match(names(spatial_zone), rownames(hep_bulks)),]





