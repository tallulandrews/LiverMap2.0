
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

#pdf("Figure1_UMAP_Cluster_and_Type.pdf", width=8, height=8)
png("Figure1_UMAP_Cluster_and_Type.png", width=8, height=8, units="in", res=300)
print( DimPlot(mergedobj, reduction="umap", group.by="Cell Type", pt.size=.1)+scale_color_manual(values=new_colour_scheme[,2])+annotate("text", x=umap_lab_pos[1,], y=umap_lab_pos[2,], label=colnames(umap_lab_pos), colour="black") )
dev.off()

# UMAP by demograph data

mergedobj@meta.data$donor_sex <- factor(metaData$donor_sex)

#pdf("Figure1_UMAP_Sex.pdf", width=8, height=8)
png("Figure1_UMAP_Sex.png", width=6, height=6, units="in", res=300)
print( DimPlot(mergedobj, reduction="umap", group.by="donor_sex", pt.size=.1)+scale_color_manual(values=colours_sex) )
dev.off()

#pdf("Figure1_UMAP_Age.pdf", width=8, height=8)
png("Figure1_UMAP_Age.png", width=6, height=6, units="in", res=300)
print( DimPlot(mergedobj, reduction="umap", group.by="donor_age_group", pt.size=.1)+scale_color_manual(values=colours_age) )
dev.off()

#pdf("Figure1_UMAP_Tech.pdf", width=8, height=8)
png("Figure1_UMAP_Tech.png", width=6, height=6, units="in", res=300)
print( DimPlot(mergedobj, reduction="umap", group.by="assay_type", pt.size=.1)+scale_color_manual(values=colours_tech) )
dev.off()

# Barplot by demographic data

bar_dat <- table(metaData$Coarse_Manual_Anno, mergedobj@meta.data$donor_sex)
bar_dat <- bar_dat/rowSums(bar_dat)
#pdf("Figure1_Barplot_Sex.pdf", width=6, height=6)
png("Figure1_Barplot_Sex.png", width=6, height=6, units="in", res=300)
par(mar=c(8, 4, 1,1))
barplot(t(bar_dat)*100, col=colours_sex, las=2, main="", xlab="", ylab="Proportion (%)")
dev.off()

bar_dat <- table(metaData$Coarse_Manual_Anno, mergedobj@meta.data$assay_type)
bar_dat <- bar_dat/rowSums(bar_dat)
#pdf("Figure1_Barplot_Tech.pdf", width=6, height=6)
png("Figure1_Barplot_Tech.png", width=6, height=6, units="in", res=300)
par(mar=c(8, 4, 1,1))
barplot(t(bar_dat)*100, col=colours_tech, las=2, main="", xlab="", ylab="Proportion (%)")
dev.off()

bar_dat <- table(metaData$Coarse_Manual_Anno, mergedobj@meta.data$donor_age_group)
bar_dat <- bar_dat/rowSums(bar_dat)
reorder <- c(3,1,2)
#pdf("Figure1_Barplot_Age.pdf", width=6, height=6)
png("Figure1_Barplot_Age.png", width=6, height=6, units="in", res=300)
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
png("Figure1_Demographics_Dots.png", width=4, height=8,units="in", res=300)
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

pdf("Figure1_technology_bars.pdf", width=2, height=2)
png("Figure1_technology_bars.png", width=2, height=2, units="in", res=300)
# cell count by technology
assay <- table(metaData$assay_type)
bars <- barplot(assay, horiz=T, xlab="N cells", ylab="", main="", xlim=c(0, 75000), col=colours_tech)
tmp <- table(metaData$assay_type, metaData$sample)
text(assay, bars[,1], paste("(",tmp,")", sep=""), pos=4)
dev.off()




