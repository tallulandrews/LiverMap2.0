
require(Seurat)
source("/cluster/home/tandrews/scripts/LiverMap2.0/Colour_Scheme.R")
source("/cluster/home/tandrews/scripts/LiverMap2.0/My_R_Scripts.R")

dir = "/cluster/projects/macparland/TA/LiverMap2.0/RawData/Analysis/EmptyDrops/Figures"

metaData <- readRDS("../Merged_EmptyOnly_obj_Map2.2_ImportedClusters_ManualAnno_MetaDataOnly.rds")

cluster_2_type <- table(metaData$Coarse_clusters, metaData$Coarse_Manual_Anno)
cluster_2_type <- colnames(cluster_2_type)[ apply(cluster_2_type, 1, function(x){which(x == max(x))}) ]
cluster_2_colour <- type_2_colour(cluster_2_type)

# Simpson Index
freq_by_sample <- table(metaData$Coarse_clusters, metaData$sample)
simpson_index <- function(freq_vec) {
	freq_vec <- freq_vec/sum(freq_vec)
	return(sum(freq_vec^2))
}

simpson_by_cluster <- apply(freq_by_sample, 1, simpson_index)
null <- simpson_index(table(metaData$sample));

pdf("Figure2_Simpson_by_cluster.pdf", width=5, height=5)
barplot(simpson_by_cluster, col=cluster_2_colour, horiz=2, las=2)
dev.off()

# Frequency boxplot
# Proportion of sample for each cell-type
# Boxplots with dots

typefreq_by_sample <- table(metaData$sample, metaData$Coarse_Manual_Anno)
typefreq_by_sample <- typefreq_by_sample/rowSums(typefreq_by_sample)*100
bxp <- boxplot(as.vector(typefreq_by_sample)~rep(colnames(typefreq_by_sample), each=nrow(typefreq_by_sample)))
reorder <- rev(order(bxp$stats[3,]))

typefreq_by_sample <- typefreq_by_sample[,reorder]

pdf("Figure2_Freq_by_Sample_boxplot.pdf", width=10, height=6)
par(mar=c(8,4,1,1))
bxp <- boxplot(as.vector(typefreq_by_sample)~factor(rep(colnames(typefreq_by_sample), each=nrow(typefreq_by_sample)), levels=colnames(typefreq_by_sample)), 
	labels=colnames(typefreq_by_sample), col=type_2_colour(colnames(typefreq_by_sample)), las=2, 
	xlab="", ylab="Frequency (%)", outline=F)
points(jitter(rep(1:15, each=34)), as.vector(typefreq_by_sample), pch=16)
dev.off()

# Sample-sample similarity by cluster
mergedobj <- readRDS("../Merged_EmptyOnly_obj_Map2.2_ImportedClusters_ManualAnno.rds")

hvgs <- mergedobj@misc$repeated_hvgs

obj_5pr <- mergedobj[rownames(mergedobj) %in% hvgs, mergedobj@meta.data$assay_type == "5pr"]
obj_3pr <- mergedobj[rownames(mergedobj) %in% hvgs, mergedobj@meta.data$assay_type == "3pr"]

dim(mergedobj)
cluster_5pr_pseudo <- get_pseudobulk(obj_5pr@assays$RNA@counts, factor(obj_5pr@meta.data$Coarse_clusters), factor(obj_5pr@meta.data$donor))
cluster_3pr_pseudo <- get_pseudobulk(obj_3pr@assays$RNA@counts, obj_3pr@meta.data$Coarse_clusters, obj_3pr@meta.data$donor)

id <- factor(unlist(lapply(strsplit(colnames(cluster_3pr_pseudo), "_"), function(x){x[[1]]})))
id2 <- factor(unlist(lapply(strsplit(colnames(cluster_5pr_pseudo), "_"), function(x){x[[1]]})))

d_5pr <- as.dist(1-cor(cluster_5pr_pseudo, method="pearson"))
d_3pr <- as.dist(1-cor(cluster_3pr_pseudo, method="pearson"))

scores <- c();
stderr <- c();

cross_scores <- c();

for (i in unique(id)) {
	score_3pr <- mean(as.matrix(1-d_3pr)[id==i, id==i], na.rm=T);
	crossscore_3pr <- mean(as.matrix(1-d_3pr)[id==i, id!=i], na.rm=T);
	score_5pr <- mean(as.matrix(1-d_5pr)[id2==i, id2==i], na.rm=T);
	crossscore_5pr <- mean(as.matrix(1-d_5pr)[id2==i, id2!=i], na.rm=T);
	std_3pr <- sd(as.matrix(1-d_3pr)[id==i, id==i], na.rm=T)/sqrt(sum(id==i, na.rm=T))
	std_5pr <- sd(as.matrix(1-d_5pr)[id2==i, id2==i], na.rm=T)/sqrt(sum(id2==i, na.rm=T))
	if (is.na(score_5pr)) {
		scores <- c(scores, score_3pr);
		stderr <- c(stderr, std_3pr)
		cross_scores <- c(cross_scores, crossscore_3pr);
	} else {
		scores <- c(scores, (score_3pr+score_5pr)/2)
		stderr <- c(stderr, sqrt( (std_3pr^2*sqrt(sum(id==i)) + std_5pr^2*sum(id2==i))/(sum(id==i)+sum(id2==i)) ))
		cross_scores <- c(cross_scores, (crossscore_3pr+crossscore_5pr)/2)
	}
}

pdf("Figure2_Similarity_by_cluster.pdf", width=5, height=5)
bar_loc <- barplot(scores, name=unique(id), xlab="Cross Donor Correlation (average)", xlim=c(0,1), las=1, horiz=T)
arrows(scores, bar_loc, scores+2*stderr,  bar_loc, angle=90)
legend("bottomright", lty=1, c("95% CI"), bty="n")
abline(v=mean(cross_scores), lty=2)
dev.off()


# Most variable genes by cell-type
obj_5pr <- mergedobj[ , mergedobj@meta.data$assay_type == "5pr"]
obj_3pr <- mergedobj[ , mergedobj@meta.data$assay_type == "3pr"]

all_cell_types <- sort(unique(obj_3pr@meta.data$Coarse_Manual_Anno))

all_scores=c()
all_means=c()
all_types=c()
all_genes=c();

all_cell_types <- all_cell_types[! all_cell_types %in% c("Eryth", "ErythHep")]

for (i in all_cell_types) {
	# 3pr
	this_3pr <- obj_3pr[,obj_3pr@meta.data$Coarse_Manual_Anno == i]
	means_3pr <- group_rowmeans(this_3pr@assays$RNA@data, this_3pr@meta.data$sample)
	btw_var <- group_rowvars(means_3pr, rep("A", ncol(means_3pr)))
	within_var <- rowMeans(group_rowvars(this_3pr@assays$RNA@data, this_3pr@meta.data$sample))
	scores_3pr <- log2(btw_var[,1]/within_var)
	if (sum(obj_5pr@meta.data$Coarse_Manual_Anno == i) > 10) {
		# 5pr
		this_5pr <- obj_5pr[,obj_5pr@meta.data$Coarse_Manual_Anno == i]
		means_5pr <- group_rowmeans(this_5pr@assays$RNA@data, this_5pr@meta.data$sample)
		btw_var <- group_rowvars(means_5pr, rep("A", ncol(means_5pr)))
		within_var <- rowMeans(group_rowvars(this_5pr@assays$RNA@data, this_5pr@meta.data$sample))
		scores_5pr <- log2(btw_var[,1]/within_var)
		# Score = log2(btw/within)
		all_means <- c(all_means, (my_rowMeans(means_3pr)+my_rowMeans(means_5pr))/2)
		all_scores <- c(all_scores, (scores_3pr+scores_5pr)/2)
	} else {
		all_means <- c(all_means, my_rowMeans(means_3pr))
		all_scores <- c(all_scores, scores_3pr)
	}
	all_types <- c(all_types, rep(i, nrow(this_3pr)))
	all_genes <- c(all_genes, rownames(obj_3pr));
}

summary(all_scores)

keep <- is.finite(all_scores)

to.label <- quantile(all_scores[keep], prob=1-30/sum(keep));
to.label <- keep & all_scores > to.label

tmp <- unique(all_types[keep])

pdf("Figure2_HighlyVariable_scatter.pdf", width=12, height=5)
plot(all_means[keep], all_scores[keep], pch=16, col=type_2_colour(all_types[keep]))
text(x=all_means[to.label], all_scores[to.label], names(all_scores)[to.label], pos=3)
legend("bottom", tmp, horiz=T, fill=type_2_colour(tmp), bty="n")
dev.off()

#expr_means_3pr <- group_rowmeans(obj_3pr@assays$RNA@data, obj_3pr@meta.data$Coarse_Manual_Anno)
#expr_means_5pr <- group_rowmeans(obj_5pr@assays$RNA@data, obj_5pr@meta.data$Coarse_Manual_Anno)

#var_group_3pr <- group_rowvars(obj_3pr@assays$RNA@data, obj_3pr@meta.data$Coarse_Manual_Anno)
#var_group_5pr <- group_rowvars(obj_5pr@assays$RNA@data, obj_5pr@meta.data$Coarse_Manual_Anno)

#overall_var_3pr <- group_rowvars(obj_3pr@assays$RNA@data, rep("all", ncol(obj_3pr)) )
#overall_var_5pr <- group_rowvars(obj_5pr@assays$RNA@data, rep("all", ncol(obj_5pr)) )

#l2fc_3pr <- log2(var_group_3pr/overall_var_3pr[,1])
#l2fc_5pr <- log2(var_group_5pr/overall_var_5pr[,1])

#l2fc_3pr[!is.finite(l2fc_3pr)] <- 0
#l2fc_5pr[!is.finite(l2fc_5pr)] <- 0

#l2fc_3pr <- l2fc_3pr[,match(colnames(l2fc_5pr), colnames(l2fc_3pr))]

# Consistent direction btw 3pr and 5pr
#disagree <- sign(l2fc_3pr) != sign(l2fc_5pr)
#l2fc_3pr[disagree] <- 0
#l2fc_5pr[disagree] <- 0

#l2fc_total <- l2fc_3pr + l2fc_5pr

#sort(rownames(l2fc_total)[which(apply(l2fc_total, 1, max) > 10)]) # Fibers (Stellate), Immunoglobulins (Bcells), CC (NK cell), CXC genes (cholangiocytes), mucus cholangiocytes


# xaxt = mean expression
# yaxt = log2( across sample var / within sample var)
# label top genes
# colour by cell-type



