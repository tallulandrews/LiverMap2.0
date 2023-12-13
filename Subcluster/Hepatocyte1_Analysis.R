source("../../scripts/LiverMap2.0/My_R_Scripts.R")
source("../../scripts/LiverMap2.0/Colour_Scheme.R")
require(fgsea)
require(Seurat)
require(ggplot2)

immune_path <- gmtPathways("../../ExternalData/MSigDb_immune_signatures_c7.all.v7.1.symbols.gmt")
Hallmark_path <- gmtPathways("../../ExternalData/MSigDbHalmarkPathways.gmt")
MSigAll <- gmtPathways("../../ExternalData/MSigDb_curated_c2.all.v7.1.symbols.gmt")
reactome <- gmtPathways("../../ExternalData/ReactomePathways.gmt")
BaderMSig <- gmtPathways("../../ExternalData/BaderLab25Aug2020/Human_MSigdb_August_01_2020_symbol.gmt.txt")
BaderReact <- gmtPathways("../../ExternalData/BaderLab25Aug2020/Human_Reactome_August_01_2020_symbol.gmt.txt")

## Read in data ##
obj <- readRDS("Hepatocyte1_varimax_Subcluster.rds")
set.seed(1920)
obj_full_subset <- obj[, sample(1:ncol(obj), 0.2*ncol(obj))] ## For DE subset to 20% of all cells.

cluster_col = "Coarse_clusters"

## Use original Clusters & make relevant dot plots

obj@meta.data$Global_Coarse_clusters = factor(obj@meta.data[,37])

spatial_marker_genes <- c("CYP3A4", "CYP2E1", "ADH1B", "GLUL", "ADH4",
	"DCXR", "FTL", "ADH1A", "GSTA1", "ADH1C",
	"GPX2", "HPD", "SULT2A1", "ALDH1L1",
	"ALDOB", "FABP1", "AGT", "FGB", "CYP2A7", "HAL", "FGA", "SDS",
	"FGG", "CYP2A6", "CYP3A5", "HAMP")
Prolif_genes <- c("MKI67", "BIRC5", "TOP2A", "CDK1", "HMGB2", "CDC20", "CDK4")
RBC_genes <- c("HBB", "HBA1", "HBA2", "HBD")

Metal <- c("MT1X", "MT1G", "MT2A", "MT1E", "MT1H")

Lipoprotein <- c("APOA2", "APOC1", "APOC3", "APOA1")
mito <- c("MT-ND3", "MTCO3", "MT-CO1", "MT-CO2", "MT-ND4", "MT-ATP6")

DotPlot(obj, group.by="Global_Coarse_clusters", features=c(spatial_marker_genes, Prolif_genes, Metal, Lipoprotein, mito))+
	theme(axis.text.x=element_text(angle = -90, hjust = 0))

DimPlot(obj, group.by="Global_Coarse_clusters")

##Pseudotime
#require(slingshot)
#require(destiny)
#require(Seurat)

#Type_DimPlot(obj, type_col = "Coarse_Manual_Anno")


#mat <- obj@assays$RNA@scale.data

##sds <- slingshot(Embeddings(obj, "umap"), clusterLabels = obj$Coarse_clusters)



## Manual Annotation & Individual Variability ##

cluster_col = "Global_Coarse_clusters"

manual_cluster_anno <- c("Portal", "Portal", "Interzonal", "Central",
	"Central", "InterProlif")

obj@meta.data$Subcluster_Manual <- manual_cluster_anno[obj@meta.data[,cluster_col]]

#saveRDS(obj@meta.data, "Hepatocyte1_fullmetadata.rds")


suppl_tab1 <- get_cluster_summary(obj, samples="sample")
suppl_tab2 <- get_cluster_summary(obj, samples="donor")
write.table(suppl_tab2, file="Hepatocyte1_cluster_summary.csv", sep=",")



### Object for Shiny ###
shiny_obj <- obj;
detect_rate <- group_rowmeans(shiny_obj@assays$RNA@counts > 0, shiny_obj@meta.data$Subcluster_Manual)
exclude_genes <- apply(detect_rate, 1, max) < 0.005
shiny_obj@meta.data <- shiny_obj@meta.data[,c("nCount_RNA","nFeature_RNA","percent.mt","Phase", "donor", "sample", "donor_sex", "donor_age", "Global_Coarse_clusters", "Subcluster_Manual")]
shiny_obj@reductions <- obj@reductions
shiny_obj@assays$RNA@scale.data <- matrix(0)
shiny_obj@assays$RNA@var.features <- c(0)
shiny_obj <- shiny_obj[!exclude_genes,]
expr_mat <- shiny_obj@assays$RNA@data
metadata <- shiny_obj@meta.data
metadata$UMAP1 <- shiny_obj@reductions$umap@cell.embeddings[,1]
metadata$UMAP2 <- shiny_obj@reductions$umap@cell.embeddings[,2]

saveRDS(shiny_obj, "ShinyApps/Hepatocyte1_Lattice_obj.rds")
saveRDS(list(expr_mat=expr_mat, metadata=metadata), "ShinyApps/Hepatocyte1_shiny.rds")












png("Hepatocyte1_Sex_Proportions.png", width=6, height=8, units="in", res=300)
barplot_freq(obj, cluster_col="Subcluster_Manual", colours=colours_sex, metadata="donor_sex")
dev.off()
png("Hepatocyte1_Age_Proportions.png", width=6, height=8, units="in", res=300)
barplot_freq(obj, cluster_col="Subcluster_Manual", colours=colours_age, metadata="donor_age_group")
dev.off()
png("Hepatocyte1_Rejection_Proportions.png", width=6, height=8, units="in", res=300)
barplot_freq(obj, cluster_col="Subcluster_Manual", colours=colours_reject, metadata="trans.rejected")
dev.off()

png("Hepatocyte1_donor_Frequency.png", width=9, height=6, units="in", res=300)
indivar_freq(obj@meta.data$Subcluster_Manual, obj@meta.data$donor, obj@meta.data)
dev.off()

indi_expr(obj, tag="Hepatocyte1") 

## Do PSEUDOTIME ##
set.seed(1920)

require(slingshot)
require(scater)
downsample <- c()
nmax=1500
for (i in unique(obj@meta.data$Global_Coarse_clusters)) {
	if (sum(obj@meta.data$Global_Coarse_clusters == i) < nmax) {
		downsample <- c(downsample, which(obj@meta.data$Global_Coarse_clusters == i))
	} else {
		subset <- sample(which(obj@meta.data$Global_Coarse_clusters == i), nmax)
		downsample <- c(downsample, subset)
	}
}

obj_subset <- obj[,downsample]
obj_subset_sce <- as.SingleCellExperiment(obj_subset)
p1 = plotExpression(obj_subset_sce, features="GLUL", x="Global_Coarse_clusters")
p2 = plotExpression(obj_subset_sce, features="CYP2E1", x="Global_Coarse_clusters")
p3 = plotExpression(obj_subset_sce, features="CYP2A6", x="Global_Coarse_clusters")
p4 = plotExpression(obj_subset_sce, features="CYP2A7", x="Global_Coarse_clusters")

p1+p2+p3+p4

require(CycleMix)
CC_genes <- HGeneSets$Tirosh

plotExpression(obj_subset_sce, features="RRM2", x="Global_Coarse_clusters")
plotExpression(obj_subset_sce, features="TOP2A", x="Global_Coarse_clusters")
plotExpression(obj_subset_sce, features="CDK1", x="Global_Coarse_clusters")
plotExpression(obj_subset_sce, features="TUBB", x="Global_Coarse_clusters")
plotExpression(obj_subset_sce, features="TUBB4B", x="Global_Coarse_clusters")

reducedDim(obj_subset_sce, "umap") <- obj_subset@reductions$umap@cell.embeddings

sim <- slingshot(obj_subset_sce, clusterLabels = 'Global_Coarse_clusters', reducedDim = 'umap')

plot(reducedDims(sim)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sim), lwd=2, col='black', type="lineages") # looks like a star??

#Note CYP2A6 is missing from subcluster HVGs.

require(destiny)
mat <- obj_subset@assays$RNA@scale.data[rownames(obj_subset) %in% obj_subset@misc$repeated_hvgs,]
dm<- destiny::DiffusionMap(as.matrix(mat))

#
#### DO DE - MM ####
# MAST + random effect DE #
require(Seurat)
require(MAST)
require(ggplot2)
obj <- obj_full_subset[,obj_full_subset@meta.data$assay_type == "5pr"]

obj_sce <- as.SingleCellExperiment(obj)
tmp <- as.matrix(obj@assays$RNA@scale.data); rownames(tmp) <- rownames(obj); colnames(tmp) <- colnames(obj)
obj_sca <- FromMatrix(tmp, colData(obj_sce), rowData(obj_sce), class="SingleCellAssay")

zlm.out2 <- zlm(~donor_sex + donor_age + (1/sample), obj_sca)
coefAndCI <- as.data.frame(summary(zlm.out2, logFC=FALSE)$datatable)

res_sex <- coefAndCI[coefAndCI[,3] == "donor_sexM" & !is.na(coefAndCI[,6]),]
res_age <- coefAndCI[coefAndCI[,3] == "donor_age" & !is.na(coefAndCI[,6]),]

max_coef.age <- aggregate(res_age$coef, by=list(res_age$primerid), function(x) {x[abs(x) == max(abs(x))]})
max_coef.sex <- aggregate(res_sex$coef, by=list(res_sex$primerid), function(x) {x[abs(x) == max(abs(x))]})


zlm.lr.sex <- lrTest(zlm.out2, "donor_sex")
zlm.lr.age <- lrTest(zlm.out2, "donor_age")

pval.sex <- zlm.lr.sex[,,3][,"hurdle"]
pval.age <- zlm.lr.age[,,3][,"hurdle"]

sex.output <- cbind(max_coef.sex, pval.sex[match( max_coef.sex[,1], names(pval.sex))])
sex.output$qval <- p.adjust(sex.output[,3], method="fdr")
age.output <- cbind(max_coef.age, pval.age[match( max_coef.age[,1], names(pval.age))])
age.output$qval <- p.adjust(age.output[,3], method="fdr")

print(sum(age.output$qval < 0.05))
print(sum(sex.output$qval < 0.05))

# >7000 age, > 4000 sex out of 10,000 total genes - same if use normalized, 3' only, 5' only, or scaled data.

saveRDS(list(sex=sex.output, age=age.output), "Hepatocyte1_MM_DE_5pr.rds")



#### DO DE ####
## Sex ##
set.seed(38201)
expr_mat <- obj@assays$RNA@scale.data
label <- obj@meta.data$donor_sex


de_sex <- t(apply(expr_mat, 1, function(x) {
		res = t.test(x[label=="F"], x[label=="M"])
		return(c(res$estimate, res$p.value))
		}))
de_sex <- cbind(de_sex, p.adjust(de_sex[,3], "fdr"))
de_sex <- cbind(de_sex, -1*log(de_sex[,3])*sign(de_sex[,1]-de_sex[,2]))
colnames(de_sex) <- c("F", "M", "pval", "qval", "score")

de_clean <- de_sex[!grepl("^RPS", rownames(de_sex)),]
de_clean <- de_clean[!grepl("^RPL", rownames(de_clean)),]

richments_de <- do_fgsea(de_clean[,5], reactome, fdr=0.05)

require(igraph)
names(V(richments_de$graph))

fix_names <- c("TCR", "Lymphoid-non-Lymphoid", "CD28", "TCR", "CD3", "2ndmessensger", "PD-1", "Interferon gamma", 
	"DAP12", "ZAP-70 synapse", "Antigen processing", "DAP12", "MHCI", "ER-Phagosome", "Interferon",
	"MHCII", "Interferon a/b", "Immune System", "actin phagocytic cup", "HIV", "PPARA", "Fatty acid", "PPARA", "Xenobiotics",
	"Protein localization", "AA", "phototransduction", "Cytocrhome P450", "Glucuronidation", "Peroxisome", "Sulfur AA", 
	"Lipoprotein", "Xenobiotics", "Glyoxylate", "Peroxisome", "oxidation", "Xenobiotics")
G <- set.vertex.attribute(richments_de$graph, "name", value=fix_names)

png("Hepatocyte1_DE_Reactome_Sex_pathways.png", width=8, height=8, units="in", res=300)
set.seed(2910)
plot(G, vertex.color=richments_de$vertex_col, vertex.size=richments_de$vertex_size)
dev.off()

## Peudobulk ##

pseudo_sex <- get_pseudobulk_means(expr_mat, obj@meta.data$donor_sex, obj@meta.data$sample)
is.F <- grepl("^F", colnames(pseudo_sex))
de_pseudo <- t(apply(pseudo_sex, 1, function(x) {
		res = t.test(x[is.F], x[!is.F])
		return(c(res$estimate, res$p.value))
		}))
de_pseudo <- cbind(de_pseudo, p.adjust(de_pseudo[,3], "fdr"))
de_pseudo <- cbind(de_pseudo, -1*log(de_pseudo[,3])*sign(de_pseudo[,1]-de_pseudo[,2]))
colnames(de_pseudo) <- c("F", "M", "pval", "qval", "score")

saveRDS(list(de=de_sex, richments=richments_de, fix_names = fix_names, 
		pseudo_de=de_pseudo, pseudo_table=pseudo_sex), 
		file="Hepatocyte1_sex_de.rds")



## Rejection ##
set.seed(1672)
expr_mat <- obj@assays$RNA@scale.data
label <- obj@meta.data$trans.rejected


de_rej <- t(apply(expr_mat, 1, function(x) {
		res = t.test(x[label=="Y"], x[label=="N"])
		return(c(res$estimate, res$p.value))
		}))
de_rej <- cbind(de_rej, p.adjust(de_rej[,3], "fdr"))
de_rej <- cbind(de_rej, -1*log(de_rej[,3])*sign(de_rej[,1]-de_rej[,2]))
colnames(de_rej) <- c("Y", "N", "pval", "qval", "score")

de_clean <- de_rej[!grepl("^RPS", rownames(de_rej)),]
de_clean <- de_clean[!grepl("^RPL", rownames(de_clean)),]

richments_de <- do_fgsea(de_clean[,5], reactome, fdr=0.05)

require(igraph)
names(V(richments_de$graph))

fix_names <- c("Oxidation", "AA", "Vitamins", "Xenobiotics", "Xenobiotics", "Liproprotein", "Complement cascade", "IGF", 
	"phototransduction", "fat-soluble vitamins", "protein phosphorylation", "Complement cascade", "lipoprotein", "Retinoid", "ETC",
	"Fatty acid", "ATP synthesis", "Bile acid", "Steroids", "Citric Acid Cycle", "RHO GTPase", "Antigen processing", 
	"Clathrin endocytosis", "ERBB2", "Clathrin endocytosis", "ERBB2", "EGFR", "FCGR", "FCGR3A Phagocytosis", "Phagocytosis", 
	"parasite", "actin phagocytic cup", "Adaptive Immune", "PI3K/AKT", "PI3K", "Negative PI3K/AKT", "", "EGFR", "DAP12", 
	"Lymphoid-non-Lymphoid")
G <- set.vertex.attribute(richments_de$graph, "name", value=fix_names)

png("Hepatocyte1_DE_Reactome_Rejection_pathways.png", width=8, height=8, units="in", res=300)
set.seed(2910)
plot(G, vertex.color=richments_de$vertex_col, vertex.size=richments_de$vertex_size)
dev.off()

## Peudobulk ##

pseudo_reject <- get_pseudobulk_means(expr_mat, obj@meta.data$trans.rejected, obj@meta.data$sample)
reject <- grepl("^Y", colnames(pseudo_reject))
notreject <- grepl("^N", colnames(pseudo_reject))
de_pseudo <- t(apply(pseudo_reject, 1, function(x) {
		res = t.test(x[reject], x[notreject])
		return(c(res$estimate, res$p.value))
		}))
de_pseudo <- cbind(de_pseudo, p.adjust(de_pseudo[,3], "fdr"))
de_pseudo <- cbind(de_pseudo, -1*log(de_pseudo[,3])*sign(de_pseudo[,1]-de_pseudo[,2]))
colnames(de_pseudo) <- c("Reject", "NotReject", "pval", "qval", "score")

saveRDS(list(de=de_rej, richments=richments_de, fix_names = fix_names, 
		pseudo_de=de_pseudo, pseudo_table=pseudo_reject), 
		file="Hepatocyte1_rejection_de.rds")

## Age ##
set.seed(3771)
expr_mat <- obj@assays$RNA@scale.data
label <- obj@meta.data$donor_age_group


de_age <- t(apply(expr_mat, 1, function(x) {
		res = t.test(x[label=="elderly"], x[label=="young"])
		return(c(res$estimate, res$p.value))
		}))
de_age <- cbind(de_age, p.adjust(de_age[,3], "fdr"))
de_age <- cbind(de_age, -1*log(de_age[,3])*sign(de_age[,1]-de_age[,2]))
colnames(de_age) <- c("elderly", "young", "pval", "qval", "score")

de_clean <- de_age[!grepl("^RPS", rownames(de_age)),]
de_clean <- de_clean[!grepl("^RPL", rownames(de_clean)),]

richments_de <- do_fgsea(de_clean[,5], reactome, fdr=0.05)

require(igraph)
names(V(richments_de$graph))

fix_names <- c("TCR", "", "Translation", "HIV", "Antigen processing", "Adaptive Immune", "HIV", "APC/CC", "Cell cycle", "CD28", "Metaphase", "Anaphase", "actin phagocytic cup",
	"Infection", "BCR", "ER-Phagosome", "ACP/CDC20", "FCGR3A", "phagocytosis", "parasite", "Cile", "SREBF", "Lipoprotein", "",
	"PPARA", "Cytochrome P450", "tRNA", "PPARA", "Clotting Cascade", "mt-tRNA", "Lipoprotein", "oxidation", "Retinoid",
	"fat-soluble vitamin", "mt-rRNA", "Xenobiotics", "phototransduction", "Vitamins", "IGF", "protein phosphorylation")
G <- set.vertex.attribute(richments_de$graph, "name", value=fix_names)

png("Hepatocyte1_DE_Reactome_Age_pathways.png", width=8, height=8, units="in", res=300)
set.seed(2910)
plot(G, vertex.color=richments_de$vertex_col, vertex.size=richments_de$vertex_size)
dev.off()

## Peudobulk ##

pseudo_age <- get_pseudobulk_means(expr_mat, obj@meta.data$donor_age_group, obj@meta.data$sample)
young <- grepl("^young", colnames(pseudo_age))
elderly <- grepl("^elderly", colnames(pseudo_age))
de_pseudo <- t(apply(pseudo_age, 1, function(x) {
		res = t.test(x[elderly], x[young])
		return(c(res$estimate, res$p.value))
		}))
de_pseudo <- cbind(de_pseudo, p.adjust(de_pseudo[,3], "fdr"))
de_pseudo <- cbind(de_pseudo, -1*log(de_pseudo[,3])*sign(de_pseudo[,1]-de_pseudo[,2]))
colnames(de_pseudo) <- c("Elderly", "Young", "pval", "qval", "score")

saveRDS(list(de=de_age, richments=richments_de, fix_names = fix_names, 
		pseudo_de=de_pseudo, pseudo_table=pseudo_age), 
		file="Hepatocyte1_age_de.rds")

# ------------------- DE Vis ------------- #
source("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/IndividualVariability.R")
de_sex <- readRDS("Hepatocyte1_sex_de.rds")
de_age <- readRDS("Hepatocyte1_age_de.rds")
de_rej <- readRDS("Hepatocyte1_rejection_de.rds")

obj <- readRDS("Hepatocyte1_varimax_Subcluster.rds")

type ="reject"

#top DE
if (type =="age") {
	de <- de_age
	up = "elderly"
	dn = "young"
}
if (type =="sex") {
	de <- de_sex
	up = "F"
	dn = "M"
}
if (type =="reject") {
	de <- de_rej
	up = "Y"
	dn = "N"
}

#de <- de_rej$de
de2 <- de$de
de2 <- de2[de2[,5] >0 & de2[,1] >0 | de2[,5] < 0 & de2[,2] >0,]

top_genes <- rownames(de2[order(de2[,5],decreasing=T),])[1:8] # F = N, F = M, F = young

png(paste("Hepatocyte1", up, "Hi_genes.png", sep="_"), width=8, height=6, units="in", res=150)

par(mfrow=c(2,4))
for(g in top_genes) {
	indivar_DE_vis(obj@assays$RNA@scale.data[g,], obj@meta.data$donor, obj@meta.data, type=type)
#	indivar_DE_vis(obj@assays$RNA@data[g,], obj@meta.data$donor, obj@meta.data)

	title(main=g)
}
dev.off()


top_genes <- rownames(de2[order(de2[,5],decreasing=F),])[1:8]

png(paste("Hepatocyte1", dn, "Hi_genes.png", sep="_"), width=8, height=6, units="in", res=150)

par(mfrow=c(2,4))
for(g in top_genes) {
	indivar_DE_vis(obj@assays$RNA@scale.data[g,], obj@meta.data$donor, obj@meta.data, type=type)
#	indivar_DE_vis(obj@assays$RNA@data[g,], obj@meta.data$donor, obj@meta.data)

	title(main=g)
}
dev.off()

#top pathways

#de = de_age

up_pathways <- de$richments$rich[de$richments$rich$NES > 0,] # F, Y, elderly
up_pathways <- up_pathways[rev(1:nrow(up_pathways)),]
dn_pathways <- de$richments$rich[de$richments$rich$NES < 0,] # M, N, young

png(paste("Hepatocyte1", up, "Hi_pathways.png", sep="_"), width=8, height=6, units="in", res=150)

par(mfrow=c(2,4))
for(i in 1:min(8, nrow(up_pathways))) {
	genes <- unlist(up_pathways[i, "leadingEdge"])
	indivar_DE_vis(obj@assays$RNA@scale.data[genes,], obj@meta.data$donor, obj@meta.data, type=type)
	title(main=up_pathways$pathway[i])
}

dev.off()

png(paste("Hepatocyte1", dn, "Hi_pathways.png", sep="_"), width=8, height=6, units="in", res=150)

par(mfrow=c(2,4))
for(i in 1:min(8, nrow(dn_pathways))) {
	genes <- unlist(dn_pathways[i, "leadingEdge"])
	indivar_DE_vis(obj@assays$RNA@scale.data[genes,], obj@meta.data$donor, obj@meta.data, type=type)
	title(main=dn_pathways$pathway[i])
}

dev.off()

