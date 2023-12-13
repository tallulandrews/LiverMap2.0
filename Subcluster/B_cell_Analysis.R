## set up enrichment ##
source("../../scripts/LiverMap2.0/My_R_Scripts.R")
source("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/IndividualVariability.R")
require(fgsea)
require(Seurat)

immune_path <- gmtPathways("../../ExternalData/MSigDb_immune_signatures_c7.all.v7.1.symbols.gmt")
Hallmark_path <- gmtPathways("../../ExternalData/MSigDbHalmarkPathways.gmt")
MSigAll <- gmtPathways("../../ExternalData/MSigDb_curated_c2.all.v7.1.symbols.gmt")
reactome <- gmtPathways("../../ExternalData/ReactomePathways.gmt")
BaderMSig <- gmtPathways("../../ExternalData/BaderLab25Aug2020/Human_MSigdb_August_01_2020_symbol.gmt.txt")
BaderReact <- gmtPathways("../../ExternalData/BaderLab25Aug2020/Human_Reactome_August_01_2020_symbol.gmt.txt")

## Read in data ##
obj <- readRDS("AntiB_varimax_Subcluster.rds")
obj_full <- obj
obj_all_genes <- readRDS("AllGenes/AntiB_harmony_Subcluster_Allgenes.rds")



cluster_col = "Coarse_clusters"

manual_cluster_anno <- c("Plasmablast", "IgK+IgG+", "IgK+IgA+", "IgL+IgG+",
	"IgL+IgA+", "Naive", "Debris")

obj@meta.data$Subcluster_Manual <- factor(manual_cluster_anno[obj@meta.data[,cluster_col]], levels=manual_cluster_anno)
obj_all_genes@meta.data <- obj@meta.data
#saveRDS(obj@meta.data, "AntiB_fullmetadata.rds")

suppl_tab1 <- get_cluster_summary(obj, samples="sample")
suppl_tab2 <- get_cluster_summary(obj, samples="donor")
write.table(suppl_tab2, file="Bcell_cluster_summary.csv", sep=",")



subtype_cols <- get_seurat_colours(obj, "Subcluster_Manual"); names(subtype_cols) <- levels(factor(obj@meta.data$Subcluster_Manual))
### Object for Shiny ###
shiny_obj <- obj_all_genes;
detect_rate <- group_rowmeans(shiny_obj@assays$RNA@counts > 0, shiny_obj@meta.data$Subcluster_Manual)
exclude_genes <- apply(detect_rate, 1, max) < 0.005
shiny_obj@meta.data <- shiny_obj@meta.data[,c("nCount_RNA","nFeature_RNA","percent.mt","Phase", "donor", "sample", "donor_sex", "donor_age", "Subcluster_Manual")]
shiny_obj@reductions <- obj@reductions
shiny_obj@assays$RNA@scale.data <- matrix(0)
shiny_obj@assays$RNA@var.features <- c(0)
shiny_obj <- shiny_obj[!exclude_genes,]

expr_mat <- shiny_obj@assays$RNA@data
metadata <- shiny_obj@meta.data
metadata$UMAP1 <- shiny_obj@reductions$umap@cell.embeddings[,1]
metadata$UMAP2 <- shiny_obj@reductions$umap@cell.embeddings[,2]
metadata$cell_colour <- subtype_cols[match(shiny_obj@meta.data$Subcluster_Manual, names(subtype_cols))]


saveRDS(shiny_obj, "ShinyApps/Bcell_LatticeObj.rds")
saveRDS(list(expr_mat=expr_mat, metadata=metadata), "ShinyApps/Bcell_shiny.rds")

# ----------------------------------- #


png("Supplementary_AntiB_map.png", width=7, height=5, units="in", res=150)
DimPlot(obj, group.by="Subcluster_Manual")
dev.off()

require(ggplot2)
label_genes <- function(g, name) { names(g)<- rep(name, length(g)); return(g)}
png("Supplementary_AntiB_markers.png", width=9, height=5, units="in", res=150)
DotPlot(obj, features=c("IGHM", "IGHD",  
				label_genes("IGKC", "LK"), 
				label_genes(c("IGLC1", "IGLC2", "IGLC3", "IGLL5"),"LL"),
				label_genes(c("IGHA1", "IGHA2"), "HA"), 
				label_genes(c("IGHG1", "IGHG2", "IGHG3", "IGHG4"), "HG"),
				label_genes(c("TOP2A", "BRIC5", "MKI67", "CDK1"), "CC")),
				group.by="Subcluster_Manual") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


## Frequency Plots ##

png("AntiB_Sex_Proportions.png", width=6, height=8, units="in", res=300)
barplot_freq(obj, cluster_col="Subcluster_Manual", colours=colours_sex, metadata="donor_sex")
dev.off()
png("AntiB_Age_Proportions.png", width=6, height=8, units="in", res=300)
barplot_freq(obj, cluster_col="Subcluster_Manual", colours=colours_age, metadata="donor_age_group")
dev.off()
png("AntiB_Rejection_Proportions.png", width=6, height=8, units="in", res=300)
barplot_freq(obj, cluster_col="Subcluster_Manual", colours=colours_reject, metadata="trans.rejected")
dev.off()

png("AntiB_donor_Frequency.png", width=9, height=6, units="in", res=300)
indivar_freq(obj@meta.data$Subcluster_Manual, obj@meta.data$donor, obj@meta.data)
dev.off()

indi_expr(obj, tag="AntiB") 



### ------ AUCell vs Guilliams ------ ###
genesets <- readRDS("../../ExternalData/GuilliamsCellPaper/Guilliams_Genesets.rds")

require(AUCell)
cells_rankings <- AUCell_buildRankings(obj_all_genes@assays$RNA@scale.data)
cells_AUC <- AUCell_calcAUC(genesets, cells_rankings, aucMaxRank=5, nCores=1)

tmp <- cells_AUC@assays@data$AUC
tmp[tmp < 0] <- 0

assignment<-apply(tmp, 2, function(x) { 

							out <- which(x==max(x)); 
							if (length(out) > 1) {return("None")} 
							else {return (rownames(cells_AUC@assays@data$AUC)[out])}
							}
			)
obj@meta.data$AUCell <- assignment
DimPlot(obj, split.by="AUCell", group.by="Subcluster_Manual")

obj@meta.data$Plasma <- cells_AUC@assays@data$AUC["HumanPlasmacell",]+cells_AUC@assays@data$AUC["Plasmacells",]
obj@meta.data$B_cell <- cells_AUC@assays@data$AUC["HumanBcell",]
obj@meta.data$pDC <- cells_AUC@assays@data$AUC["pDC",] + cells_AUC@assays@data$AUC["HumanpDC",]


png("Bcell_Guilliams_AUCell.png", width=6, height=8/3, units="in", res=300)
FeaturePlot(obj, feature=c("Plasma", "B_cell"))
dev.off()



#### DO DE - MM - NEBULA ####

for(assay_type in c("5pr", "3pr")) {

require(nebula)
tofit <- obj_all_genes@assays$RNA@counts[,obj@meta.data$assay_type==assay_type & !is.na(obj@meta.data$donor_bmi)];
tofit <- tofit[rowSums(tofit > 0) > 20,]
metadata <- obj@meta.data[!is.na(obj@meta.data$donor_bmi) & obj@meta.data$assay_type==assay_type,]
predictors <- model.matrix(~metadata$donor_age+metadata$donor_sex+metadata$donor_bmi)
colnames(predictors) <- c("(Intercept)", "Age", "SexM", "BMI")

test <- nebula(tofit, metadata$sample, pred=predictors, offset=colSums(tofit))
output <- test$summary; 
output$q_Age <- p.adjust(output$p_Age)
output$q_Sex <- p.adjust(output$p_Sex)
output$q_BMI <- p.adjust(output$p_BMI)

write.table(output, paste("Bcell_Nebula", assay_type, "DE.csv", sep="_"), sep=",", row.names=FALSE, col.names=TRUE, quote=TRUE)
}




#### DO DE - MM ####
# MAST + random effect DE #
require(Seurat)
require(MAST)
require(ggplot2)
obj <- obj_full[,obj_full@meta.data$assay_type == "3pr"]

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

saveRDS(list(sex=sex.output, age=age.output), "AntiB_MM_DE_3pr.rds")


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

fix_names <- c("Xenobiotics", "Oxidation", "IGF", "Cytochrome P450", "Protein phosphorylation",
	"Integrin", "Plasma lipoprotein", "Clotting cascade", "Vitamin Met", "Bile",
	"", "", "", "BRAF", "", "Retinoid", "RAF-BRAF", "RAS", "", "", "FCGR3A phagocytosis",
	"Parasite", "Interleukin-2", "Anti-inflammatory", "Leishmania survival",
	"Interferon gamma", "IL10", "PD-1", "TCR phosphorylation",
	"Phagocytosis", "DAP12", "calcium mobilization", "Interferon a/b",
	"MAPK", "2nd messenger", "", "BCR activation", "Lymphoid/NonLymphoid interaction",
	"FCGR activation")
G <- set.vertex.attribute(richments_de$graph, "name", value=fix_names)

png("AntiB_DE_Reactome_Sex_pathways.png", width=8, height=8, units="in", res=300)
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
		file="AntiB_sex_de.rds")



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

fix_names <- c("IGF", "Complement cascade", "", "protein phosphorylation", 
	"Plasma lipoprotein", "Vitamin Met", "Steroid Met", "Scavenger Receptors",
	"Oxidation", "Fatty acid", "Xenobiotic", "Phototransduction", "",
	"Lipid Met", "Retinoid", "Cytochrome P450", "Clotting cascade", 
	"Plasma lipoprotein", "", "PPARA", "FCGR3A phagocytosis", "",
	"Parasite", "RHO", "Splicing", "EPHB", "DAP12", "RHO Effector",
	"ROBO", "TB Infection", "Axon guidance", "HIV", "", "TCR", "EPH-Ephrin",
	"", "", "Splicing", "", "")
G <- set.vertex.attribute(richments_de$graph, "name", value=fix_names)

png("AntiB_DE_Reactome_Rejection_pathways.png", width=8, height=8, units="in", res=300)
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
		file="AntiB_rejection_de.rds")

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

fix_names <- c("IGF", "protein phosphorylation", "Complement cascade", "", 
	"Steroid Met", "Lipid Met", "Plasma lipoprotein", "Oxidation", "", "PPARA",
	"Vitamin Met", "Fatty acid", "Xenobiotic", "", "Clotting cascade", "Scavenger Receptors",
	"phototransduction", "Cholesterol SREBP", "Retinoid", "", "Cell Cycle",
	"", "M Phase", "Translation", "", "Mitotic", "TCR signaling", "RNA Met",
	"Translation", "SLITs/ROBOs", "Anaphase", "BCR signaling", "", "ROBO", 
	"Splicing", "TCR", "HIV", "", "", "" )
G <- set.vertex.attribute(richments_de$graph, "name", value=fix_names)

png("AntiB_DE_Reactome_Age_pathways.png", width=8, height=8, units="in", res=300)
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
		file="AntiB_age_de.rds")

# ------------------- DE Vis ------------- #
source("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/IndividualVariability.R")
de_sex <- readRDS("AntiB_sex_de.rds")
de_age <- readRDS("AntiB_age_de.rds")
de_rej <- readRDS("AntiB_rejection_de.rds")

obj <- readRDS("AntiB_varimax_Subcluster.rds")

type ="sex"

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

png(paste("AntiB", up, "Hi_genes.png", sep="_"), width=8, height=6, units="in", res=150)

par(mfrow=c(2,4))
for(g in top_genes) {
	indivar_DE_vis(obj@assays$RNA@scale.data[g,], obj@meta.data$donor, obj@meta.data, type=type)
#	indivar_DE_vis(obj@assays$RNA@data[g,], obj@meta.data$donor, obj@meta.data)

	title(main=g)
}
dev.off()


top_genes <- rownames(de2[order(de2[,5],decreasing=F),])[1:8]

png(paste("AntiB", dn, "Hi_genes.png", sep="_"), width=8, height=6, units="in", res=150)

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

png(paste("AntiB", up, "Hi_pathways.png", sep="_"), width=8, height=6, units="in", res=150)

par(mfrow=c(2,4))
for(i in 1:min(8, nrow(up_pathways))) {
	genes <- unlist(up_pathways[i, "leadingEdge"])
	indivar_DE_vis(obj@assays$RNA@scale.data[genes,], obj@meta.data$donor, obj@meta.data, type=type)
	title(main=up_pathways$pathway[i])
}

dev.off()

png(paste("AntiB", dn, "Hi_pathways.png", sep="_"), width=8, height=6, units="in", res=150)

par(mfrow=c(2,4))
for(i in 1:min(8, nrow(dn_pathways))) {
	genes <- unlist(dn_pathways[i, "leadingEdge"])
	indivar_DE_vis(obj@assays$RNA@scale.data[genes,], obj@meta.data$donor, obj@meta.data, type=type)
	title(main=dn_pathways$pathway[i])
}

dev.off()


###------------------ Figure for Paper ----------------- ###

source("../../scripts/LiverMap2.0/My_R_Scripts.R")
source("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/IndividualVariability.R")
source("C:/Users/tandrews/Documents/UHNSonya/scripts/LiverMap2.0/Colour_Scheme.R")
require(Seurat)
require(ggplot2)

Contamination_genes <- c("ALB", "SERPINA1", "APOA1", "FGA", "CYP3A5", "CYP2D6", "ASGR1")
Prolif_genes <- c("MKI67", "BIRC5", "TOP2A", "CDK1", "HMGB2", "CDC20", "CDK4")
RBC_genes <- c("HBB", "HBA1", "HBA2", "HBD")

AntiB_genes_dot <- c("IGHA1", "IGHA2", 
			"IGHG1", "IGHG2", "IGHG3", "IGHG4", 
			"IGLC1", "IGLC2", "IGLC3",
			"IGKC", "JCHAIN",
			"IGHGP", "IGHD", "IGHM", "IGLL5", 
			"CD79A", "CD79B",
			"RRM2", "CCND2", "CDKN3", "CDC20", "TUBA1B",
			"HMGN2", "HMGB1", "HMGB2",
			"MKI67", "CD19", "CD38")
#			"SERPINA1", "SAA1", "FGA", "FGB", "APOA2")

## Read in data ##
obj <- readRDS("AntiB_varimax_Subcluster.rds")

cluster_col = "Coarse_clusters"

manual_cluster_anno <- c("Plasmablast", "IgK+IgG+", "IgK+IgA+", "IgL+IgG+",
	"IgL+IgA+", "Naive", "Doublet")

obj@meta.data$Manual_Subcluster_Anno <- manual_cluster_anno[obj@meta.data[,cluster_col]]

obj <- obj[,obj@meta.data$Manual_Subcluster_Anno != "Doublet"]

# Subclustering plot
pdf("Figure_AntiB_UMAP.pdf", width=12, height=12)
DimPlot(obj, group.by="Manual_Subcluster_Anno", label=T)
dev.off()

# dotplot of markers
pdf("Figure_AntiB_Dotplot.pdf", width=9, height=4)
DotPlot(obj, features=unique(c(AntiB_genes_dot)), group.by="Manual_Subcluster_Anno")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# Freq by sample
pdf("Figure_AntiB_Indifreq.pdf", width=9, height=4)
indivar_freq(obj@meta.data$Manual_Subcluster_Anno, obj@meta.data$donor, obj@meta.data)
dev.off()


# Cross donor correlation similarity - with horizontal line for global similarity between different cell-types
pdf("Figure_AntiB_CrossSim.pdf", width=9, height=4)
require(cluster)
mat <- get_pseudobulk_means(obj@assays$RNA@data, obj@meta.data[,"Manual_Subcluster_Anno"], obj@meta.data[,"donor"])
ns <- table(obj@meta.data[,"donor"])
d <- as.dist(1-cor(mat, method="pearson"))
id <- factor(unlist(lapply(strsplit(colnames(mat), "_"), function(x){x[[1]]})))

scores <- c();
stderr <- c();
for (i in unique(id)) {
	scores <- c(scores, mean(as.matrix(1-d)[id==i, id==i]))
	stderr <- c(stderr, sd(as.matrix(1-d)[id==i, id==i])/sqrt(sum(id==i))) # Fix 4 May 2021
}
colours <- get_seurat_colours(obj, "Manual_Subcluster_Anno")
par(mar=c(6,4,1,1))
bar_loc <- barplot(scores, name=unique(id), ylab="Cross Donor Correlation (average)", ylim=c(0,1), las=2, col=colours)
arrows(bar_loc, scores, bar_loc, scores+2*stderr, angle=90)
legend("topright", lty=1, c("95% CI"), bty="n")
dev.off()


# sex / age DE genes/pathways plot
de_sex <- readRDS("AntiB_sex_de.rds")
de_age <- readRDS("AntiB_age_de.rds")
de_rej <- readRDS("AntiB_rejection_de.rds")

#top_de_genes <- de_sex$pseudo_de[order(de_sex$pseudo_de[,"pval"]),]

gene = "TNFRSF17"
indivar_DE_vis(obj@assays$RNA@scale.data[gene,], obj@meta.data$donor, obj@meta.data, type="sex"); title(main=gene)
pval =  signif(de_sex$pseudo_de[gene,"pval"], digits=1); legend("topright", paste("p =", pval), bty="n", col="white")


gene = "CXCR3"
indivar_DE_vis(obj@assays$RNA@scale.data[gene,], obj@meta.data$donor, obj@meta.data, type="sex"); title(main=gene)
pval =  signif(de_sex$pseudo_de[gene,"pval"], digits=1); legend("topright", paste("p =", pval), bty="n", col="white")



### Supplementary - Key Markers ###
Spatial <- readRDS("C:/Users/tandrews/Documents/UHNSonya/Spatial/C73_C1/C73_C1_Processed_Spatial_Seurat.rds")

key_genes = c("JCHAIN", "CD79A", "CD79B", "IGKC", "IGHG1", "IGHA1", "IGLC2")
SpatialFeaturePlot(Spatial, features=key_genes)


for (gene in key_genes) {
	# UMAP with gene level
	FeaturePlot(obj, features=gene)
	# Spatial of gene
	SpatialFeaturePlot(Spatial, features=gene))
	# HPA of gene
}


































#------------------------ OLD -------------#

## Rejection ##
varimax_comp1 = 1
varimax_comp2 = -3
exclude_cluster = c("6")

expr_mat <- obj@assays$RNA@scale.data
label <- obj@meta.data$trans.rejected

de <- t(apply(expr_mat, 1, function(x) {
		res = t.test(x[label=="Y"], x[label=="N"])
		return(c(res$estimate, res$p.value))
		}))
de <- cbind(de, p.adjust(de[,3], "fdr"))
de <- cbind(de, -1*log(de[,3])*sign(de[,1]-de[,2]))
colnames(de) <- c("Y", "N", "pval", "qval", "score")

de_clean <- de[!grepl("^RPS", rownames(de)),]
de_clean <- de_clean[!grepl("^RPL", rownames(de_clean)),]

richments_de <- do_fgsea(de_clean[,5], reactome, fdr=0.05)

require(igraph)
names(V(richments_de$graph))

fix_names <- c("Insulin-like Growth Factor", "Complement cascade", "", "Protein Phosphorylation",
	"Lipoprotein", "Vitamins", "Steroids", "Scavenger receptors", "Oxidation", "Fatty acid",
	"Xenobiotic Met", "Phototransduction", "", "Lipids", "Cytochrome P450", "Retinoid",
	"Clotting Cascade", "", "", "PPARA", "RHO GTPases", "", 
	"Phagocytosis", "Parasite infection", "", "DAP12", "EPHB Signaling", "",
	"ROBO", "TB infection", "Axon guidance", "HIV", "", "",
	"EPH-Ephrin", "HIV", "TCR signaling", "", "", "Splicing")
G <- set.vertex.attribute(richments_de$graph, "name", value=fix_names)

png("Bcell_DE_Reactome_Rejection_pathways.png", width=8, height=8, units="in", res=300)
set.seed(2910)
plot(G, vertex.color=richments_de$vertex_col, vertex.size=richments_de$vertex_size)
dev.off()



varimax1_cell_scores <- obj@reductions$varimax@cell.embeddings[,abs(varimax_comp1)]
varimax2_cell_scores <- obj@reductions$varimax@cell.embeddings[,abs(varimax_comp2)]

vmax_gene_cor_1 <- simil(matrix(obj@reductions$varimax@cell.embeddings[,abs(varimax_comp1)], nrow=1), obj@assays$RNA@scale.data, method="Pearson")
vmax_gene_cor_1 <- as.vector(vmax_gene_cor_1[,order(vmax_gene_cor_1)]); names(vmax_gene_cor_1) <- rownames(obj@assays$RNA@scale.data)
vmax_gene_cor_1_clean <- vmax_gene_cor_1[!grepl("^RPS", names(vmax_gene_cor_1))]
vmax_gene_cor_1_clean <- vmax_gene_cor_1_clean[!grepl("^RPL", names(vmax_gene_cor_1_clean))]
## These are worse than before ##

varimax_gene_score1 = rowSums(obj@assays$RNA@scale.data * 
		(varimax1_cell_scores * sign(varimax_comp1)))
varimax_gene_score2 = rowSums(obj@assays$RNA@scale.data * 
		(varimax2_cell_scores * sign(varimax_comp2))) # this is all hepatocyte stuff.



#richments <- do_fgsea(varimax_gene_score1, Hallmark_path, fdr=0.05)
richments <- do_fgsea(varimax_gene_score1, reactome, fdr=0.05)
#richments <- do_fgsea(varimax_gene_score1, BaderMSig, fdr=0.05)
#richments <- do_fgsea(varimax_gene_score1, BaderReact, fdr=0.05)

require(igraph)
names(V(richments$graph))
fix_names <- c("Oxidations", "MT-rRNA", "Xenobiotic Met", "MT-tRNA", 
	"Clotting Cascade", "Phototransduction", "Platelets", "Platelet Ca2+",
	"Fatty acids", "Retinoid", "Lipoprotein", "Vitamins", "",
	"Insulin-like Growth Factor", "Cytochrome P450", "Lipids", "",
	"Bile", "NR1H2/3 signaling", "Arachidonic acid", "Chaperones", 
	"MT-translation", "Hedgehog ligand", "FCGR phagocytosis", 
	"APC Mitosis", "BCR activation", "APC Mitosis", "SLITs/ROBOs",
	"", "NF-kB activation", "Phagocytic actin", "CFTR Disease",
	"Th2 Macrophage", "", "Unfolded Protein Response", "IL10 synthesis",
	"Phagocytic phospolipids", "", "Membrane Targetting", "Translation")
G <- set.vertex.attribute(richments$graph, "name", value=fix_names)

png("Bcell_Varimax_Reactome_Rejection_pathways.png", width=8, height=8, units="in", res=300)
set.seed(2910)
plot(G, vertex.color=richments$vertex_col, vertex.size=richments$vertex_size)
dev.off()



#Phagocytosis genes#
# term: "Role of phospholipids in phagocytosis", 
#  "Fcgamma receptor (FCGR) dependent phagocytosis", 
#  "Regulation of actin dynamics for phagocytic cup formation"
tmp <- data.frame(richments$rich)
tmpde <- data.frame(richments_de$rich)
phago_genes <- unique( c(unlist(tmp[tmp[,1] == "Role of phospholipids in phagocytosis",8]), 
  unlist(tmp[tmp[,1] == "Fcgamma receptor (FCGR) dependent phagocytosis",8]),
  unlist(tmp[tmp[,1] == "Regulation of actin dynamics for phagocytic cup formation",8]),
  unlist(tmpde[tmpde[,1] == "FCGR3A-mediated phagocytosis",8]),
  unlist(tmpde[tmpde[,1] == "Leishmania phagocytosis",8]),
  unlist(tmpde[tmpde[,1] == "RHO GTPases Activate WASPs and WAVEs",8])
	))
actual_phago_genes <- phago_genes <- unique( c(
  unlist(tmpde[tmpde[,1] == "FCGR3A-mediated phagocytosis",8]),
  unlist(tmpde[tmpde[,1] == "Leishmania phagocytosis",8]),
  unlist(tmpde[tmpde[,1] == "RHO GTPases Activate WASPs and WAVEs",8])
	))

score <- Matrix::colSums(obj@assays$RNA@scale.data[rownames(obj@assays$RNA@scale.data) %in% actual_phago_genes,])

boxplot(score ~ obj@meta.data$trans.rejected, col=c("salmon", "grey85", "dodgerblue"), 
		names=c("Not Rejected", "Minor Rejection", "Rejection"), notch=T, xlab="", ylab="Phagocytic Expression")

require(Seurat)
plot_obj <- obj
plot_obj@assays$RNA@data <- plot_obj@assays$RNA@scale.data
VlnPlot(plot_obj, features=phago_genes[13:24], group.by="trans.rejected")

c(tmp[tmp[,1] == "FCGR3A-mediated IL10 synthesis",8])


## Age ##
varimax_comp = c(1, 

## Gender ##