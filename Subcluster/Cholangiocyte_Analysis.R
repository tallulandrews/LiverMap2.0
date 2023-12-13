source("../../scripts/LiverMap2.0/My_R_Scripts.R")
source("../../scripts/LiverMap2.0/Colour_Scheme.R")
source("SubColour_Scheme.R")
source("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/IndividualVariability.R")
require(fgsea)

immune_path <- gmtPathways("../../ExternalData/MSigDb_immune_signatures_c7.all.v7.1.symbols.gmt")
Hallmark_path <- gmtPathways("../../ExternalData/MSigDbHalmarkPathways.gmt")
MSigAll <- gmtPathways("../../ExternalData/MSigDb_curated_c2.all.v7.1.symbols.gmt")
reactome <- gmtPathways("../../ExternalData/ReactomePathways.gmt")
BaderMSig <- gmtPathways("../../ExternalData/BaderLab25Aug2020/Human_MSigdb_August_01_2020_symbol.gmt.txt")
BaderReact <- gmtPathways("../../ExternalData/BaderLab25Aug2020/Human_Reactome_August_01_2020_symbol.gmt.txt")

## Read in data ##
obj <- readRDS("Cholangiocyte_varimax_Subcluster.rds")
obj_full <- obj
obj_all_genes <- readRDS("AllGenes/Cholangiocyte_harmony_Subcluster_Allgenes.rds")


cluster_col = "Coarse_clusters"

manual_cluster_anno <- c("Cholangiocyte", "Lipoprotein", "Lipoprotein", "Cholangiocyte", "Mucosal")


obj@meta.data$Subcluster_Manual <- manual_cluster_anno[obj@meta.data[,cluster_col]]
obj_all_genes@meta.data$Subcluster_Manual <- manual_cluster_anno[obj@meta.data[,cluster_col]]


#saveRDS(obj@meta.data, "Cholangiocyte_fullmetadata.rds")


suppl_tab1 <- get_cluster_summary(obj, samples="sample")
suppl_tab2 <- get_cluster_summary(obj, samples="donor")
write.table(suppl_tab2, file="Cholangiocyte_cluster_summary.csv", sep=",")




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

saveRDS(shiny_obj, "ShinyApps/Cholangiocyte_LatticeObj.rds")
saveRDS(list(expr_mat=expr_mat, metadata=metadata), "ShinyApps/Cholangiocyte_shiny.rds")


## UMAP ##

## DotPlot ##
label_genes <- function(x, this_name) { names(x) <- rep(this_name, length(x)); return(x)}


Cholangiocyte_genes <- label_genes(c("EPCAM", "MUC1", "MUC20", "MUC3A", "KRT19", "KRT9",
				"KRT18", "KRT8", "TROP2", "TROP1",
				 "FGFR2", "TM4SF4", "CLDN1",
				"ANXA4", "MKI67", "TFF3"), "Chol")

Cholangiocyte_genes_dot <- c(label_genes(c("MUC1", "MUC5B", "MUC3A", "TFF3", "SCGB3A1", "SPINK1", "LGALS2", "PIGR", "SLPI", "LYZ"), "Mucus"),
		"ANXA4", "FXYD2", "RPL3", "EEF1A1", label_genes(c("KRT7", "KRT8", "KRT18", "KRT19"), "Keratins"), "DEFB1", "GNB2L1", "AMBP", "NEAT1",  
		label_genes(c("APOC3", "APOA2", "APOA1", "APOC1"), "ApoLipo"), "HP", "MT2A", "AGXT", "EPCAM", "CLDN1")

require(ggplot2)
png("Cholangiocyte_Marker_Dotplot.png", width = 8*1.2, height=4*1.2, units="in", res=300)
DotPlot(obj, features=Cholangiocyte_genes_dot, group.by="Subcluster_Manual")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

png("Cholangiocyte_FinalSubcluster_UMAP.png", width=8*0.65, height=7*0.4, units="in", res=300)
DimPlot(obj, group.by="Subcluster_Manual", label=F)
dev.off()

colours <- get_seurat_colours(obj, "Subcluster_Manual")

png("Cholangiocyte_frequency.png", width=5*0.65, height=7*0.4, units="in", res=300)
par(mar=c(4,1,1,1))
barplot(table(obj@meta.data$Subcluster_Manual), col=colours, horiz=T, las=1, xlab="Cells")
dev.off()

## Frequency Plots ##

png("Cholangiocyte_Sex_Proportions.png", width=6, height=8, units="in", res=300)
barplot_freq(obj, cluster_col="Subcluster_Manual", colours=colours_sex, metadata="donor_sex")
dev.off()
png("Cholangiocyte_Age_Proportions.png", width=6, height=8, units="in", res=300)
barplot_freq(obj, cluster_col="Subcluster_Manual", colours=colours_age, metadata="donor_age_group")
dev.off()
png("Cholangiocyte_Rejection_Proportions.png", width=6, height=8, units="in", res=300)
barplot_freq(obj, cluster_col="Subcluster_Manual", colours=colours_reject, metadata="trans.rejected")
dev.off()

png("Cholangiocyte_donor_Frequency.png", width=9, height=6, units="in", res=300)
indivar_freq(obj@meta.data$Subcluster_Manual, obj@meta.data$donor, obj@meta.data)
dev.off()

indi_expr(obj, tag="Cholangiocyte") 


## DE
for (i in unique(obj@meta.data$Subcluster_Manual)) {
	de <- FindMarkers(obj_all_genes, group.by="Subcluster_Manual", ident.1=i,  logfc.threshold=-Inf)
	write.table(de, file=paste("SubCluster_DE_Tables/Cholangiocyte_", i, "_DE.csv", sep=""), sep=",", row.names=T, col.names=T, quote=T)
}



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

write.table(output, paste("Cholangiocyte_Nebula", assay_type, "DE.csv", sep="_"), sep=",", row.names=FALSE, col.names=TRUE, quote=TRUE)
}

# Subclusters
for (subtype in unique(obj@meta.data$Subcluster_Manual)) {

	assay_type = "3pr"

	require(nebula)
	tofit <- obj_all_genes@assays$RNA@counts[, obj@meta.data$Subcluster_Manual == subtype & obj@meta.data$assay_type==assay_type & !is.na(obj@meta.data$donor_bmi)];
	tofit <- tofit[rowSums(tofit > 0) > 20,]
	metadata <- obj@meta.data[obj@meta.data$Subcluster_Manual == subtype & !is.na(obj@meta.data$donor_bmi) & obj@meta.data$assay_type==assay_type,]
	predictors <- model.matrix(~metadata$donor_age+metadata$donor_sex+metadata$donor_bmi)
	colnames(predictors) <- c("(Intercept)", "Age", "SexM", "BMI")
	
	test <- nebula(tofit, metadata$sample, pred=predictors, offset=colSums(tofit))
	output <- test$summary; 
	output$q_Age <- p.adjust(output$p_Age)
	output$q_Sex <- p.adjust(output$p_Sex)
	output$q_BMI <- p.adjust(output$p_BMI)

	write.table(output, paste("Cholangiocyte_Nebula", assay_type, "Subcluster", subtype, "DE.csv", sep="_"), sep=",", row.names=FALSE, col.names=TRUE, quote=TRUE)

}






#### DO DE - MM ####
# MAST + random effect DE #
require(Seurat)
require(MAST)
require(ggplot2)
obj <- obj_full[,obj_full@meta.data$assay_type == "5pr"]

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

saveRDS(list(sex=sex.output, age=age.output), "Cholangiocyte_MM_DE_5pr.rds")

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

fix_names <- c("Complement cascade", "Bile", "", "Vitamin Met", "Protein phosphorylation",
	"Xenobiotics", "Retinoid", "IGF", "Phototransduction", "", "", "Oxidation",
	"Clotting cascade", "LDL clearance", "Arachidonic acid", "heme scavenging",
	"Plasma lipoprotein", "Scavenger receptors", "Cytochrome P450", 
	"Keratinization", "tRNA", "Antigen presentation", "mt-rRNA", "mt-tRNA")
G <- set.vertex.attribute(richments_de$graph, "name", value=fix_names)

png("Cholangiocyte_DE_Reactome_Sex_pathways.png", width=8, height=8, units="in", res=300)
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
		file="Cholangiocyte_sex_de.rds")



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

fix_names <- c("Complement cascade", "", "Cytochrome P450", "Xenobiotics", 
	"Oxidation", "IGF", "Arachidonic acid", "Vitamin Met", "Protein phosphorylation",
	"Retinoid", "Clotting cascade", "Trigger complement", "phototransduction",
	"Keratinization", "Cell-cell junction", "Steroid Met", "Bile",
	"Platelet", "Plasma lipoprotein", "Fatty acid",
	"HIV", "Interferon", "Phagocytic cup", "Antigen processing",
	"EPH-Ephrin", "Interferon a/b", "Interferon gamma", "ER-Phagosome", "TP53", "FCGR3A phagocytosis",
	"", "Parasite", "Lymphoid/non-Lymphoid", "Antigen Presentation", "Nef",
	"RHO", "TB Infection", "", "DAP12", "")
G <- set.vertex.attribute(richments_de$graph, "name", value=fix_names)

png("Cholangiocyte_DE_Reactome_Rejection_pathways.png", width=8, height=8, units="in", res=300)
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
		file="Cholangiocyte_rejection_de.rds")

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

fix_names <- c("heme scavenging", "Complement cascade", "", "Interferon", "DAP12", "Interferon gamma", 
	"Interferon a/b", "Nef", "", "O-glycosylation mucin", "Antigen Presentation")
G <- set.vertex.attribute(richments_de$graph, "name", value=fix_names)

png("Cholangiocyte_DE_Reactome_Age_pathways.png", width=8, height=8, units="in", res=300)
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
		file="Cholangiocyte_age_de.rds")



# ------------------- DE Vis ------------- #
source("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/IndividualVariability.R")
de_sex <- readRDS("Cholangiocyte_sex_de.rds")
de_age <- readRDS("Cholangiocyte_age_de.rds")
de_rej <- readRDS("Cholangiocyte_rejection_de.rds")

obj <- readRDS("Cholangiocyte_varimax_Subcluster.rds")

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

png(paste("Cholangiocyte", up, "Hi_genes.png", sep="_"), width=8, height=6, units="in", res=150)

par(mfrow=c(2,4))
for(g in top_genes) {
	indivar_DE_vis(obj@assays$RNA@scale.data[g,], obj@meta.data$donor, obj@meta.data, type=type)
#	indivar_DE_vis(obj@assays$RNA@data[g,], obj@meta.data$donor, obj@meta.data)

	title(main=g)
}
dev.off()


top_genes <- rownames(de2[order(de2[,5],decreasing=F),])[1:8]

png(paste("Cholangiocyte", dn, "Hi_genes.png", sep="_"), width=8, height=6, units="in", res=150)

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

png(paste("Cholangiocyte", up, "Hi_pathways.png", sep="_"), width=8, height=6, units="in", res=150)

par(mfrow=c(2,4))
for(i in 1:min(8, nrow(up_pathways))) {
	genes <- unlist(up_pathways[i, "leadingEdge"])
	indivar_DE_vis(obj@assays$RNA@scale.data[genes,], obj@meta.data$donor, obj@meta.data, type=type)
	title(main=up_pathways$pathway[i])
}

dev.off()

png(paste("Cholangiocyte", dn, "Hi_pathways.png", sep="_"), width=8, height=6, units="in", res=150)

par(mfrow=c(2,4))
for(i in 1:min(8, nrow(dn_pathways))) {
	genes <- unlist(dn_pathways[i, "leadingEdge"])
	indivar_DE_vis(obj@assays$RNA@scale.data[genes,], obj@meta.data$donor, obj@meta.data, type=type)
	title(main=dn_pathways$pathway[i])
}

dev.off()




