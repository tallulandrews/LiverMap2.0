source("../../scripts/LiverMap2.0/My_R_Scripts.R")
source("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/IndividualVariability.R")
source("SubColour_Scheme.R")
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
obj <- readRDS("Endo_varimax_Subcluster.rds")
obj_full <- obj
obj_all_genes <- readRDS("AllGenes/Endo_harmony_Subcluster_Allgenes.rds")



cluster_col = "Coarse_clusters"


manual_cluster_anno <- c("cvLSEC", "cvLSEC", "cvLSEC", "ppLSEC", "cvLSEC", "cvEndo", 
	"ppLSEC", "interLSEC", "HepContam", "VasEndo")


manual_cluster_anno <- c("cvLSEC", "cvLSEC", "cvLSEC", "ppLSEC", "cvLSEC", "cvEndo", 
	"ppLSEC", "interLSEC", "HepContam", "Arterial")


obj@meta.data$Subcluster_Manual <- manual_cluster_anno[obj@meta.data[,cluster_col]]
obj_all_genes@meta.data$Subcluster_Manual <- obj@meta.data$Subcluster_Manual

#saveRDS(obj@meta.data, "Endo_fullmetadata.rds")

suppl_tab1 <- get_cluster_summary(obj, samples="sample")
suppl_tab2 <- get_cluster_summary(obj, samples="donor")
write.table(suppl_tab2, file="Endo_cluster_summary.csv", sep=",")




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

saveRDS(shiny_obj, "ShinyApps/Endo_LatticeObj.rds")
saveRDS(list(expr_mat=expr_mat, metadata=metadata), "ShinyApps/Endo_shiny.rds")



# ---------------------------#

# Marker Genes
label_genes <- function(x, this_name) { names(x) <- rep(this_name, length(x)); return(x)}

Endo_genes_dot <- c(label_genes(c("FCN2", "FCN3", "STAB1", "STAB2", "CLEC1B", "CLEC4G", "CLEC4M", "CRHBP"), "cvLSEC"),
			label_genes(c("SPARCL1", "TM4SF1", "CLEC14A", "ID1", "IGFBP7", "MGP", "ADRIF", "S100A6", "AQP1"), "ppLSEC"), "VWF", 
			label_genes(c("RSPO3", "ACKR1", "WNT2"), "cvEndo"), "INMT", "PLAC8",
			label_genes(c("PODXL", "PLVAP", "CD34"), "Arterial"), "GSN", "RBP7",
			label_genes(c("ENG", "PECAM1",  "DNASE1L3", "TIMP3", "LIFR", "C7"), "Endo"), "RAMP3"
			 )
#"CCL21", "CST3", "ALB", "APOA1","PDPN"

obj_all_genes@meta.data$cluster_anno <- paste(obj_all_genes@meta.data$Subcluster_Manual, " (", obj@meta.data[,cluster_col],")", sep="")
png("ForFigure_Endo_Dotplot.png", width=10, height=4, units="in", res=300)
DotPlot(obj_all_genes, features=Endo_genes_dot, group.by="cluster_anno") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


DimPlot(obj, group.by="Subcluster_Manual")

png("ForFigure_Endo_Dimplot.png", width=5, height=4.5, units="in", res=300)
DimPlot(obj, group.by=cluster_col)
dev.off()

de <- run_wilcox(obj_all_genes, binned=obj@meta.data$Subcluster_Manual, ident.1="HepContam")


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

obj@meta.data$cvEndo <- cells_AUC@assays@data$AUC["cvEndo",]
obj@meta.data$pvEndo <- cells_AUC@assays@data$AUC["pvEndo",]
obj@meta.data$LSEC <- cells_AUC@assays@data$AUC["LSEC",]
obj@meta.data$VSMC <- cells_AUC@assays@data$AUC["HumanVSMC",]
obj@meta.data$Endo <- cells_AUC@assays@data$AUC["HumanEndo",]
obj@meta.data$Fibroblast <- cells_AUC@assays@data$AUC["HumanFibro",]


png("Endo_Guilliams_AUCell.png", width=6, height=8, units="in", res=300)
FeaturePlot(obj, feature=c("cvEndo", "pvEndo", "LSEC", "VSMC", "Endo", "Fibroblast"))
dev.off()




## Frequency Plots ##
png("Endo_Sex_Proportions.png", width=6, height=8, units="in", res=300)
barplot_freq(obj, cluster_col="Subcluster_Manual", colours=colours_sex, metadata="donor_sex")
dev.off()
png("Endo_Age_Proportions.png", width=6, height=8, units="in", res=300)
barplot_freq(obj, cluster_col="Subcluster_Manual", colours=colours_age, metadata="donor_age_group")
dev.off()
png("Endo_Rejection_Proportions.png", width=6, height=8, units="in", res=300)
barplot_freq(obj, cluster_col="Subcluster_Manual", colours=colours_reject, metadata="trans.rejected")
dev.off()

png("Endo_donor_Frequency.png", width=9, height=6, units="in", res=300)
indivar_freq(obj@meta.data$Subcluster_Manual, obj@meta.data$donor, obj@meta.data)
dev.off()

indi_expr(obj, tag="Endo") 


## DE
for (i in unique(obj@meta.data$Subcluster_Manual)) {
	de <- FindMarkers(obj_all_genes, group.by="Subcluster_Manual", ident.1=i,  logfc.threshold=-Inf)
	write.table(de, file=paste("SubCluster_DE_Tables/Endo_", i, "_DE.csv", sep=""), sep=",", row.names=T, col.names=T, quote=T)
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

	write.table(output, paste("Endo_Nebula", assay_type, "DE.csv", sep="_"), sep=",", row.names=FALSE, col.names=TRUE, quote=TRUE)
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

saveRDS(list(sex=sex.output, age=age.output), "Endo_MM_DE_5pr.rds")

## Sex ##

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

fix_names <- c("Protein Phosphorylation", "Oxidation", "EPH-Ephrin", "IGF", "AA metabolism",
	"Cytokine Signaling", "RNA Met", "Adaptive Immune", "ETC & ATP", "FCGR", "CD28", "PD-1",
	"Immunological synapse", "CD3 TCR", "Lymphoid-non-Lymphoid interactions", "mt-tRNA", "mt-rRNA", 
	"TCR", "", "tRNA", "TCR", "rRNA")
G <- set.vertex.attribute(richments_de$graph, "name", value=fix_names)

png("Endo_DE_Reactome_Sex_pathways.png", width=8, height=8, units="in", res=300)
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
		file="Endo_sex_de.rds")





## Rejection ##
set.seed(2891)
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

fix_names <- c("Xenobiotics", "Oxidation", "GPCR", "GPCR", "G alpha",
	"Interleukin signaling", "Infection", "Cytokine Signaling", "Vascular wall", "Anti-inflammatory",
	"Leishmania survival", "GPCR ligand", "", "MB infection", "Rhodospin-like", "DAP12", 
	"Interferon gamma", "FCGR", "Chemokine", "Innate Immune", "Immunological synapse",
	"PD-1", "Neutrophil", "TCR Phosphorylation", "Adaptive Immune", "CD28", "2nd messenger", 
	"TCR signaling", "Lymphoid/Non-Lymphoid", "TCR")
G <- set.vertex.attribute(richments_de$graph, "name", value=fix_names)

png("Endo_DE_Reactome_Rejection_pathways.png", width=8, height=8, units="in", res=300)
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
		file="Endo_rejection_de.rds")



## Age ##
set.seed(7881)
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

fix_names <- c("Complement cascade", "", "lipoprotein", "Scavenger Receptors", "lipoprotein", "", "C4 + C2", 
	"IGF", "Steroids", "mt-rRNA", "cholesterol", "initial complement", "Heme scavenging", "mt-tRNA", 
	"protein phosphorylation", "SREBF", "NR1H2/3", "Vitamins", "Glycerophospholipid", "Cholesterol", 
	"MyD88", "", "TLR2", "TLR1-TLR2", "TLR6-TLR2", "pathogen DNA sensors", "HIV", "TCR", "AUF1",
	"mRNA stability", "TLR9", "NfkB/TLR", "MyD88", "TLR7/8", "MAPK6/MAPK4", "Interleukin17", "TCR",
	"MAPK", "Oxidative-Senescence", "Senescence-Secretory"
	)
G <- set.vertex.attribute(richments_de$graph, "name", value=fix_names)

png("Endo_DE_Reactome_Age_pathways.png", width=8, height=8, units="in", res=300)
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
		file="Endo_age_de.rds")



# ------------------- DE Vis ------------- #
source("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/IndividualVariability.R")
de_sex <- readRDS("Endo_sex_de.rds")
de_age <- readRDS("Endo_age_de.rds")
de_rej <- readRDS("Endo_rejection_de.rds")

obj <- readRDS("Endo_varimax_Subcluster.rds")

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

png(paste("Endo", up, "Hi_genes.png", sep="_"), width=8, height=6, units="in", res=150)

par(mfrow=c(2,4))
for(g in top_genes) {
	indivar_DE_vis(obj@assays$RNA@scale.data[g,], obj@meta.data$donor, obj@meta.data, type=type)
#	indivar_DE_vis(obj@assays$RNA@data[g,], obj@meta.data$donor, obj@meta.data)

	title(main=g)
}
dev.off()


top_genes <- rownames(de2[order(de2[,5],decreasing=F),])[1:8]

png(paste("Endo", dn, "Hi_genes.png", sep="_"), width=8, height=6, units="in", res=150)

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

png(paste("Endo", up, "Hi_pathways.png", sep="_"), width=8, height=6, units="in", res=150)

par(mfrow=c(2,4))
for(i in 1:min(8, nrow(up_pathways))) {
	genes <- unlist(up_pathways[i, "leadingEdge"])
	indivar_DE_vis(obj@assays$RNA@scale.data[genes,], obj@meta.data$donor, obj@meta.data, type=type)
	title(main=up_pathways$pathway[i])
}

dev.off()

png(paste("Endo", dn, "Hi_pathways.png", sep="_"), width=8, height=6, units="in", res=150)

par(mfrow=c(2,4))
for(i in 1:min(8, nrow(dn_pathways))) {
	genes <- unlist(dn_pathways[i, "leadingEdge"])
	indivar_DE_vis(obj@assays$RNA@scale.data[genes,], obj@meta.data$donor, obj@meta.data, type=type)
	title(main=dn_pathways$pathway[i])
}

dev.off()


