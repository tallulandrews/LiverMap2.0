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
obj <- readRDS("Hepatocyte2_varimax_Subcluster.rds")
obj_full <- obj
obj_all_genes <- readRDS("AllGenes/Hepatocyte2_harmony_Subcluster_Allgenes.rds")


cluster_col = "Coarse_clusters"


manual_cluster_anno <- c("PortalHep", "PortalHep","PortalHep","RiboHi", "CentralHep",
		"CentralHep", "Interzonal", "TSMB4X+", "CentralHep", "Erythroid", "Erythroblasts", 
		"Erythroid", "Female")


obj@meta.data$Subcluster_Manual <- manual_cluster_anno[obj@meta.data[,cluster_col]]

#saveRDS(obj@meta.data, "Hepatocyte2_fullmetadata.rds")


suppl_tab1 <- get_cluster_summary(obj, samples="sample")
suppl_tab2 <- get_cluster_summary(obj, samples="donor")
write.table(suppl_tab2, file="Hepatocyte2_cluster_summary.csv", sep=",")


### Object for Shiny ###
shiny_obj <- obj;
detect_rate <- group_rowmeans(shiny_obj@assays$RNA@counts > 0, shiny_obj@meta.data$Subcluster_Manual)
exclude_genes <- apply(detect_rate, 1, max) < 0.005
shiny_obj@meta.data <- shiny_obj@meta.data[,c("nCount_RNA","nFeature_RNA","percent.mt","Phase", "donor", "sample", "donor_sex", "donor_age","Subcluster_Manual")]
shiny_obj@reductions <- obj@reductions
shiny_obj@assays$RNA@scale.data <- matrix(0)
shiny_obj@assays$RNA@var.features <- c(0)
shiny_obj <- shiny_obj[!exclude_genes,]
expr_mat <- shiny_obj@assays$RNA@data
metadata <- shiny_obj@meta.data
metadata$UMAP1 <- shiny_obj@reductions$umap@cell.embeddings[,1]
metadata$UMAP2 <- shiny_obj@reductions$umap@cell.embeddings[,2]

saveRDS(shiny_obj, "ShinyApps/Hepatocyte2_Lattice_obj.rds")
saveRDS(list(expr_mat=expr_mat, metadata=metadata), "ShinyApps/Hepatocyte2_shiny.rds")





png("Hepato2_Sex_Proportions.png", width=6, height=8, units="in", res=300)
barplot_freq(obj, cluster_col="Subcluster_Manual", colours=colours_sex, metadata="donor_sex")
dev.off()
png("Hepato2_Age_Proportions.png", width=6, height=8, units="in", res=300)
barplot_freq(obj, cluster_col="Subcluster_Manual", colours=colours_age, metadata="donor_age_group")
dev.off()
png("Hepato2_Rejection_Proportions.png", width=6, height=8, units="in", res=300)
barplot_freq(obj, cluster_col="Subcluster_Manual", colours=colours_reject, metadata="trans.rejected")
dev.off()

png("Hepato2_donor_Frequency.png", width=9, height=6, units="in", res=300)
indivar_freq(obj@meta.data$Subcluster_Manual, obj@meta.data$donor, obj@meta.data)
dev.off()

indi_expr(obj, tag="Hepato2")


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

write.table(output, paste("Hepatocyte2_Nebula", assay_type, "DE.csv", sep="_"), sep=",", row.names=FALSE, col.names=TRUE, quote=TRUE)
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

saveRDS(list(sex=sex.output, age=age.output), "Hepatocyte2_MM_DE_5pr.rds")


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

fix_names <- c("Adaptive Immune", "TCR signaling", "Lymphoid-non-Lymphoid", "TCR", "CD28", "MHCII",
	"CD3", "IFNg", "PD-1", "ZAP-70 immuno synapse", "Chemokine", "", "2nd messenger", 
	"DAP12", "Interleukin10", "ROS/RNS in phagocytes", "Antigen presentation", "Nef",
	"", "RHO GTPases", "DAP12", "FCGR3A-pheocytosis", "", "Parasite", "Infection",
	"", "MHCI", "RHO Effectors", "", "lipids", "", "lipoprotein", "lipoprotein", "Retinoid", "vitamins", "",
	"IGF", "phototransduction", "phosphorylation", "Plasma Liproproetin")
G <- set.vertex.attribute(richments_de$graph, "name", value=fix_names)

png("Hepatocyte2_DE_Reactome_Sex_pathways.png", width=8, height=8, units="in", res=300)
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
		file="Hepatocyte2_sex_de.rds")


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

fix_names <- c("Porphyrins", "RNA", "tRNA", "rRNA", "Lymphoid-non-Lymphoid",
		"Complex I", "Citric acid cycle", "Electron trasnport chain", "", "rRNA", "tRNA")

G <- set.vertex.attribute(richments_de$graph, "name", value=fix_names)

png("Hepatocyte2_DE_Reactome_Rejection_pathways.png", width=8, height=8, units="in", res=300)
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
		file="Hepatocyte2_rejection_de.rds")



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

fix_names <- c("mt-rRNA", "mt-tRNA", "CD28", "CD3", "PD-1", "Lymphoid-non-Lymphoid",
	"ZAP70 immunosynapse", "TCR", "2nd messenger", "MHCII", "FCGR3A-phagocytosis",
	"", "Parasite", "IFNg", "actin-phagocytosis", "Chemokine", "MHCI",
	"TCR", "FCGR", "RHO GTPases", "vitamins", "phototransduction", "Retinoid",
	"", "Plasma lipoprotein")
G <- set.vertex.attribute(richments_de$graph, "name", value=fix_names)

png("Hepatocyte2_DE_Reactome_Age_pathways.png", width=8, height=8, units="in", res=300)
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
		file="Hepatocyte2_age_de.rds")



# ------------------- DE Vis ------------- #
source("C:/Users/tandrews/Documents/UHNSonya/Map2.2_Empty/Subcluster/IndividualVariability.R")
de_sex <- readRDS("Hepatocyte2_sex_de.rds")
de_age <- readRDS("Hepatocyte2_age_de.rds")
de_rej <- readRDS("Hepatocyte2_rejection_de.rds")

obj <- readRDS("Hepatocyte2_varimax_Subcluster.rds")

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

png(paste("Hepatocyte2", up, "Hi_genes.png", sep="_"), width=8, height=6, units="in", res=150)

par(mfrow=c(2,4))
for(g in top_genes) {
	indivar_DE_vis(obj@assays$RNA@scale.data[g,], obj@meta.data$donor, obj@meta.data, type=type)
#	indivar_DE_vis(obj@assays$RNA@data[g,], obj@meta.data$donor, obj@meta.data)

	title(main=g)
}
dev.off()


top_genes <- rownames(de2[order(de2[,5],decreasing=F),])[1:8]

png(paste("Hepatocyte2", dn, "Hi_genes.png", sep="_"), width=8, height=6, units="in", res=150)

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

png(paste("Hepatocyte2", up, "Hi_pathways.png", sep="_"), width=8, height=6, units="in", res=150)

par(mfrow=c(2,4))
for(i in 1:min(8, nrow(up_pathways))) {
	genes <- unlist(up_pathways[i, "leadingEdge"])
	indivar_DE_vis(obj@assays$RNA@scale.data[genes,], obj@meta.data$donor, obj@meta.data, type=type)
	title(main=up_pathways$pathway[i])
}

dev.off()

png(paste("Hepatocyte2", dn, "Hi_pathways.png", sep="_"), width=8, height=6, units="in", res=150)

par(mfrow=c(2,4))
for(i in 1:min(8, nrow(dn_pathways))) {
	genes <- unlist(dn_pathways[i, "leadingEdge"])
	indivar_DE_vis(obj@assays$RNA@scale.data[genes,], obj@meta.data$donor, obj@meta.data, type=type)
	title(main=dn_pathways$pathway[i])
}

dev.off()



