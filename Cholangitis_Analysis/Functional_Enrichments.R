#### ----- Set Up ----- ####

require("Seurat")
require("RColorBrewer")
require("ggplot2")
source("~/scripts/LiverMap2.0/Colour_Scheme.R")
source("~/scripts/LiverMap2.0/My_R_Scripts.R")

## ColourScheme ##
Healthy_cols = brewer.pal(9, "Blues")
PSC_cols = brewer.pal(9, "Reds")
PBC_cols = brewer.pal(9, "Purples")

my_plots <- function(seur_obj, outname, integration_name) {
	png(paste(outname, integration_name, "umap_sample.png", sep="_"),
                width=9, height =6, units="in", res=300)
                print(DimPlot(seur_obj, reduction="umap", pt.size=0.1,
                          group.by="orig.ident", label=TRUE))
        dev.off();
	png(paste(outname, integration_name, "umap_pheno.png", sep="_"),
                width=9, height =6, units="in", res=300)
                print(DimPlot(seur_obj, reduction="umap", pt.size=0.1,
                          group.by="Phenotype", label=TRUE))
        dev.off();
        png(paste(outname, integration_name, "barplot_sample.png", sep="_"),
                width=9, height =6, units="in", res=300)
                simpson_plot(seur_obj, samples="orig.ident", clusters="seurat_clusters")
        dev.off();
        png(paste(outname, integration_name, "barplot_pheno.png", sep="_"),
                width=9, height =6, units="in", res=300)
                simpson_plot(seur_obj, samples="Phenotype", clusters="seurat_clusters")
        dev.off();
        png(paste(outname, integration_name, "umap_autoanno.png", sep="_"),
                 width=9, height =6, units="in", res=300)
                 print(Type_DimPlot(seur_obj,type_col="marker_labs", cluster_col="marker_labs"))
        dev.off();
}

dir="/cluster/projects/macparland/TA/PostReview_Disease_vs_Healthy_Map"

args <- commandArgs(trailingOnly=TRUE)

##### Input #####
if (args[1] == "SC") {
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
} else if (args[1] == "SN") {

##### Input #####
	PSC = c(
        "/cluster/projects/macparland/TA/PSC/Processed/PSC005_Section5_TST_3pr_EmptyOnly.rds",
        "/cluster/projects/macparland/TA/PSC/Processed/PSC018_Section3_Nuclei_EmptyOnly.rds",
        "/cluster/projects/macparland/TA/PSC/Processed/PSC011_3pr_sn_TST_EmptyOnly.rds",
        "/cluster/projects/macparland/TA/PSC/Processed/PSC010_Section4_TST_3pr_EmptyOnly.rds",
        "/cluster/projects/macparland/TA/PSC/Processed/MacParland_SingleNuc_PSC014_Section3_3prV3_take2_EmptyOnly.rds",
        "/cluster/projects/macparland/TA/PSC/Processed/PSC019_Section6_TST_3pr_EmptyOnly.rds",
        "/cluster/projects/macparland/TA/PSC/Processed/PSC012_Section7_TST_3pr_EmptyOnly.rds"
        )

	PBC = c("/cluster/projects/macparland/TA/PBC/Processed/PBC001_Section1_TST_3pr_EmptyOnly.rds", 
	"/cluster/projects/macparland/TA/PBC/Processed/PBC004_SN_3pr_EmptyOnly.rds")

	Healthy = c("/cluster/projects/macparland/TA/LiverMap2.0/Cleaned/C41_TST_EmptyOnly.rds",
            "/cluster/projects/macparland/TA/LiverMap2.0/Cleaned/C58_TST_EmptyOnly.rds",
            "/cluster/projects/macparland/TA/LiverMap2.0/Cleaned/C70_TST_EmptyOnly.rds",
            "/cluster/projects/macparland/TA/LiverMap2.0/Cleaned/C72_TST_EmptyOnly.rds")

	all_samples <- c(PSC, PBC, Healthy)


	outname <- "SN_Integrated_Map";
} else {
	print("Error:Must provide either SC or SN as argument!")
	exit()
}

print(outname)

sample_names <- strsplit(all_samples, "/");
sample_names <- sapply(sample_names, function(x){x[length(x)]})
sample_names <- sub("_EmptyOnly.rds", "", sample_names)

###### Case vs Control ######

all_disease_de <- list()
files <- Sys.glob(paste(outname,"*edgeR_pseudobulkDE_rerun_all.csv", sep=""))
for (f in files) {
        tab <- read.table(f, sep=",")
        tmp <- unlist(strsplit(f, "_"))
        type <- tmp[length(tmp)-4]
        all_disease_de[[type]] <- tab;
}
outname2 = outname

# PBC vs Healthy
all_disease_de <- list()
files <- Sys.glob(paste(outname,"*edgeR_pseudobulkDE_PBC-Healthy_all.csv", sep=""))
for (f in files) {
        tab <- read.table(f, sep=",")
        tmp <- unlist(strsplit(f, "_"))
        type <- tmp[length(tmp)-4]
        all_disease_de[[type]] <- tab;
}
outname2 = paste(outname, "PBC-Healthy", sep="_")

#### Pathway Analysis ####

source("~/scripts/LiverMap2.0/My_R_Scripts.R")
require(fgsea)
require(gProfileR)

immune_path <- gmtPathways("/cluster/projects/macparland/TA/ExternalData/MSigDb_immune_signatures_c7.all.v7.1.symbols.gmt")
Hallmark_path <- gmtPathways("/cluster/projects/macparland/TA/ExternalData/MSigDbHalmarkPathways.gmt")
MSigAll <- gmtPathways("/cluster/projects/macparland/TA/ExternalData/MSigDb_curated_c2.all.v7.1.symbols.gmt")
reactome <- gmtPathways("/cluster/projects/macparland/TA/ExternalData/ReactomePathways.gmt")
BaderMSig <- gmtPathways("/cluster/projects/macparland/TA/ExternalData/BaderLab25Aug2020/Human_MSigdb_August_01_2020_symbol.gmt.txt")
BaderReact <- gmtPathways("/cluster/projects/macparland/TA/ExternalData/BaderLab25Aug2020/Human_Reactome_August_01_2020_symbol.gmt.txt")


cell_type_pathways <- read.table("/cluster/home/tandrews/scripts/PostReview_Disease_vs_Healthy_Map/Supplementary_Table_2_Marker_Gene_Lists.csv", sep=",", header=T)
cell_type_as_paths <- list()
for(general_type in unique(cell_type_pathways[,3])) {
	cell_type_as_paths[[general_type]] <- as.character(cell_type_pathways[cell_type_pathways[,3]==general_type,1])
}

extract_pathway_dat <- function(res, this_pathway, this_name) {
        genes <- unlist(res$rich[unlist(res$rich[,1]) == this_pathway,8])
        names(genes) <- rep(this_name, length(genes));

        bar_point <- unlist(abs(log10(res$rich[unlist(res$rich[,1]) == this_pathway,"padj"]))); names(bar_point)<- this_name;
	direction <- sign(paths$rich$NES[unlist(res$rich[,1]) == this_pathway])
        return(list(genes=genes, point=bar_point, dir=direction))
}

all_fsea_out <- list()
for (type in names(all_disease_de)) {
	set.seed(101);
	scores <- all_disease_de[[type]][,1]; names(scores) <- rownames(all_disease_de[[type]]);
	if (length(scores) == 1) {next;}
	paths <- do_fgsea(sort(scores), pathways=BaderMSig)
	if (!is.null(paths)) {
		all_fsea_out[[type]] <- paths
	}
}

saveRDS(all_fsea_out, file=paste(outname2, "fgsea_allDE_out.rds", sep="_"))

PSC_fsea_out <- all_fsea_out
PBC_fsea_out <- all_fsea_out
# check cell-type enriched paths

#all_fsea_out2 <- list()
#for (type in names(all_disease_de)) {
#        set.seed(101);
#        scores <- all_disease_de[[type]][,1]; names(scores) <- rownames(all_disease_de[[type]]);
#        paths <- do_fgsea(sort(scores), pathways=cell_type_as_paths)
#        if (!is.null(paths)) {
#                all_fsea_out2[[type]] <- paths
#        }
#}

saveRDS(all_fsea_out2, file=paste(outname2, "fgsea_celltype_out.rds", sep="_"))



##### ---- Individual Pathway Plots ---- #####
require(pheatmap)
require(RColorBrewer)
## SC
outname <- "SC_Integrated_Map";
outname2 = paste(outname, "PBC-Healthy", sep="_")
#outname2=outname

all_paths <- unlist(sapply(all_fsea_out, function(x){unlist(x$rich$pathway)}))
sort(table(all_paths))

all_fsea_out=readRDS(paste(outname2, "fgsea_allDE_out.rds", sep="_"))
cell_types <- names(all_fsea_out); cell_types <- cell_types[!grepl("Doublet", cell_types)]

pathways <- c("HALLMARK_COMPLEMENT%MSIGDB_C2%HALLMARK_COMPLEMENT",
	"HALLMARK_OXIDATIVE_PHOSPHORYLATION%MSIGDB_C2%HALLMARK_OXIDATIVE_PHOSPHORYLATION",
	"HALLMARK_TNFA_SIGNALING_VIA_NFKB%MSIGDB_C2%HALLMARK_TNFA_SIGNALING_VIA_NFKB",
	"HALLMARK_XENOBIOTIC_METABOLISM%MSIGDB_C2%HALLMARK_XENOBIOTIC_METABOLISM",
	"ST_T_CELL_SIGNAL_TRANSDUCTION%MSIGDB_C2%ST_T_CELL_SIGNAL_TRANSDUCTION",
	"PID_IL12_2PATHWAY%MSIGDB_C2%PID_IL12_2PATHWAY",
	"HALLMARK_INTERFERON_GAMMA_RESPONSE%MSIGDB_C2%HALLMARK_INTERFERON_GAMMA_RESPONSE",
	"HALLMARK_INTERFERON_ALPHA_RESPONSE%MSIGDB_C2%HALLMARK_INTERFERON_ALPHA_RESPONSE",
	"HALLMARK_IL2_STAT5_SIGNALING%MSIGDB_C2%HALLMARK_IL2_STAT5_SIGNALING",
	"PID_CXCR4_PATHWAY%MSIGDB_C2%PID_CXCR4_PATHWAY",
	"HALLMARK_BILE_ACID_METABOLISM%MSIGDB_C2%HALLMARK_BILE_ACID_METABOLISM",
	"PID_IL2_1PATHWAY%MSIGDB_C2%PID_IL2_1PATHWAY",
	"HALLMARK_FATTY_ACID_METABOLISM%MSIGDB_C2%HALLMARK_FATTY_ACID_METABOLISM",
	"HALLMARK_ALLOGRAFT_REJECTION%MSIGDB_C2%HALLMARK_ALLOGRAFT_REJECTION",
	"HALLMARK_KRAS_SIGNALING_UP%MSIGDB_C2%HALLMARK_KRAS_SIGNALING_UP",
	"HALLMARK_IL6_JAK_STAT3_SIGNALING%MSIGDB_C2%HALLMARK_IL6_JAK_STAT3_SIGNALING",
	"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION%MSIGDB_C2%HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
	"HALLMARK_P53_PATHWAY%MSIGDB_C2%HALLMARK_P53_PATHWAY",
	"HALLMARK_APOPTOSIS%MSIGDB_C2%HALLMARK_APOPTOSIS",
	"HALLMARK_MITOTIC_SPINDLE%MSIGDB_C2%HALLMARK_MITOTIC_SPINDLE",
	"HALLMARK_MYOGENESIS%MSIGDB_C2%HALLMARK_MYOGENESIS",
	"HALLMARK_INFLAMMATORY_RESPONSE%MSIGDB_C2%HALLMARK_INFLAMMATORY_RESPONSE",
	"HALLMARK_HEME_METABOLISM%MSIGDB_C2%HALLMARK_HEME_METABOLISM",
	"NABA_MATRISOME%MSIGDB_C2%NABA_MATRISOME", 

	"PID_HNF3B_PATHWAY%MSIGDB_C2%PID_HNF3B_PATHWAY", 
	"BIOCARTA_CTLA4_PATHWAY%MSIGDB_C2%BIOCARTA_CTLA4_PATHWAY", 
	"PID_IL8_CXCR2_PATHWAY%MSIGDB_C2%PID_IL8_CXCR2_PATHWAY", 
	"BIOCARTA_TH1TH2_PATHWAY%MSIGDB_C2%BIOCARTA_TH1TH2_PATHWAY", 
	"PID_IL27_PATHWAY%MSIGDB_C2%PID_IL27_PATHWAY", 
	"PID_IL23_PATHWAY%MSIGDB_C2%PID_IL23_PATHWAY", 
	"PID_IL4_2PATHWAY%MSIGDB_C2%PID_IL4_2PATHWAY", 
	"BIOCARTA_IL10_PATHWAY%MSIGDB_C2%BIOCARTA_IL10_PATHWAY", 
	"BIOCARTA_CELLCYCLE_PATHWAY%MSIGDB_C2%BIOCARTA_CELLCYCLE_PATHWAY", 
	"PID_NFAT_TFPATHWAY%MSIGDB_C2%PID_NFAT_TFPATHWAY" 
	
	)

names(pathways) <- c("Complement cascade", "Oxidative phosphorylation", "TNFa signaling", "Xenobiotic Metabolism", 
			"T cell signal transduction", "IL12 pathway", "IFNg response", "IFNa response",
			"IL2_STAT5 signaling", "CXCR4 pathway", "Bile Acid Metabolism", "IL2 pathway", 
			"Fatty acid metabolism", "Allograft rejection", "KRAS signaling", "IL6_JAK_STAT3 signaling",
			"Epithelial Mesenchymal Transition", "P53 pathway", "Apoptosis", "Mitotic spindle", "Myogenesis",
			"Inflammatory response", "Heme Metabolism", "Matrisome", 

			"HNF3B pathway", "CTLA4 pathway", "IL8-CXCR2 pathway", "Th1Th2 pathway", "IL27 pathway", "IL23 pathway", 
			"IL4 pathway", "IL10 pathway", "Cell Cycle", "NFAT pathway"

			)


out_mat <- matrix(NA, nrow=length(cell_types), ncol=length(pathways))
for (i in 1:length(cell_types)) {
        this_type <- cell_types[i]
        rich <- all_fsea_out[[this_type]]$rich
        out_mat[i,] <- unlist(rich[match(pathways, rich$pathway),"NES"])
}

rownames(out_mat) <- cell_types
colnames(out_mat) <- names(pathways)

my_heatmap_colours <- colorRampPalette(c(brewer.pal(n=9, "RdPu")[6:9], "black", rev(brewer.pal(n = 7, name = "RdYlBu")[1:4])))(200)

out_mat[is.na(out_mat)] <- 0
require(pheatmap)
legend_thing <- seq(from=floor(min(out_mat)), to=floor(max(out_mat)), length=5)
png(paste(outname2, "DEall_pathway_enrichment_heatmap.png", sep="_"), width=7, height=6, units="in", res=300)
breaks = seq(from=-1*(max(abs(out_mat))), to=max(abs(out_mat)), length=length(my_heatmap_colours)+1)
pheatmap(t(out_mat), legend_breaks = c(legend_thing, max(out_mat)), main = "", legend_labels = c(legend_thing, "NES\n"), scale="none",
		col=my_heatmap_colours, breaks=breaks)
dev.off()
write.table(out_mat, file=paste(outname2, "DEall_pathway_enrichment_heatmap.csv", sep="_"), sep=",", row.names=T, col.names=T)

out_mat_PSC <- out_mat
out_mat_PSC[is.na(out_mat_PSC)] <- 0

common_types <- intersect(rownames(out_mat_PSC), rownames(out_mat_PBC))
out_mat_PSC[match(common_types, rownames(out_mat_PSC)),] - out_mat_PBC[match(common_types, rownames(out_mat_PBC)),]

# Cleaned:
PBC_matched <- out_mat_PBC[match(rownames(out_mat_PSC), rownames(out_mat_PBC)),];
PBC_matched[is.na(PBC_matched)] <- 0
PSC_filtered <- out_mat_PSC
PSC_filtered[sign(PSC_filtered) == sign(PBC_matched)] <- 0

PBC_filtered <- PBC_matched
PBC_filtered[sign(PBC_filtered) == sign(out_mat_PSC)] <- 0

png(paste("SC_PSC_DEall_filtered_pathway_enrichment_heatmap.png", sep="_"), width=7, height=6, units="in", res=300)
legend_thing <- seq(from=floor(min(out_mat_PSC)), to=floor(max(out_mat_PSC)), length=5)
breaks = seq(from=-1*(max(abs(out_mat_PSC))), to=max(abs(out_mat_PSC)), length=length(my_heatmap_colours)+1)
pheatmap(t(PSC_filtered), legend_breaks = c(legend_thing, max(out_mat_PSC)), main = "", legend_labels = c(legend_thing, "NES\n"), scale="none",
                col=my_heatmap_colours, breaks=breaks)
dev.off()

png(paste("SC_PBC_DEall_filtered_pathway_enrichment_heatmap.png", sep="_"), width=7, height=6, units="in", res=300)
legend_thing <- seq(from=floor(min(out_mat_PSC)), to=floor(max(out_mat_PSC)), length=5)
breaks = seq(from=-1*(max(abs(out_mat_PSC))), to=max(abs(out_mat_PSC)), length=length(my_heatmap_colours)+1)
pheatmap(t(PBC_filtered), legend_breaks = c(legend_thing, max(out_mat_PSC)), main = "", legend_labels = c(legend_thing, "NES\n"), scale="none",
                col=my_heatmap_colours, breaks=breaks)
dev.off()


## cNK cell details ##

require(pheatmap)
require(RColorBrewer)
## SC
outname <- "SC_Integrated_Map";
outname2 = paste(outname, "PBC-Healthy", sep="_")
#outname2=outname

PSC_paths_fsea_out=readRDS(paste(outname, "fgsea_allDE_out.rds", sep="_"))
PBC_paths_fsea_out=readRDS(paste(outname2, "fgsea_allDE_out.rds", sep="_"))

PSC_cNK <- PSC_paths_fsea_out[["cNK"]]
PBC_cNK <- PBC_paths_fsea_out[["cNK"]]


all_disease_de <- list()
files <- Sys.glob(paste(outname,"*edgeR_pseudobulkDE_rerun_all.csv", sep=""))
for (f in files) {
        tab <- read.table(f, sep=",")
        tmp <- unlist(strsplit(f, "_"))
        type <- tmp[length(tmp)-4]
        all_disease_de[[type]] <- tab;
}
all_disease_de_PSC <- all_disease_de

all_disease_de <- list()
files <- Sys.glob(paste(outname,"*edgeR_pseudobulkDE_PBC-Healthy_all.csv", sep=""))
for (f in files) {
        tab <- read.table(f, sep=",")
        tmp <- unlist(strsplit(f, "_"))
        type <- tmp[length(tmp)-4]
        all_disease_de[[type]] <- tab;
}
all_disease_de_PBC <- all_disease_de

PSC_cNK_de <- all_disease_de_PSC[["cNK"]]
PBC_cNK_de <- all_disease_de_PBC[["cNK"]]

# NOTHING CAME OF THIS, I don't know why PBC doesn't have significant upregulation of inflammation pathways in cNKs because they upregulate the same genes????


# PBC vs Healthy
all_disease_de <- list()
files <- Sys.glob(paste(outname,"*edgeR_pseudobulkDE_PBC-Healthy_all.csv", sep=""))
for (f in files) {
        tab <- read.table(f, sep=",")
        tmp <- unlist(strsplit(f, "_"))
        type <- tmp[length(tmp)-4]
        all_disease_de[[type]] <- tab;
}
outname2 = paste(outname, "PBC-Healthy", sep="_")

g =c("ETV5", "IFNG", "CD28", "CD3G", "FOS", "CD86", "JUN", "IL18RAP", "IL2", "TGFB1", "CD3D", "STAT4", "PRF1", "CD80", "IL18", "PPP3CA", "TBX21", "CD247", "IRF1", "HLA-DRB1", "PPP3R1")

PSC_cNK_de[g,]
PBC_cNK_de[g,]

####### ---------------- SN --------------- ##########
outname <- "SN_Integrated_Map";
outname2 = "SN_Integrated_Map_PBC-Healthy"

all_fsea_out=readRDS(paste(outname2, "fgsea_allDE_out.rds", sep="_"))

all_paths <- unlist(sapply(all_fsea_out, function(x){unlist(x$rich$pathway)}))
sort(table(all_paths))

cell_types <- names(all_fsea_out); cell_types <- cell_types[!grepl("Doublet", cell_types)]
pathways <- c("HALLMARK_BILE_ACID_METABOLISM%MSIGDB_C2%HALLMARK_BILE_ACID_METABOLISM",
		"HALLMARK_HYPOXIA%MSIGDB_C2%HALLMARK_HYPOXIA",
		"HALLMARK_TNFA_SIGNALING_VIA_NFKB%MSIGDB_C2%HALLMARK_TNFA_SIGNALING_VIA_NFKB",
		"HALLMARK_COMPLEMENT%MSIGDB_C2%HALLMARK_COMPLEMENT",
		"NABA_MATRISOME%MSIGDB_C2%NABA_MATRISOME",
		"HALLMARK_INTERFERON_GAMMA_RESPONSE%MSIGDB_C2%HALLMARK_INTERFERON_GAMMA_RESPONSE",
		"HALLMARK_IL2_STAT5_SIGNALING%MSIGDB_C2%HALLMARK_IL2_STAT5_SIGNALING",
		"PID_CXCR4_PATHWAY%MSIGDB_C2%PID_CXCR4_PATHWAY",
		"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION%MSIGDB_C2%HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
		"HALLMARK_APOPTOSIS%MSIGDB_C2%HALLMARK_APOPTOSIS",
		"HALLMARK_MYOGENESIS%MSIGDB_C2%HALLMARK_MYOGENESIS",
		"HALLMARK_INFLAMMATORY_RESPONSE%MSIGDB_C2%HALLMARK_INFLAMMATORY_RESPONSE",
		"HALLMARK_FATTY_ACID_METABOLISM%MSIGDB_C2%HALLMARK_FATTY_ACID_METABOLISM",
		"HALLMARK_KRAS_SIGNALING_UP%MSIGDB_C2%HALLMARK_KRAS_SIGNALING_UP",

		"HALLMARK_ALLOGRAFT_REJECTION%MSIGDB_C2%HALLMARK_ALLOGRAFT_REJECTION",
		"HALLMARK_APICAL_JUNCTION%MSIGDB_C2%HALLMARK_APICAL_JUNCTION",
		"HALLMARK_IL6_JAK_STAT3_SIGNALING%MSIGDB_C2%HALLMARK_IL6_JAK_STAT3_SIGNALING",
		"WNT_SIGNALING%MSIGDB_C2%WNT_SIGNALING",
		"PID_IL4_2PATHWAY%MSIGDB_C2%PID_IL4_2PATHWAY"
		)
names(pathways) <- c("Bile Acid Metabolism", "Hypoxia", "TNFa signaling", "Complement Cascade", "Matrisome", "IFNg response",
			"IL2_STAT5 Signaling","CXCR4 pathway", "Epithelial Mesenchymal Transistion", "Apoptosis",
			"Myogenesis", "Infammatory Response", "Fatty Acid Metabolism",
			"KRAS signaling",

			"Allograft rejection", "Apical junction", "IL6_JAK_STAT3 signaling", "Wnt signaling", "IL4 pathway")

out_mat <- matrix(NA, nrow=length(cell_types), ncol=length(pathways))
for (i in 1:length(cell_types)) {
	this_type <- cell_types[i]
	rich <- all_fsea_out[[this_type]]$rich
	out_mat[i,] <- unlist(rich[match(pathways, rich$pathway),"NES"])
}

rownames(out_mat) <- cell_types
colnames(out_mat) <- names(pathways)

my_heatmap_colours <- colorRampPalette(c(brewer.pal(n=9, "RdPu")[6:9], "black", rev(brewer.pal(n = 7, name = "RdYlBu")[1:4])))(200)

out_mat[is.na(out_mat)] <- 0
require(pheatmap)
legend_thing <- seq(from=floor(min(out_mat)), to=floor(max(out_mat)), length=5)
png(paste(outname2, "allDE_pathway_enrichment_heatmap.png", sep="_"), width=6, height=6, units="in", res=300)
breaks = seq(from=-1*(max(abs(out_mat))), to=max(abs(out_mat)), length=length(my_heatmap_colours)+1)
pheatmap(t(out_mat), legend_breaks = c(legend_thing, max(out_mat)), main = "", legend_labels = c(legend_thing, "NES\n"),
		col=my_heatmap_colours, breaks=breaks)
dev.off()
write.table(out_mat, file=paste(outname2, "DEall_pathway_enrichment_heatmap.csv", sep="_"), sep=",", row.names=T, col.names=T)

out_mat_PSC <- out_mat
out_mat_PBC <- out_mat
	
out_mat_PSC[is.na(out_mat_PSC)] <- 0
out_mat_PBC[is.na(out_mat_PBC)] <- 0

common_types <- intersect(rownames(out_mat_PSC), rownames(out_mat_PBC))
out_mat_PSC[match(common_types, rownames(out_mat_PSC)),] - out_mat_PBC[match(common_types, rownames(out_mat_PBC)),]

# Cleaned:
PBC_matched <- out_mat_PBC[match(rownames(out_mat_PSC), rownames(out_mat_PBC)),];
PBC_matched[is.na(PBC_matched)] <- 0
rownames(PBC_matched) <- rownames(out_mat_PSC)
PSC_filtered <- out_mat_PSC
PSC_filtered[sign(PSC_filtered) == sign(PBC_matched)] <- 0

PBC_filtered <- PBC_matched
PBC_filtered[sign(PBC_filtered) == sign(out_mat_PSC)] <- 0

png(paste("SN_PSC_DEall_filtered_pathway_enrichment_heatmap.png", sep="_"), width=6, height=6, units="in", res=300)
legend_thing <- seq(from=floor(min(out_mat_PSC)), to=floor(max(out_mat_PSC)), length=5)
breaks = seq(from=-1*(max(abs(out_mat_PSC))), to=max(abs(out_mat_PSC)), length=length(my_heatmap_colours)+1)
pheatmap(t(PSC_filtered), legend_breaks = c(legend_thing, max(out_mat_PSC)), main = "", legend_labels = c(legend_thing, "NES\n"), scale="none",
                col=my_heatmap_colours, breaks=breaks)
dev.off()

png(paste("SN_PBC_DEall_filtered_pathway_enrichment_heatmap.png", sep="_"), width=6, height=6, units="in", res=300)
legend_thing <- seq(from=floor(min(out_mat_PSC)), to=floor(max(out_mat_PSC)), length=5)
breaks = seq(from=-1*(max(abs(out_mat_PSC))), to=max(abs(out_mat_PSC)), length=length(my_heatmap_colours)+1)
pheatmap(t(PBC_filtered), legend_breaks = c(legend_thing, max(out_mat_PSC)), main = "", legend_labels = c(legend_thing, "NES\n"), scale="none",
                col=my_heatmap_colours, breaks=breaks)
dev.off()





########################################################################################################


names(V(paths$graph))

this_names <- rev(c("IFNa Response", "IFNg Response", "Bile Metabolism", "EMT", "Xenobiotic Metabolism", "Coagulation", "ECM", "Matrisome Associated", "Matrisome")) # lrNK
this_names <- rev(c("IL12", "CD8 TCR", "IFNa Response", "Allograft Rejection", "TCR pathway", "Downstream CD8", "IFNg Response", "CXCR4 Pathway", "MYC targets", "Adipogenesis", "Myogenesis", "EMT", "Estrogen response", "Fatty Acid Metabolism", "Bile Metabolism", "Matrisome Associated", "Matrisome", "ECM", "Coagulation", "Xenobiotic Metabolism" )) # cNK
this_names <- c("Matrisome", "Coagulation", "Matrisome Associated", "ECM", "Xenobiotic Metablism", "Estrogen Response", "Bile Acid Metablism", "IFNg Response", "Allograft Rejection") # CD8+T
this_names <- c("ECM", "Coagulation", "Xenobiotic Metabolism","Matrisome associated", "Matrisome") # CD4+T
this_names <- c("Xenobiotic Metabolism", "ECM Regulators", "Coagulation", "Matrisome") # Monocyte
this_names <- c("Inflammatory Response", "Matrisome Associated", "Allograft Rejection", "Matrisome", "NABA Secreted Factors") # C-Hepato
this_names <- c("Matrisome", "Matrisome Associated", "IFNg Response", "Allograft Rejection", "IFNa Response") # P-Hepato
this_names <- c("Matrisome", "Matrisome Associated") # Kupffer
this_names <- c("ECM", "Coagulation", "Xenobiotic Metabolism", "Matrisome associated", "Matrisome", "TNFA Signalling") # cvLSEC

bars <- c();
for(i in 1:length(V(paths$graph))) {
	p = names(V(paths$graph))[i]
	n = this_names[i]
	summary <- extract_pathway_dat(paths, p, n)
	bars <- c(bars, summary$point*summary$dir)
}
png(paste(outname, type, "pathways.png", sep="_"), width=8, height=8, units="in", res=150)
V(paths$graph)$label <- this_names
plot(paths$graph, vertex.color=paths$vertex_col)
dev.off()
png(paste(outname, type, "pathways2.png", sep="_"), width=8, height=8, units="in", res=150)
par(mar=c(4,10,1,1))
barplot(bars, names=this_names, las=2, horiz=TRUE, xlab="log10(p-value)", col=paths$vertex_col)
dev.off()




#Kupffer - downregulation of Complement Cascade
#lrNK - up regulation of KLRC1/2/D1 - inhibitory molecules & HLA-I


tmp<- all_disease_de[["Kupffer"]]
sort(rownames(tmp)[tmp[,1] < 0])

### Cell-Coll Communication ###
# Version Conflicts
#library(CellChat)
#library(patchwork)
#options(stringsAsFactors = FALSE)

#objs <- readRDS(paste(outname, "Objects_forCellCellInterations.rds", sep="_"))

#healthy_cellchat <- createCellChat(object=objs$healthy_5pr, group.by="full_annotation")
#healthy_cellchat@DB <- CellChatDB.human
#healthy_cellchat <- subsetData(healthy_cellchat)
#psc_cellchat <- createCellChat(object=objs$psc, group.by="cell_type")
#psc_cellchat@DB <- CellChatDB.human
#psc_cellchat <- subsetData(psc_cellchat)

#cell_chat_pipeline <- function(cellchat_obj) {
#	cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
#	cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
#	cellchat_obj <- computeCommunProb(cellchat_obj, type="truncatedMean", trim=0.05)
#	cellchat_obj <- computeCommunProbPathway(cellchat_obj, type="truncatedMean", trim=0.05)
#	cellchat_obj <- aggregateNet(cellchat_obj)
#	return(cellchat_obj)
#}

#healthy_cellchat <- cell_chat_pipeline(healthy_cellchat)
#psc_cellchat <- cell_chat_pipeline(psc_cellchat)






