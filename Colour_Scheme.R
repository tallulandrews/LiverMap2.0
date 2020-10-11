# Function to get the same colours as used by default in Seurat DimPlot
get_seurat_colours <- function(obj, group.by) {
	require(scales)
	identities <- obj@meta.data[,group.by]
	if (class(identities) != "factor") {
		identities <- factor(identities)
	}
	identities <- levels(identities)

	my_color_palette <- hue_pal()(length(identities))
	return(my_color_palette)
#usage:
# TSNEPlot(object = object, do.return = T) + 
# scale_color_manual(values = my_color_palette)

}


## Map Subtypes to general cell-type labels ##
simplify_annotations <- function(annotations, types=c("B", "Mac", "T", "Hep")) {
	simplified <- as.character(annotations)
	if ("B" %in% types) {
	simplified[simplified %in% c(
		"AntibodysecretingBcells", 
		"MatureBcells", "B_cells", "any_B_cell", "Anti_B")] <- "Bcells"
	}
	if ("Mac" %in% types) {
	simplified[simplified %in% c(
		"InflamatoryMacrophages", "Non-inflammatoryMacrophages")] <- "Macrophages"
	}
	if ("T" %in% types) {
	simplified[simplified %in% c(
		"CD3abTcells", "gdTcells1", "gdTcells2", "gd_T_1", "gd_T_2", "CD3_T")] <- "Tcells"
	}
	if ("Hep" %in% types) {
	simplified[simplified %in% c(
		"PericentralHep", "UnidentifiedHep", "PeriportalHep",
		"interzonalHep", "Hep", "PortalHep2", "PortaHep1", "PortHep", "CentralHep1", "CentralHep")] <- "Hepatocyte"
	}
	return(simplified);
}

Cell_type_colours <- rbind(
	c("ambiguous", "#7F7F7F"), #grey50
	c("Eryth", "#8B0000"), #darkred
	c("Hepatocyte", "#e41a1c"), #bright red
	c("OtherHep", "#e41a1c"), #bright red
	c("PortalHep", "#084081"),
	c("interHep", "#0868ac"),
	c("CentralHep", "#2b8cbe"),
	c("Stellate", "#fec44f"),
	c("Cholangiocyte", "#ec7014"),
	c("Portalendo", "#993494"),
	c("Macrophage", "#377eb8"), # bright blue
	c("NonInfMac", "#377eb8"),
	c("InfMac", "#41b6c4"),
	c("Tcell", "#4daf4a"), #green
	c("gdTcells1", "#41ae76"),
	c("gdTcells2", "#00441b"),
	c("CD3abTcells", "#4daf4a"),
	c("Bcell", "#f781bf"), #violet
        c("AntiBcell", "#f78abf"), #violet
        c("MatBcell", "#6a51a3"), #violet
	c("NKcells", "#a65628"),
	c("LSECs", "#ffff33"),
	c("cvLSECs", "#ff7f00"),
	c("PortalLSECs", "#4d004b"))

colnames(Cell_type_colours) <- c("type", "colour")


map_cell_types <- function(types) {
	type <- as.character(types)
	tab <- unique(matrix(c("ambiguous", "ambiguous",
		"None", "ambiguous",
		"AntibodysecretingBcells", "AntiBcell",
		"AntiBcell", "AntiBcell",
		"Bcells", "Bcell",
		"B_cells", "Bcell",
		"any_B_cell", "Bcell",
		"Bcell", "Bcell",
		"CD3abTcells", "CD3abTcells",
		"CD3abTcell", "CD3abTcells",
		"NaiveTcell", "CD3abTcells",
		"CentralvenousLSECs", "cvLSECs",
		"CV_LSECs", "cvLSECs",
		"cvLSECs", "cvLSECs",
		"cvLSEC", "cvLSECs",
		"LSECs", "LSECs",
		"Cholangiocytes", "Cholangiocyte",
		"Cholangiocyte", "Cholangiocyte",
		"Cholan", "Cholangiocyte",
		"Chol", "Cholangiocyte",
		"Erythoidcells", "Eryth",
		"Erythoid", "Eryth",
		"Erythroid", "Eryth",
		"Erythroidcells", "Eryth",
		"Eryth", "Eryth",
		"gdTcells1", "gdTcells1",
		"gd_T_1", "gdTcells1",
		"gd_T_2", "gdTcells2",
		"gdTcell", "gdTcells2",
		"Hepatocyte", "Hepatocyte",
		"Hepatocytes", "Hepatocyte",
		"Hep", "Hepatocyte",
		"inflamatoryMacrophages", "InfMac",
		"InflamatoryMacrophages", "InfMac",
		"InflammatoryMacrophages", "InfMac",
		"inflammatoryMacrophages", "InfMac",
		"InfMac", "InfMac",
		"Infl_Mac", "InfMac",
		"interzonalHep", "interHep",
		"interHep", "interHep",
		"InterHep", "interHep",
		"Macrophage", "Macrophage",
		"Mac", "Macrophage",
		"Macrophages", "Macrophage",
		"MatureBcells", "MatBcell",
		"MatBcell", "MatBcell",
		"NK-likecells", "NKcells",
		"NK_like", "NKcells",
		"NKTcell", "NKcells",
		"NKTcells", "NKcells",
		"NKcells", "NKcells",
		"Non-inflammatoryMacrophages", "NonInfMac",
		"NonInfl_Mac", "NonInfMac",
		"NonInflMac1", "NonInfMac",
		"NonInflMac22", "NonInfMac",
		"NonInfMac", "NonInfMac",
		"PericentralHep", "CentralHep",
		"CentralHep", "CentralHep",
		"CentralHep1", "CentralHep",
		"PeriportalHep", "PortalHep",
		"PortaHep1", "PortalHep",
		"PortaHep2", "PortalHep",
		"PortalHep2", "PortalHep",
		"PortalHep1", "PortalHep",
		"PortalHep", "PortalHep",
		"PortHep", "PortalHep",
		"PeriportalLSECs", "PortalLSECs",
		"PeriLSECs", "PortalLSECs",
		"PortalLSECs", "PortalLSECs",
		"pLSECs", "PortalLSECs",
		"pLSEC", "PortalLSECs",
		"PortLSECs", "PortalLSECs",
		"Periportalendothelialcells", "Portalendo",
		"Portal_Endo", "Portalendo",
		"PortalEndo", "Portalendo",
		"Portal_Endo", "Portalendo",
		"PortEndo", "Portalendo",
		"Portalendothelialcells", "Portalendo",
		"Portalendo", "Portalendo",
		"Stellatecells", "Stellate",
		"Stellate", "Stellate",
		"Tcell", "Tcell",
		"Tcells", "Tcell",
		"gdTcells2", "gdTcells2",
		"UnidentifiedHep", "OtherHep",
		"OtherHep", "OtherHep",
		"Bcell", "Bcell",
		"Cholan", "Cholangiocyte",
		"PortHep", "PortalHep",
		"Anti_B", "AntiBcell",
		"any_B_cell", "Bcell",
		"any_T_cell", "Tcell",
		"B_cells", "Bcell",
		"CD3_T", "CD3abTcells",
		"Bcells", "Bcell"), ncol=2, byrow=TRUE))

	colnames(tab) <- c("orig", "tidy");
	return(tab[match(types,tab[,1]),2])
}
	


type_2_colour <- function(type) {
	type <- as.character(map_cell_types(type))
	colour <- Cell_type_colours[match(type, Cell_type_colours[,"type"]),"colour"]
	return(colour)
}


Type_DimPlot <- function(myseur, type_col="consistent_labs", reduction="umap", cluster_col="seurat_clusters") {
	require(ggplot2)
	agg_coord_by_cluster <- function(coords, clusters) {
		x <- split(seq(nrow(coords)), clusters)
		result <- sapply(x, function(a) apply(coords[a,],2,median))
		return(result)
	}

	umap_lab_pos <- agg_coord_by_cluster(Reductions(myseur, reduction)@cell.embeddings, myseur@meta.data[,cluster_col])

	# UMAP + Ref scmap anno
	new_colour_scheme <- Cell_type_colours[order(Cell_type_colours[,1]),]
	myseur@meta.data[,type_col] <- map_cell_types(myseur@meta.data[,type_col])
	new_colour_scheme <- new_colour_scheme[new_colour_scheme[,1] %in% myseur@meta.data[,type_col],]

	DimPlot(myseur, reduction=reduction, group.by=type_col, pt.size=.1)+scale_color_manual(values=new_colour_scheme[,2])+annotate("text", x=umap_lab_pos[1,], y=umap_lab_pos[2,], label=colnames(umap_lab_pos), colour="grey35")
}
